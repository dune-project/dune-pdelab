// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cmath>
#include <cstddef>
#include <deque>
#include <iostream>
#include <string>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/geometrytype.hh>
#include <dune/common/misc.hh>
#include <dune/common/mpihelper.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/common/tuples.hh>

#include <dune/grid/common/quadraturerules.hh>
#include <dune/grid/io/file/vtk/common.hh>
#include <dune/grid/io/file/vtk/vtksequencewriter.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/pdelab/backend/istlmatrixbackend.hh>
#include <dune/pdelab/backend/istlsolverbackend.hh>
#include <dune/pdelab/backend/istlvectorbackend.hh>
#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/common/geometrywrapper.hh>
#include <dune/pdelab/common/vtkexport.hh>
#include <dune/pdelab/finiteelementmap/conformingconstraints.hh>
#include <dune/pdelab/finiteelementmap/p1fem.hh>
#include <dune/pdelab/function/const.hh>
#include <dune/pdelab/gridfunctionspace/constraints.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/gridoperatorspace/localmatrix.hh>
#include <dune/pdelab/localoperator/defaultimp.hh>
#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/localoperator/idefault.hh>
#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/localoperator/poisson.hh>
#include <dune/pdelab/localoperator/scaled.hh>
#include <dune/pdelab/multistep/gridoperatorspace.hh>
#include <dune/pdelab/multistep/method.hh>
#include <dune/pdelab/multistep/parameter.hh>
#include <dune/pdelab/newton/newton.hh>

//===============================================================
//===============================================================
// Solve 1D homogenous wave equation
//   \partial_t^2 u - c \partial_x^2 u = 0 in \Omega,
//                                   u = 0 on \partial\Omega
//===============================================================
//===============================================================

//===============================================================
// Define parameter functions f,g,j and \partial\Omega_D/N
//===============================================================

////////////////////////////////////////////////////////////////////////
//
// boundary grid function selecting boundary conditions
//

template<typename GV>
class B
  : public Dune::PDELab::BoundaryGridFunctionBase
          < Dune::PDELab::
            BoundaryGridFunctionTraits<GV,int,1,
                                       Dune::FieldVector<int,1> >,
            B<GV> >,
    public Dune::PDELab::InstationaryFunctionDefaults
{
  const GV& gv;

public:
  typedef Dune::PDELab::BoundaryGridFunctionTraits
    <GV,int,1,Dune::FieldVector<int,1> > Traits;
  typedef Dune::PDELab::BoundaryGridFunctionBase<Traits,B<GV> > BaseT;

  B (const GV& gv_) : gv(gv_) {}

  template<typename I>
  inline void evaluate (const Dune::PDELab::IntersectionGeometry<I>& ig,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    y = 1; // Dirichlet
  }

  //! get a reference to the GridView
  inline const GV& getGridView () const
  {
    return gv;
  }
};

////////////////////////////////////////////////////////////////////////
//
// function for Dirichlet boundary conditions and initialization
//

template<typename GV, typename RF, typename Time>
class G
  : public Dune::PDELab::AnalyticGridFunctionBase
           < Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
             G<GV,RF,Time> >
{
  Time time;

public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,G<GV,RF,Time> > BaseT;

  G (const GV& gv) : BaseT(gv), time(0) {}
  void setTime(Time time_) { time = time_; }

  inline void evaluateGlobal (const typename Traits::DomainType& x,
                              typename Traits::RangeType& y) const
  {
    typename Traits::DomainType center;
    for (int i=0; i<GV::dimension; i++) center[i] = 0.5;
    center[0] += time;
    center -= x;
    y = std::exp(-center.two_norm2()/(2*Dune::SQR(0.05)));
  }
};

//===============================================================
// Local Operators
//===============================================================

template<typename Time>
class R1
  : public Dune::PDELab::LocalOperatorDefaultFlags
  , public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<Time>
{ };

template<typename Time>
class R2
  : public Dune::PDELab::FullVolumePattern
  , public Dune::PDELab::LocalOperatorDefaultFlags
  , public Dune::PDELab::JacobianBasedAlphaVolume<R2<Time> >
  , public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<Time>
{
  unsigned qorder;

public:
  // pattern assembly flags
  enum { doPatternVolume = true };
  enum { doAlphaVolume = true };

  R2(int qorder_ = 2) : qorder(qorder_) { }

  template<typename EG, typename LFSU, typename X, typename LFSV,
           typename R>
  void jacobian_volume(const EG& eg,
                       const LFSU& lfsu, const X& x, const LFSV& lfsv,
                       Dune::PDELab::LocalMatrix<R>& mat) const
  {
    // domain and range field type
    typedef typename LFSU::Traits::LocalFiniteElementType LFEU;
    typedef typename LFEU::Traits::LocalBasisType LBU;
    typedef typename LBU::Traits::RangeType RangeU;

    typedef typename LFSV::Traits::LocalFiniteElementType LFEV;
    typedef typename LFEV::Traits::LocalBasisType LBV;
    typedef typename LBV::Traits::RangeType RangeV;

    typedef typename LBU::Traits::DomainFieldType DF;
    static const unsigned dimD = LBU::Traits::dimDomain;
    typedef typename LBU::Traits::DomainType Domain;

    // select quadrature rule
    typedef Dune::QuadratureRule<DF,dimD> QR;
    typedef Dune::QuadratureRules<DF,dimD> QRs;
    Dune::GeometryType gt = eg.geometry().type();
    const QR& rule = QRs::rule(gt,qorder);

    // loop over quadrature points
    for(typename QR::const_iterator it=rule.begin();
        it!=rule.end(); ++it) {
      std::vector<RangeU> phiu(lfsu.size());
      lfsu.localFiniteElement().localBasis()
        .evaluateFunction(it->position(),phiu);

      std::vector<RangeV> phiv(lfsv.size());
      lfsv.localFiniteElement().localBasis()
        .evaluateFunction(it->position(),phiv);

      R factor = it->weight()
        * eg.geometry().integrationElement(it->position());

      for(unsigned i = 0; i < lfsu.size(); ++i)
        for(unsigned j = 0; j < lfsv.size(); ++j)
          mat(i,j) += factor * (phiu[i] * phiv[j]);
    }
  }
};

//===============================================================
// Problem setup and solution
//===============================================================

// generate a P1 function and output it
template<typename GV, typename FEM, typename CON, int qorder>
void wave (const GV& gv, const FEM& fem, typename GV::ctype dt,
           typename GV::ctype c, unsigned steps, std::string filename)
{
  // constants and types
  typedef typename GV::Grid::ctype DF;
  typedef typename FEM::Traits::LocalFiniteElementType::Traits::
    LocalBasisType::Traits::RangeFieldType RF;

  // make function space
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,
    Dune::PDELab::ISTLVectorBackend<1> > GFS;
  GFS gfs(gv,fem);

  // make constraints map and initialize it from a function
  typedef typename GFS::template ConstraintsContainer<RF>::Type C;
  C cg;
  cg.clear();
  typedef B<GV> BType;
  BType b(gv);
  Dune::PDELab::constraints(b,gfs,cg);

  // make coefficent Vector and initialize it from a function
  typedef typename GFS::template VectorContainer<RF>::Type V;
  std::deque<Dune::shared_ptr<V> > oldvalues;
  typedef G<GV,RF,DF> GType;
  GType g(gv);
  oldvalues.push_front(Dune::shared_ptr<V>(new V(gfs)));
  g.setTime(-dt);
  Dune::PDELab::interpolate(g,gfs,*oldvalues.front());
  oldvalues.push_front(Dune::shared_ptr<V>(new V(gfs)));
  g.setTime(0);
  Dune::PDELab::interpolate(g,gfs,*oldvalues.front());

  // make parameters
  Dune::PDELab::CentralDifferencesParameters<DF> msParams;

  // make grid function operator
  typedef Dune::PDELab::ConstGridFunction<GV, RF, 1> ZeroFunc;
  ZeroFunc zeroFunc(gv, 0);

  typedef Dune::PDELab::InstationaryPoisson<DF, ZeroFunc, BType, ZeroFunc,
    qorder> Poisson;
  Poisson poisson(zeroFunc, b, zeroFunc);

  typedef Dune::PDELab::ScaledLocalOperator<Poisson, RF, DF> R0;
  R0 r0(poisson, c*c);
  R1<DF> r1;
  R2<DF> r2(qorder);

  typedef Dune::tuple<R0, R1<DF>, R2<DF> > LOPs;
  typedef Dune::tuple<R0&, R1<DF>&, R2<DF>&> LOPRefs;

  typedef Dune::PDELab::MultiStepGridOperatorSpace<DF,V,GFS,GFS,
    LOPs,C,C,Dune::PDELab::ISTLBCRSMatrixBackend<1,1> > MGOS;
  MGOS mgos(msParams, gfs,cg,gfs,cg, LOPRefs(r0, r1, r2));

// #if HAVE_SUPERLU
//   typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LS;
//   LS ls(false);
// #else
  typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
  LS ls(5000,false);
  //#endif

  // <<<7>>> make Newton for time-dependent problem
  typedef Dune::PDELab::Newton<MGOS,LS,V> PDESOLVER;
  PDESOLVER tnewton(mgos,ls);
  tnewton.setReassembleThreshold(0.0);
  tnewton.setVerbosityLevel(0);
  tnewton.setReduction(0.9);
  tnewton.setMinLinearReduction(1e-9);

  // <<<8>>> time-stepper
  Dune::PDELab::MultiStepMethod<DF,MGOS,PDESOLVER,V,V>
    msMethod(msParams,mgos,tnewton);
  msMethod.setVerbosityLevel(2);

  // output grid function with VTKWriter
  Dune::VTKSequenceWriter<GV> vtkwriter(gv,filename,"","");

  // make discrete function object
  typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
  {
    DGF dgf(gfs,*oldvalues[1]);
    vtkwriter.addVertexData
      (new Dune::PDELab::VTKGridFunctionAdapter<DGF>(dgf,"solution"));
    vtkwriter.write(-dt,Dune::VTKOptions::binary);
    vtkwriter.clear();
  }
  {
    DGF dgf(gfs,*oldvalues[0]);
    vtkwriter.addVertexData
      (new Dune::PDELab::VTKGridFunctionAdapter<DGF>(dgf,"solution"));
    vtkwriter.write(0,Dune::VTKOptions::binary);
    vtkwriter.clear();
  }

  DF time = 0;

  for(unsigned step = 0; step < steps; ++step) {
    Dune::shared_ptr<V> xnew(new V(*oldvalues.front()));

    dt = msMethod.apply(time, dt, oldvalues, *xnew);
    time += dt;

    oldvalues.pop_back();
    oldvalues.push_front(xnew);

    // output grid function with VTKWriter
    DGF dgf(gfs,*xnew);
    vtkwriter.addVertexData
      (new Dune::PDELab::VTKGridFunctionAdapter<DGF>(dgf,"solution"));
    vtkwriter.write(time,Dune::VTKOptions::binary);
    vtkwriter.clear();
  }
}

//===============================================================
// Main program with grid setup
//===============================================================

int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper::instance(argc, argv);

    // YaspGrid P1 1D test
    {
      typedef Dune::YaspGrid<1> Grid;
      // make grid
      Dune::FieldVector<double,1> L(1.0);
      Dune::FieldVector<int,1> N(1);
      Dune::FieldVector<bool,1> B(false);
      Grid grid(L,N,B,0);
      grid.globalRefine(7);

      // get view
      typedef Grid::LeafGridView GV;
      const GV& gv=grid.leafView();

      // make finite element map
      typedef GV::Grid::ctype DF;
      typedef Dune::PDELab::P1LocalFiniteElementMap<DF,double, 1> FEM;
      FEM fem;

      // solve problem
      wave<GV,FEM,Dune::PDELab::ConformingDirichletConstraints,2>
        (gv,fem,1.0/(1<<9), 1.0, 1000 ,"testmultistep_yasp_P1_1d");
    }

    // test passed
    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
    throw;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
    throw;
  }
}
