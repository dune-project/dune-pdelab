// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include <sstream>
#include <vector>
#include <map>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/pdelab/finiteelementmap/p0fem.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/common/vtkexport.hh>
#include <dune/pdelab/localoperator/laplacedirichletccfv.hh>
#include <dune/pdelab/backend/backendselector.hh>
#include <dune/pdelab/backend/istlvectorbackend.hh>
#include <dune/pdelab/backend/istl/bcrsmatrixbackend.hh>
#include <dune/pdelab/backend/seqistlsolverbackend.hh>

#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/ordering/singlecodimleafordering.hh>
#include <dune/pdelab/gridfunctionspace/vtk.hh>


#include"gridexamples.hh"


// define some grid functions to interpolate from
template<typename GV, typename RF>
class G
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
                                                  G<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,G<GV,RF> > BaseT;

  G (const GV& gv) : BaseT(gv) {}
  inline void evaluateGlobal (const typename Traits::DomainType& x,
							  typename Traits::RangeType& y) const
  {
    typename Traits::DomainType center;
    for (int i=0; i<GV::dimension; i++) center[i] = 0.5;
    center -= x;
	y = exp(-center.two_norm2());
  }
};

// define some boundary grid functions to define boundary conditions
template<typename GV>
class B
  : public Dune::PDELab::BoundaryGridFunctionBase<Dune::PDELab::BoundaryGridFunctionTraits<GV,int,1,
                                                                                           Dune::FieldVector<int,1> >,
                                                  B<GV> >
{
  GV gv;

public:
  typedef Dune::PDELab::BoundaryGridFunctionTraits<GV,int,1,Dune::FieldVector<int,1> > Traits;
  typedef Dune::PDELab::BoundaryGridFunctionBase<Traits,B<GV> > BaseT;

  B (const GV& gv_) : gv(gv_) {}

  template<typename I>
  inline void evaluate (const Dune::PDELab::IntersectionGeometry<I>& ig,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    y = 1; // all is Dirichlet boundary
  }

  //! get a reference to the GridView
  inline const GV& getGridView () const
  {
    return gv;
  }
};


template<class GV>
void test (const GV& gv)
{
  typedef typename GV::Grid::ctype DF;
  typedef double RF;
  const int dim = GV::dimension;

  // instantiate finite element maps
  Dune::GeometryType gt;
  gt.makeCube(dim);
  typedef Dune::PDELab::P0LocalFiniteElementMap<DF,RF,dim> FEM;
  FEM fem(gt); // works only for cubes

#ifdef TEST_SIMPLIFIED_INFRASTRUCTURE
  typedef Dune::PDELab::SingleCodimMapper Mapper;
#else
  typedef Dune::PDELab::DefaultLeafOrderingTag Mapper;
#endif

  // make function space
  typedef Dune::PDELab::GridFunctionSpace<
    GV,
    FEM,
    Dune::PDELab::NoConstraints,
    Dune::PDELab::ISTLVectorBackend<>,
    Mapper> GFS;
  GFS gfs(gv,fem);
  gfs.name("u");
  gfs.setDataSetType(GFS::Output::cellData);

  typedef G<GV,RF> GType;
  GType g(gv);

  // make grid function operator
  typedef Dune::PDELab::LaplaceDirichletCCFV<GType> LO;
  LO lo(g);

  typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
  MBE mbe(GV::dimension == 2 ? 5 : 7); // P0 on quadrilaterals

  typedef Dune::PDELab::GridOperator<
    GFS,GFS,LO,
    MBE,
    RF,RF,RF> GO;
  GO go(gfs,gfs,lo,mbe);

  // make coefficent Vector and initialize it from a function
  typedef typename GO::Traits::Domain V;
  V x0(gfs);
  x0 = 0.0;

  Dune::PDELab::interpolate(g,gfs,x0);

  // represent operator as a matrix
  typedef typename GO::Traits::Jacobian M;
  M m(go);
  m = 0.0;

  go.jacobian(x0,m);

  // evaluate residual w.r.t initial guess
  V r(gfs);
  r = 0.0;
  go.residual(x0,r);

  typedef typename M::Container ISTLM;
  typedef typename V::Container ISTLV;

  // make ISTL solver
  Dune::MatrixAdapter<ISTLM,ISTLV,ISTLV> opa(m.base());
  //  typedef Dune::PDELab::OnTheFlyOperator<V,V,GOS> ISTLOnTheFlyOperator;
  //  ISTLOnTheFlyOperator opb(gos);
  //  Dune::SeqSSOR<M,V,V> ssor(m,1,1.0);
  Dune::SeqILU0<ISTLM,ISTLV,ISTLV> ilu0(m.base(),1.0);
  //  Dune::Richardson<V,V> richardson(1.0);
  Dune::CGSolver<ISTLV> solvera(opa,ilu0,1E-10,5000,2);
  //  Dune::CGSolver<V> solverb(opb,richardson,1E-10,5000,2);
  Dune::InverseOperatorResult stat;

  // solve the jacobian system
  r *= -1.0; // need -residual
  V x(gfs,0.0);
  solvera.apply(x.base(),r.base(),stat);
  x += x0;

  // output grid function with VTKWriter
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::nonconforming);
  Dune::PDELab::addSolutionToVTKWriter(vtkwriter,gfs,x);

  std::stringstream vtu_name;

#ifdef TEST_SIMPLIFIED_INFRASTRUCTURE
  vtu_name << "testlaplacedirichletccfv-simplified";
#else
  vtu_name << "testlaplacedirichletccfv";
#endif

  vtu_name << "-" << dim << "D";

  vtkwriter.write(vtu_name.str(),Dune::VTK::ascii);
}

int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper::instance(argc, argv);

    // 2D
    {
      std::cout << "2D tests" << std::endl;
      // need a grid in order to test grid functions
      Dune::FieldVector<double,2> L(1.0);
      Dune::array<int,2> N(Dune::fill_array<int,2>(1));
      std::bitset<2> B(false);
      Dune::YaspGrid<2> grid(L,N,B,0);
      grid.globalRefine(3);

      test(grid.leafGridView());
    }

    // 3D
    {
      std::cout << "3D tests" << std::endl;
      // need a grid in order to test grid functions
      Dune::FieldVector<double,3> L(1.0);
      Dune::array<int,3> N(Dune::fill_array<int,3>(1));
      std::bitset<3> B(false);
      Dune::YaspGrid<3> grid(L,N,B,0);
      grid.globalRefine(3);

      test(grid.leafGridView());
    }

	// test passed
	return 0;

  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
	return 1;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
	return 1;
  }
}
