// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifdef HAVE_CONFIG_H
#include "config.h"     
#endif

#include <iostream>
#include <string>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/smartpointer.hh>

#include <dune/grid/albertagrid/agrid.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include <dune/istl/operators.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>

#include "../backend/istlmatrixbackend.hh"
#include "../backend/istlsolverbackend.hh"
#include "../backend/istlvectorbackend.hh"
#include "../common/function.hh"
#include "../common/geometrywrapper.hh"
#include "../common/vtkexport.hh"
#include "../finiteelementmap/conformingconstraints.hh"
#include "../finiteelementmap/edges03dfem.hh"
#include "../gridfunctionspace/constraints.hh"
#include "../gridfunctionspace/gridfunctionspace.hh"
#include "../gridfunctionspace/gridfunctionspaceutilities.hh"
#include "../gridfunctionspace/interpolate.hh"
#include "../gridoperatorspace/gridoperatorspace.hh"
#include "../localoperator/electrostatic.hh"

#include"stddomains.hh"


//===============================================================
//===============================================================
// Solve the Electrostatic "wave" equation
//           rot(1/mu * rot E) = 0 in \Omega, 
//                 n x (n x E) = g on \partial\Omega_D
//===============================================================
//===============================================================

//===============================================================
// Define parameter function mu
//===============================================================

// function for defining the source term
template<typename GV, typename RF>
class Mu
  : public Dune::PDELab::AnalyticGridFunctionBase<
      Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
      Mu<GV,RF>
    >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,Mu<GV,RF> > BaseT;

  Mu (const GV& gv) : BaseT(gv) {}

  inline void
  evaluateGlobal (const typename Traits::DomainType& x, 
                  typename Traits::RangeType& y) const
  {
    y = 1.0;
  }
};

// boundary grid function selecting boundary conditions 
template<typename GV>
class B
  : public Dune::PDELab::BoundaryGridFunctionBase<
      Dune::PDELab::BoundaryGridFunctionTraits<
        GV,
        int,1,Dune::FieldVector<int,1>
      >,
      B<GV>
    >
{
  const GV& gv;

public:
  typedef Dune::PDELab::BoundaryGridFunctionTraits<GV,int,1,Dune::FieldVector<int,1> > Traits;
  typedef Dune::PDELab::BoundaryGridFunctionBase<Traits,B<GV> > BaseT;

  B (const GV& gv_) : gv(gv_) {}

  template<typename I>
  inline void
  evaluate (const Dune::PDELab::IntersectionGeometry<I>& ig, 
            const typename Traits::DomainType& x,
            typename Traits::RangeType& y) const
  {
    y = 1; // Dirichlet
  }

  //! get a reference to the GridView
  inline const GV& getGridView ()
  {
    return gv;
  }
};

// function for Dirichlet boundary conditions and initialization
template<typename GV, typename RF>
class G
  : public Dune::PDELab::AnalyticGridFunctionBase<
      Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,3>,
      G<GV,RF>
    >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,3> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,G<GV,RF> > BaseT;

  G (const GV& gv)
    : BaseT(gv)
  {
    prescribedE[0] = 1;
    prescribedE[1] = 0;
    prescribedE[2] = 0;
  }

  inline void
  evaluateGlobal (const typename Traits::DomainType& x, 
                  typename Traits::RangeType& y) const
  {
    y = prescribedE;
  }

private:
  typename Traits::RangeType prescribedE;
};

//===============================================================
// Problem setup and solution 
//===============================================================

// generate a P1 function and output it
template<typename GV, typename FEM, typename CON, int q> 
void electrostatic (const GV& gv, const FEM& fem, std::string filename)
{
  // constants and types
  typedef typename GV::Grid::ctype DF;
  typedef typename FEM::Traits::LocalFiniteElementType::Traits::
    LocalBasisType::Traits::RangeFieldType R;

  // make function space
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,
    Dune::PDELab::ISTLVectorBackend<1> > GFS; 
  GFS gfs(gv,fem);

  // make constraints map and initialize it from a function
  typedef typename GFS::template ConstraintsContainer<R>::Type C;
  C cg;
  cg.clear();
  typedef B<GV> BType;
  BType b(gv);
  Dune::PDELab::constraints(b,gfs,cg);

  // make coefficent Vector and initialize it from a function
  typedef typename GFS::template VectorContainer<R>::Type V;
  V x0(gfs);
  x0 = 0.0;
  typedef G<GV,R> GType;
  GType g(gv);
  Dune::PDELab::interpolateGlobal(g,gfs,x0);
  Dune::PDELab::set_nonconstrained_dofs(cg,0.0,x0);

#if 0
  // make grid function operator
  typedef Mu<GV,R> MuType;
  MuType mu(gv);
  typedef Dune::PDELab::Electrostatic<MuType,q> LOP; 
  LOP lop(mu);
  typedef Dune::PDELab::GridOperatorSpace<GFS,GFS,
    LOP,C,C,Dune::PDELab::ISTLBCRSMatrixBackend<1,1> > GOS;
  GOS gos(gfs,cg,gfs,cg,lop);

  // represent operator as a matrix
  typedef typename GOS::template MatrixContainer<R>::Type M;
  M m(gos);
  m = 0.0;
  gos.jacobian(x0,m);
  //  Dune::printmatrix(std::cout,m.base(),"global stiffness matrix","row",9,1);

  // evaluate residual w.r.t initial guess
  V r(gfs);
  r = 0.0;
  gos.residual(x0,r);

  // make ISTL solver
  Dune::MatrixAdapter<M,V,V> opa(m);
  typedef Dune::PDELab::OnTheFlyOperator<V,V,GOS> ISTLOnTheFlyOperator;
  ISTLOnTheFlyOperator opb(gos);
  Dune::SeqSSOR<M,V,V> ssor(m,1,1.0);
  Dune::SeqILU0<M,V,V> ilu0(m,1.0);
  Dune::Richardson<V,V> richardson(1.0);

//   typedef Dune::Amg::CoarsenCriterion<Dune::Amg::SymmetricCriterion<M,
//     Dune::Amg::FirstDiagonal> > Criterion;
//   typedef Dune::SeqSSOR<M,V,V> Smoother;
//   typedef typename Dune::Amg::SmootherTraits<Smoother>::Arguments SmootherArgs;
//   SmootherArgs smootherArgs;
//   smootherArgs.iterations = 2;
//   int maxlevel = 20, coarsenTarget = 100;
//   Criterion criterion(maxlevel, coarsenTarget);
//   criterion.setMaxDistance(2);
//   typedef Dune::Amg::AMG<Dune::MatrixAdapter<M,V,V>,V,Smoother> AMG;
//   AMG amg(opa,criterion,smootherArgs,1,1);

  Dune::CGSolver<V> solvera(opa,ilu0,1E-10,5000,2);
  Dune::CGSolver<V> solverb(opb,richardson,1E-10,5000,2);
  Dune::InverseOperatorResult stat;

  // solve the jacobian system
  r *= -1.0; // need -residual
#endif
  V x(gfs,0.0);
//  solvera.apply(x,r,stat);
  x += x0;

  // make discrete function object
  typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
  DGF dgf(gfs,x);
  
  // output grid function with VTKWriter
  Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,1);
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(dgf,"solution"));
  vtkwriter.write(filename,Dune::VTKOptions::ascii);
}

//===============================================================
// Main program with grid setup
//===============================================================

int main(int argc, char** argv)
{
  try{
    // 77 is special and means "test was skipped".  Return that if non of the
    // supported grids were available
    int result = 77;

#ifdef HAVE_ALBERTA
    {
      typedef Dune::AlbertaGrid<3,3> Grid;
      Dune::SmartPointer<Grid> grid = Dune::PDELab::UnitTetrahedronMaker<Grid>::create();
      //grid->globalRefine(1);

      // get view
      typedef Grid::LeafGridView GV;
      const GV& gv=grid->leafView(); 

      // make finite element map
      typedef Grid::ctype DF;
      typedef Dune::PDELab::EdgeS03DLocalFiniteElementMap<GV,double> FEM;
      FEM fem(gv);
  
      // solve problem
      electrostatic
        <GV,FEM,Dune::PDELab::ConformingDirichletConstraints,2>
        (gv,fem,"electrostatic_alberta");
    }
    result = 0;
#endif // HAVE_ALBERTA

	return result;
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
