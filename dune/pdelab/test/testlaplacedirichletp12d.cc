// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifdef HAVE_CONFIG_H
#include "config.h"     
#endif
#include<iostream>
#include<vector>
#include<map>
#include<dune/common/mpihelper.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/static_assert.hh>
#include<dune/grid/yaspgrid.hh>
#include<dune/istl/bvector.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/preconditioners.hh>
#include<dune/istl/io.hh>
//#include<dune/istl/paamg/amg.hh>

#include"../finiteelementmap/p0fem.hh"
#include"../finiteelementmap/p12dfem.hh"
#include"../finiteelementmap/pk2dfem.hh"
#include"../gridfunctionspace/gridfunctionspace.hh"
#include"../gridfunctionspace/gridfunctionspaceutilities.hh"
#include"../gridfunctionspace/interpolate.hh"
#include"../gridfunctionspace/constraints.hh"
#include"../common/function.hh"
#include"../common/vtkexport.hh"
#include"../gridoperatorspace/gridoperatorspace.hh"
#include"../localoperator/laplacedirichletp12d.hh"
#include"../backend/istlvectorbackend.hh"
#include"../backend/istlmatrixbackend.hh"
#include"../backend/istlsolverbackend.hh"

#include"gridexamples.hh"


// define some grid functions to interpolate from
template<typename GV, typename RF>
class F
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
                                                  F<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,F<GV,RF> > BaseT;

  F (const GV& gv) : BaseT(gv) {}
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
  const GV& gv;

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
  inline const GV& getGridView ()
  {
    return gv;
  }
};

// generate a P1 function and output it
template<class GV> 
void testp1 (const GV& gv)
{
  typedef typename GV::Grid::ctype DF;
  const int dim = GV::dimension;

  // instantiate finite element maps
  typedef Dune::PDELab::P12DLocalFiniteElementMap<DF,double> P1FEM;
  P1FEM p1fem;
  
  // make function space
  typedef Dune::PDELab::GridFunctionSpace<GV,P1FEM,
    Dune::PDELab::P12DConstraints,Dune::PDELab::ISTLVectorBackend<1> > P1GFS; 
  P1GFS p1gfs(gv,p1fem);

  // make constraints map and initialize it from a function
  typedef typename P1GFS::template ConstraintsContainer<double>::Type P1C;
  P1C p1cg;
  p1cg.clear();
  typedef B<GV> BType;
  BType b(gv);
  Dune::PDELab::constraints(b,p1gfs,p1cg);

  // make coefficent Vector and initialize it from a function
  typedef typename P1GFS::template VectorContainer<double>::Type P1V;
  P1V x0(p1gfs);
  x0 = 0.0;
  typedef F<GV,double> FType;
  FType f(gv);
  Dune::PDELab::interpolate(f,p1gfs,x0);
  Dune::PDELab::set_nonconstrained_dofs(p1cg,0.0,x0);

  // make grid function operator
  Dune::PDELab::LaplaceDirichletP12D la;
  typedef Dune::PDELab::GridOperatorSpace<P1GFS,P1GFS,
    Dune::PDELab::LaplaceDirichletP12D,P1C,P1C,Dune::PDELab::ISTLBCRSMatrixBackend<1,1> > P1OP;
  P1OP p1op(p1gfs,p1cg,p1gfs,p1cg,la);

  // represent operator as a matrix
  typedef typename P1OP::template MatrixContainer<double>::Type P1M;
  P1M p1m(p1op);
  p1m = 0.0;
  p1op.jacobian(x0,p1m);
  //  Dune::printmatrix(std::cout,p1m.base(),"global stiffness matrix","row",9,1);

  // evaluate residual w.r.t initial guess
  P1V r(p1gfs);
  r = 0.0;
  p1op.residual(x0,r);

  // make ISTL solver
  Dune::MatrixAdapter<P1M,P1V,P1V> opa(p1m);
  typedef Dune::PDELab::OnTheFlyOperator<P1V,P1V,P1OP> ISTLOnTheFlyOperator;
  ISTLOnTheFlyOperator opb(p1op);
  Dune::SeqSSOR<P1M,P1V,P1V> ssor(p1m,1,1.0);
  Dune::SeqILU0<P1M,P1V,P1V> ilu0(p1m,1.0);
  Dune::Richardson<P1V,P1V> richardson(1.0);

//   typedef Dune::Amg::CoarsenCriterion<Dune::Amg::SymmetricCriterion<P1M,
//     Dune::Amg::FirstDiagonal> > Criterion;
//   typedef Dune::SeqSSOR<P1M,P1V,P1V> Smoother;
//   typedef typename Dune::Amg::SmootherTraits<Smoother>::Arguments SmootherArgs;
//   SmootherArgs smootherArgs;
//   smootherArgs.iterations = 2;
//   int maxlevel = 20, coarsenTarget = 100;
//   Criterion criterion(maxlevel, coarsenTarget);
//   criterion.setMaxDistance(2);
//   typedef Dune::Amg::AMG<Dune::MatrixAdapter<P1M,P1V,P1V>,P1V,Smoother> AMG;
//   AMG amg(opa,criterion,smootherArgs,1,1);

  Dune::CGSolver<P1V> solvera(opa,ilu0,1E-10,5000,2);
  Dune::CGSolver<P1V> solverb(opb,richardson,1E-10,5000,2);
  Dune::InverseOperatorResult stat;

  // solve the jacobian system
  r *= -1.0; // need -residual
  P1V x(p1gfs,0.0);
  solvera.apply(x,r,stat);
  x += x0;

  // make discrete function object
  typedef Dune::PDELab::DiscreteGridFunction<P1GFS,P1V> P1DGF;
  P1DGF p1dgf(p1gfs,x);
  
  // output grid function with VTKWriter
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<P1DGF>(p1dgf,"p1"));
  vtkwriter.write("testlaplacedirichletp12d",Dune::VTKOptions::ascii);
}

int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

#if HAVE_UG
 	UGUnitSquare uggrid;
  	uggrid.globalRefine(4);
    testp1(uggrid.leafView());
#endif

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
