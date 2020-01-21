// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include <vector>
#include <map>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/pdelab.hh>

#include "gridexamples.hh"


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

// generate a P1 function and output it
template<class GV>
void testp1 (const GV& gv)
{
  typedef typename GV::Grid::ctype DF;

  // instantiate finite element maps
  typedef Dune::PDELab::PkLocalFiniteElementMap<GV,DF,double,1> FEM;
  FEM fem(gv);

  // make function space
  typedef Dune::PDELab::GridFunctionSpace<
    GV,
    FEM,
    Dune::PDELab::ConformingDirichletConstraints,
    Dune::PDELab::ISTLVectorBackend<1>
    > GFS;
  GFS gfs(gv,fem);

  // make constraints map and initialize it from a function
  typedef typename GFS::template ConstraintsContainer<double>::Type C;
  C cg;
  cg.clear();
  typedef B<GV> BType;
  BType b(gv);
  Dune::PDELab::constraints(b,gfs,cg);

  // make coefficent Vector and initialize it from a function
  using V = Dune::PDELab::Backend::Vector<GFS,double>;
  V x0(gfs);
  x0 = 0.0;
  typedef G<GV,double> GType;
  GType g(gv);
  Dune::PDELab::interpolate(g,gfs,x0);
  Dune::PDELab::set_nonconstrained_dofs(cg,0.0,x0);

  // make grid function operator
  Dune::PDELab::LaplaceDirichletP12D la;
  typedef Dune::PDELab::GridOperator<
    GFS,GFS,
    Dune::PDELab::LaplaceDirichletP12D,C,C,Dune::PDELab::ISTLBCRSMatrixBackend<1,1> > GOS;
  GOS gos(gfs,cg,gfs,cg,la);

  // represent operator as a matrix
  typedef typename GOS::template MatrixContainer<double>::Type M;
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
  V x(gfs,0.0);
  solvera.apply(x,r,stat);
  x += x0;

  // make discrete function object
  typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
  DGF dgf(gfs,x);

  // output grid function with VTKWriter
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
  vtkwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF> >(dgf,"p1"));
  vtkwriter.write("testlaplacedirichletp12d",Dune::VTK::ascii);
}

int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper::instance(argc, argv);

#if HAVE_UG
    std::shared_ptr<Dune::UGGrid<2> > uggrid(TriangulatedUnitSquareMaker<Dune::UGGrid<2> >::create());
  	uggrid->globalRefine(3);
    testp1(uggrid->leafGridView());
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
