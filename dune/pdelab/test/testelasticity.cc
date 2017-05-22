// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <iostream>
#include <vector>
#include <map>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/common/vtkexport.hh>
#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/localoperator/linearelasticity.hh>

#include <dune/pdelab/gridfunctionspace/vectorgridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/vtk.hh>

static const double szX = 10.0;

// define some grid functions to interpolate from
template<typename GV, typename RF>
class G
    : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,GV::dimension>,
                                                  G<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,GV::dimension> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,G<GV,RF> > BaseT;

  G (const GV& gv) : BaseT(gv) {}
  inline void evaluateGlobal (const typename Traits::DomainType& x,
                              typename Traits::RangeType& y) const
  {
    y = 0.0;
    if (x[0] == szX)
    {
        y[0] = 0.1;
    }
  }
};

// define some boundary grid functions to define boundary conditions
template<typename GV>
class ModelProblem
  : public Dune::PDELab::LinearElasticityParameterInterface<
  Dune::PDELab::LinearElasticityParameterTraits<GV, double>,
  ModelProblem<GV> >
{
public:

  typedef Dune::PDELab::LinearElasticityParameterTraits<GV, double> Traits;

  ModelProblem(typename Traits::RangeType G,
    typename Traits::RangeFieldType l,
    typename Traits::RangeFieldType m) :
    G_(G), lambda_(l), mu_(m)
  {}

  void
  f (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
    typename Traits::RangeType & y) const
  {
    y = G_;
  }

  template<typename I>
  bool isDirichlet(const I & ig,
    const typename Traits::IntersectionDomainType & coord
    ) const
  {
    typename Traits::DomainType xg = ig.geometry().global( coord );
    return (xg[0] == 0.0); // || xg[0] == szX);  // Dirichlet b.c. on left & right boundary
  }

  void
  u (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
    typename Traits::RangeType & y) const
  {
    y = 0.0;
  }

  typename Traits::RangeFieldType
  lambda (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return lambda_;
  }

  typename Traits::RangeFieldType
  mu (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return mu_;
  }

private:

  typename Traits::RangeType G_;
  typename Traits::RangeFieldType lambda_;
  typename Traits::RangeFieldType mu_;

};

// generate a P1 function and output it
template<class GV>
void testp1 (const GV& gv, double mu, double lambda, double constG)
{
  using Dune::PDELab::Backend::native;

  typedef typename GV::Grid::ctype DF;

  const int dim = GV::dimension;

  // instantiate finite element maps
  typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,double,1> FEM;
  FEM fem(gv);

  // make function space
  typedef Dune::PDELab::ConformingDirichletConstraints Constraints;
  typedef Dune::PDELab::ISTL::VectorBackend<> ComponentVectorBackend;

  typedef Dune::PDELab::DefaultLeafOrderingTag Mapper;

  typedef Dune::PDELab::LexicographicOrderingTag OrderingTag;
  typedef Dune::PDELab::VectorGridFunctionSpace<
    GV,
    FEM,
    dim,
    Dune::PDELab::ISTL::VectorBackend<>,
    ComponentVectorBackend,
    Constraints,
    OrderingTag,
    Mapper
    > GFS;
  GFS gfs(gv,fem);
  gfs.name("displacement");

  // model description
  typedef ModelProblem<GV> Param;
  Dune::FieldVector<double, dim> G(0.0); G[dim-1] = -constG;
  Param param(G, mu, lambda);

  // make constraints map and initialize it from a function
  typedef typename GFS::template ConstraintsContainer<double>::Type C;
  C cg;
  cg.clear();
  Dune::PDELab::constraints(param,gfs,cg);

  std::cout << gfs.size() << " DOFs\n";
  std::cout << cg.size() << " constraint DOFs\n";

  // make local operator
  typedef Dune::PDELab::LinearElasticity<Param> LO;
  LO lo(param);

  typedef Dune::PDELab::ISTL::BCRSMatrixBackend<> MBE;
  MBE mbe(9); // 2D Q1

  // make grid operator
  typedef Dune::PDELab::GridOperator<
    GFS,GFS,LO,
    MBE,
    double,double,double,
    C,C> GOS;
  GOS gos(gfs,cg,gfs,cg,lo,mbe);

  // make coefficent Vector and initialize it from a function
  typedef typename GOS::Traits::Domain V;
  V x0(gfs);
  x0 = 0.0;
  typedef Dune::PDELab::LinearElasticityDirichletExtensionAdapter<Param> Displacement;
  Displacement u_fnkt(gv, param);
  Dune::PDELab::interpolate(u_fnkt,gfs,x0);
  Dune::PDELab::set_nonconstrained_dofs(cg,0.0,x0);

  // represent operator as a matrix
  typedef typename GOS::Traits::Jacobian M;
  M m(gos);
  m = 0.0;
  gos.jacobian(x0,m);
  if (gfs.size() <= 16)
    Dune::printmatrix(std::cout,native(m),"global stiffness matrix","row",9,3);

  // evaluate residual w.r.t initial guess
  V r(gfs);
  r = 0.0;
  gos.residual(x0,r);

  // make ISTL solver
  typedef typename M::Container ISTL_M;
  typedef typename V::Container ISTL_V;
  Dune::MatrixAdapter<ISTL_M,ISTL_V,ISTL_V> opa(native(m));
  Dune::SeqILU0<ISTL_M,ISTL_V,ISTL_V> ilu0(native(m),1e-2);

  Dune::CGSolver<ISTL_V> solver(opa,ilu0,1E-20,5000,2);
  Dune::InverseOperatorResult stat;

  // solve the jacobian system
  r *= -1.0; // need -residual
  V x(gfs,0.0);
  Dune::PDELab::set_nonconstrained_dofs(cg,1.0,x);
  solver.apply(native(x),native(r),stat);
  x += x0;

  // output grid function with VTKWriter
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
  Dune::PDELab::addSolutionToVTKWriter(vtkwriter,gfs,x);
  vtkwriter.write("testelasticity",Dune::VTK::ascii);
}

int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper::instance(argc, argv);

    double mu = 100.0;
    double lambda = 10000.0;
    double g = 1.0;
    int level=2;

    if (argc > 1)
        mu = atof(argv[1]);
    if (argc > 2)
        lambda = atof(argv[2]);
    if (argc > 3)
        g = atof(argv[3]);
    if (argc > 4)
        level = atoi(argv[4]);

    std::cout << "mu     = " << mu << "\n"
              << "lambda = " << lambda << "\n"
              << "g = " << g << std::endl;

    Dune::FieldVector<double,2> L(1); L[0] = szX;
    std::array<int,2> N(Dune::fill_array<int,2>(1)); N[0] = szX;
    Dune::YaspGrid<2> grid(L,N);
    grid.globalRefine(level);

    testp1(grid.leafGridView(), mu, lambda, g);

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
