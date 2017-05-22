// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <iostream>
#include <vector>
#include <map>
#include <string>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/memory/blocked_allocator.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/gridfunctionspace/vtk.hh>
#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/pdelab/finiteelementmap/pkfem.hh>
#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/localoperator/convectiondiffusionfem.hh>

#include "gridexamples.hh"

//===============================================================
//===============================================================
// Solve the Poisson equation
//           - \Delta u = f in \Omega,
//                    u = g on \partial\Omega_D
//  -\nabla u \cdot \nu = j on \partial\Omega_N
//===============================================================
//===============================================================

//===============================================================
// Define parameter functions f,g,j and \partial\Omega_D/N
//===============================================================

// function for defining the source term

template<typename GV, typename RF>
struct Problem
  : public Dune::PDELab::ConvectionDiffusionModelProblem<GV,RF>
{

  using Base = Dune::PDELab::ConvectionDiffusionModelProblem<GV,RF>;

public:

  using Traits = typename Base::Traits;

  template<typename E, typename X>
  auto f(const E& e, const X& xl) const
  {
    auto x = e.geometry().global(xl);
    return (x[0]>0.25 && x[0]<0.375 && x[1]>0.25 && x[1]<0.375) ? 50.0 : 0.0;
  }

  template<typename I, typename X>
  auto bctype(const I& i, const X& xl) const
  {
    auto x = i.geometry().global(xl);
    if ((x[1]<1E-6 || x[1]>1.0-1E-6) || (x[0]>1.0-1E-6 && x[1]>0.5+1E-6))
      {
        return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
      }
    return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
  }

  template<typename E, typename X>
  auto g(const E& e, const X& xl) const
  {
    auto x = e.geometry().global(xl);
    auto center = x;
    center = 0.5;
    center -= x;
    using std::exp;
    return exp(-center.two_norm2());
  }

  template<typename I, typename X>
  auto j(const I& i, const X& xl) const
  {
    auto x = i.geometry().global(xl);
    if (x[1]<1E-6 || x[1]>1.0-1E-6)
      return 0.0;
    else
      return -5.0;
  }

};


//===============================================================
// Problem setup and solution
//===============================================================

// generate a P1 function and output it
template<typename GV, typename FEM, typename CON, int q>
void poisson (const GV& gv, const FEM& fem, std::string filename, int chunk_size)
{
  // constants and types
  typedef typename FEM::Traits::FiniteElementType::Traits::
    LocalBasisType::Traits::RangeFieldType R;

  using Dune::PDELab::Backend::native;

  typedef Dune::Memory::blocked_cache_aligned_allocator<R,std::size_t,8> Allocator;


  typedef Dune::PDELab::istl::FlatVectorBackend<Allocator> FVB;
  typedef Dune::PDELab::istl::FlatMatrixBackend<Allocator> FMB;

  // make function space
  typedef Dune::PDELab::GridFunctionSpace<
    GV,
    FEM,
    CON,
    FVB
    > GFS;
  GFS gfs(gv,fem);
  gfs.name("solution");

  // make local operator
  using Problem = ::Problem<GV,R>;
  Problem problem;

  using LOP = Dune::PDELab::ConvectionDiffusionFEM<Problem,FEM>;
  LOP lop(problem);

  auto bctype = Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<Problem>(gv,problem);

  auto g = Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<Problem>(gv,problem);

  // make constraints map and initialize it from a function
  typedef typename GFS::template ConstraintsContainer<R>::Type C;
  C cg;
  Dune::PDELab::constraints(bctype,gfs,cg);

  // make grid operator
  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,
                                     FMB,
                                     double,double,double,
                                     C,C> GridOperator;
  GridOperator gridoperator(gfs,cg,gfs,cg,lop,FMB(9));

  // make coefficent Vector and initialize it from a function
  // There is some weird shuffling around here - please leave it in,
  // it's there to test the copy constructor and assignment operator of the
  // matrix wrapper
  typedef typename GridOperator::Traits::Domain DV;
  DV x0(gfs,Dune::PDELab::Backend::unattached_container());
  {
    DV x1(gfs);
    DV x2(x1);
    x2 = 0.0;
    x0 = x1;
    x0 = x2;
  }
  x0 = 0.0;
  Dune::PDELab::interpolate(g,gfs,x0);
  Dune::PDELab::set_nonconstrained_dofs(cg,0.0,x0);
  native(x0).setChunkSize(chunk_size);


  // represent operator as a matrix
  // There is some weird shuffling around here - please leave it in,
  // it's there to test the copy constructor and assignment operator of the
  // matrix wrapper
  typedef typename GridOperator::Traits::Jacobian M;
  M m(gridoperator);
  m = 0.0;
  /*
    {
    M m1(gridoperator);
    M m2(m1);
    m2 = 0.0;
    m = m1;
    m = m2;
    }
  */
  gridoperator.jacobian(x0,m);
  //  Dune::printmatrix(std::cout,m.base(),"global stiffness matrix","row",9,1);

  // evaluate residual w.r.t initial guess
  typedef typename GridOperator::Traits::Range RV;
  RV r(gfs);
  native(r).setChunkSize(chunk_size);
  r = 0.0;
  gridoperator.residual(x0,r);
  native(x0).setChunkSize(chunk_size);

  std::cout << "norms start" << std::endl;
  std::cout << native(x0).two_norm2() << std::endl;
  std::cout << native(r).two_norm2() << std::endl;
  std::cout << "norms end" << std::endl;

  native(m).umv(native(r),native(x0));

#if 0
  // make ISTL solver
  Dune::MatrixAdapter<Native<M>, Native<DV>, Native<RV> > opa(native(m));
  //typedef Dune::PDELab::OnTheFlyOperator<typename DV::BaseT,typename RV::BaseT,GridOperator> ISTLOnTheFlyOperator;
  //ISTLOnTheFlyOperator opb(gridoperator);
  Dune::SeqSSOR<Native<M>, Native<DV>, Native<RV> > ssor(native(m),1,1.0);
  Dune::SeqILU0<Native<M>, Native<DV>, Native<RV> > ilu0(native(m),1.0);
  Dune::Richardson<Native<DV>, Native<RV> > richardson(1.0);

  Dune::CGSolver<Native<DV> > solvera(opa,ilu0,1E-10,5000,2);
  // FIXME: Use ISTLOnTheFlyOperator in the second solver again
  Dune::CGSolver<Native<DV> > solverb(opa,richardson,1E-10,5000,2);
  Dune::InverseOperatorResult stat;

  // solve the jacobian system
  r *= -1.0; // need -residual
  DV x(gfs,0.0);
  solvera.apply(native(x),native(r),stat);
  x += x0;

  // output grid function with VTKWriter
#endif
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
  Dune::PDELab::addSolutionToVTKWriter(vtkwriter,gfs,x0);
  vtkwriter.write(filename,Dune::VTK::ascii);
}

//===============================================================
// Main program with grid setup
//===============================================================

int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper::instance(argc, argv);

    int refinement = argc > 1 ? atoi(argv[1]) : 3;
    int chunk_size = argc > 2 ? atoi(argv[2]) : 16;

    // YaspGrid Q1 2D test
    {
      // make grid
      Dune::FieldVector<double,2> L(1.0);
      std::array<int,2> N = {{1,1}};
      Dune::YaspGrid<2> grid(L,N);
      grid.globalRefine(refinement);

      // get view
      typedef Dune::YaspGrid<2>::LeafGridView GV;
      auto gv=grid.leafGridView();

      // make finite element map
      typedef Dune::PDELab::QkLocalFiniteElementMap<GV,GV::ctype,double,1> FEM;
      FEM fem(gv);

      // solve problem
      poisson<GV,FEM,Dune::PDELab::ConformingDirichletConstraints,2>(gv,fem,"poisson_yasp_Q1_2d",chunk_size);
    }

    // YaspGrid Q2 2D test
    {
      // make grid
      Dune::FieldVector<double,2> L(1.0);
      std::array<int,2> N = {{1,1}};
      Dune::YaspGrid<2> grid(L,N);
      grid.globalRefine(refinement);

      // get view
      typedef Dune::YaspGrid<2>::LeafGridView GV;
      auto gv=grid.leafGridView();

      // make finite element map
      typedef Dune::PDELab::QkLocalFiniteElementMap<GV,GV::ctype,double,2> FEM;
      FEM fem(gv);

      // solve problem
      for (int i = 0; i < 10; ++i)
        poisson<GV,FEM,Dune::PDELab::ConformingDirichletConstraints,2>(gv,fem,"poisson_yasp_Q2_2d",chunk_size);
    }
#if 0
    // YaspGrid Q2 3D test
    {
      // make grid
      Dune::FieldVector<double,3> L(1.0);
      Dune::FieldVector<int,3> N(1);
      Dune::FieldVector<bool,3> B(false);
      Dune::YaspGrid<3> grid(L,N,B,0);
      grid.globalRefine(3);

      // get view
      typedef Dune::YaspGrid<3>::LeafGridView GV;
      const GV& gv=grid.leafView();

      // make finite element map
      typedef GV::Grid::ctype DF;
      typedef Dune::PDELab::Q1LocalFiniteElementMap<DF,double,3> FEM;
      FEM fem;

      // solve problem
      poisson<GV,FEM,Dune::PDELab::ConformingDirichletConstraints,2>(gv,fem,"poisson_yasp_Q1_3d");
    }

    // UG Pk 2D test
#if HAVE_UG
    {
      // make grid
      std::shared_ptr<Dune::UGGrid<2> > grid(TriangulatedUnitSquareMaker<Dune::UGGrid<2> >::create());
      grid->globalRefine(4);

      // get view
      typedef Dune::UGGrid<2>::LeafGridView GV;
      const GV& gv=grid->leafView();

      // make finite element map
      typedef GV::Grid::ctype DF;
      typedef double R;
      const int k=3;
      const int q=2*k;
      typedef Dune::PDELab::Pk2DLocalFiniteElementMap<GV,DF,double,k> FEM;
      FEM fem(gv);

      // solve problem
      poisson<GV,FEM,Dune::PDELab::ConformingDirichletConstraints,q>(gv,fem,"poisson_UG_Pk_2d");
    }
#endif

#if HAVE_ALBERTA
    {
      // make grid
      AlbertaUnitSquare grid;
      grid.globalRefine(8);

      // get view
      typedef AlbertaUnitSquare::LeafGridView GV;
      const GV& gv=grid.leafView();

      // make finite element map
      typedef GV::Grid::ctype DF;
      typedef double R;
      const int k=3;
      const int q=2*k;
      typedef Dune::PDELab::Pk2DLocalFiniteElementMap<GV,DF,double,k> FEM;
      FEM fem(gv);

      // solve problem
      poisson<GV,FEM,Dune::PDELab::ConformingDirichletConstraints,q>(gv,fem,"poisson_Alberta_Pk_2d");
    }
#endif

#if HAVE_ALUGRID
    {
      // make grid
      ALUUnitSquare grid;
      grid.globalRefine(4);

      // get view
      typedef ALUUnitSquare::LeafGridView GV;
      const GV& gv=grid.leafView();

      // make finite element map
      typedef GV::Grid::ctype DF;
      typedef double R;
      const int k=3;
      const int q=2*k;
      typedef Dune::PDELab::Pk2DLocalFiniteElementMap<GV,DF,double,k> FEM;
      FEM fem(gv);

      // solve problem
      poisson<GV,FEM,Dune::PDELab::ConformingDirichletConstraints,q>(gv,fem,"poisson_ALU_Pk_2d");
    }
#endif
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
