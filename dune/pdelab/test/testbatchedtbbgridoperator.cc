// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <algorithm>
#include <iostream>
#include <memory>
#include <ostream>
#include <string>
#include <thread>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/utility/partitioning/ranged.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/pdelab/backend/istl/bcrsmatrixbackend.hh>
#include <dune/pdelab/backend/istlsolverbackend.hh>
#include <dune/pdelab/backend/istlvectorbackend.hh>
#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/common/vtkexport.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/finiteelementmap/p0fem.hh>
#include <dune/pdelab/finiteelementmap/pkfem.hh>
#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/gridfunctionspace/vtk.hh>
#include <dune/pdelab/gridoperator/tbb.hh>
#include <dune/pdelab/localoperator/laplacedirichletp12d.hh>
#include <dune/pdelab/localoperator/poisson.hh>

#include"gridexamples.hh"

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
    if (x[0]>0.25 && x[0]<0.375 && x[1]>0.25 && x[1]<0.375)
      y = 50.0;
    else
      y = 0.0;
  }
};


// boundary grid function selecting boundary conditions
class ConstraintsParameters
  : public Dune::PDELab::DirichletConstraintsParameters
{

public:

  template<typename I>
  bool isDirichlet(const I & ig, const Dune::FieldVector<typename I::ctype, I::dimension-1> & x) const
  {
    Dune::FieldVector<typename I::ctype,I::dimension>
      xg = ig.geometry().global(x);

    if (xg[1]<1E-6 || xg[1]>1.0-1E-6)
      {
        return false;
      }
    if (xg[0]>1.0-1E-6 && xg[1]>0.5+1E-6)
      {
        return false;
      }
    return true;
  }

  template<typename I>
  bool isNeumann(const I & ig, const Dune::FieldVector<typename I::ctype, I::dimension-1> & x) const
  {
    return !isDirichlet(ig,x);
  }

};


// function for Dirichlet boundary conditions and initialization
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

// function for defining the flux boundary condition
template<typename GV, typename RF>
class J
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
                                                  J<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,J<GV,RF> > BaseT;

  J (const GV& gv) : BaseT(gv) {}
  inline void evaluateGlobal (const typename Traits::DomainType& x,
                              typename Traits::RangeType& y) const
  {
    if (x[1]<1E-6 || x[1]>1.0-1E-6)
      {
        y = 0;
        return;
      }
    if (x[0]>1.0-1E-6 && x[1]>0.5+1E-6)
      {
        y = -5.0;
        return;
      }
  }
};

//===============================================================
// Problem setup and solution
//===============================================================

// generate a P1 function and output it
template<typename GV, typename FEM, typename CON>
bool poisson (const GV& gv, const FEM& fem, std::string filename, int q)
{
  // constants and types
  typedef typename FEM::Traits::FiniteElementType::Traits::
    LocalBasisType::Traits::RangeFieldType R;

  // make function space
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,
    Dune::PDELab::ISTLVectorBackend<> > GFS;
  GFS gfs(gv,fem);
  gfs.name("solution");

  // make constraints map and initialize it from a function
  typedef typename GFS::template ConstraintsContainer<R>::Type C;
  C cg;
  cg.clear();
  ConstraintsParameters constraintsparameters;
  Dune::PDELab::constraints(constraintsparameters,gfs,cg);

  // make local operator
  typedef G<GV,R> GType;
  GType g(gv);
  typedef F<GV,R> FType;
  FType f(gv);
  typedef J<GV,R> JType;
  JType j(gv);
  typedef Dune::PDELab::Poisson<FType,ConstraintsParameters,JType> LOP;
  LOP lop(f,constraintsparameters,j, q);

#ifdef OLD_BACKEND
  typedef Dune::PDELab::ISTLMatrixBackend MBE;
  MBE mbe;
#else
  typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
  MBE mbe(27); // 27 is too large / correct for all test cases, so should work fine
#endif

  // make grid operator
  auto concurrency = 2;
  // auto concurrency = std::max(std::thread::hardware_concurrency(), 1u);
  typedef Dune::RangedPartitioning<GV, 0> Partitioning;
  auto partitions = concurrency;
  std::cout << "Using " << partitions << " partitions." << std::endl;
  typedef Dune::PDELab::BatchedTBBGridOperator<Partitioning,GFS,GFS,LOP,
                                               MBE,
                                               double,double,double,
                                               C,C> GridOperator;
  GridOperator gridoperator(gfs,cg,gfs,cg,lop,mbe);
  gridoperator.assembler().setPartitioning
    (std::make_shared<Partitioning>
     (gv, partitions));
  gridoperator.assembler().setVerbosity(1);

  // make coefficent Vector and initialize it from a function
  // There is some weird shuffling around here - please leave it in,
  // it's there to test the copy constructor and assignment operator of the
  // matrix wrapper
  typedef typename GridOperator::Traits::Domain DV;
  DV x0(gfs,Dune::PDELab::tags::unattached_container());
  {
    DV x1(gfs);
    DV x2(x1);
    x2 = 0.0;
    x0 = x1;
    x0 = x2;
  }

  Dune::PDELab::interpolate(g,gfs,x0);
  Dune::PDELab::set_nonconstrained_dofs(cg,0.0,x0);


  // represent operator as a matrix
  // There is some weird shuffling around here - please leave it in,
  // it's there to test the copy constructor and assignment operator of the
  // matrix wrapper
  typedef typename GridOperator::Traits::Jacobian M;
  M m;
  {
    Dune::Timer patternTimer;
    M m1(gridoperator);
    std::cout << "pattern creation:" << patternTimer.elapsed() << std::endl;
#ifndef OLD_BACKEND
    std::cout << m1.patternStatistics() << std::endl;
#endif
    M m2(m1);
    m2 = 0.0;
    m = m1;
    m = m2;
  }
  gridoperator.jacobian(x0,m);
  //  Dune::printmatrix(std::cout,m.base(),"global stiffness matrix","row",9,1);

  // evaluate residual w.r.t initial guess
  typedef typename GridOperator::Traits::Range RV;
  RV r(gfs);
  r = 0.0;
  gridoperator.residual(x0,r);

  using Dune::PDELab::Backend::Native;

  // make ISTL solver
  Dune::MatrixAdapter<Native<M>,Native<DV>,Native<RV> > opa(native(m));
  typedef Dune::PDELab::OnTheFlyOperator<DV,RV,
                                         GridOperator> ISTLOnTheFlyOperator;
  ISTLOnTheFlyOperator opb(gridoperator);
  Dune::SeqSSOR<Native<M>,Native<DV>,Native<RV> > ssor(native(m),1,1.0);
  Dune::SeqILU0<Native<M>,Native<DV>,Native<RV > ilu0(native(m),1.0);
  Dune::Richardson<Native<DV>,Native<RV> > richardson(1.0);
  Dune::Richardson<DV,RV> richardsonb(1.0);

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

  Dune::CGSolver<Native<DV> > solvera(opa,ilu0,1E-10,5000,2);
  Dune::CGSolver<DV> solverb(opb,richardsonb,1E-10,5000,2);
  Dune::InverseOperatorResult stat;

  std::cout << "Solve with assembled matrix" << std::endl;
  DV xa(gfs,0.0);
  {
    RV ra = r;
    ra *= -1.0; // need -residual
    solvera.apply(native(xa),native(ra),stat);
  }
  if(!stat.converged)
  {
    std::cerr << "Error: matrix solver did not converge" << std::endl;
    return false;
  }
  xa += x0;

  std::cout << "Solve with on-the-fly assembly" << std::endl;
  DV xb(gfs,0.0);
  {
    RV rb = r;
    rb *= -1.0; // need -residual
    solverb.apply(xb,rb,stat);
  }
  if(!stat.converged)
  {
    std::cerr << "Error: on-the-fly solver did not converge" << std::endl;
    return false;
  }
  xb += x0;

  // output grid function with VTKWriter
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
  Dune::PDELab::addSolutionToVTKWriter(vtkwriter,gfs,xa,
        Dune::PDELab::vtk::DefaultFunctionNameGenerator("matrix"));
  Dune::PDELab::addSolutionToVTKWriter(vtkwriter,gfs,xb,
        Dune::PDELab::vtk::DefaultFunctionNameGenerator("on_the_fly"));
  vtkwriter.write(filename,Dune::VTK::ascii);

  return true;
}

//===============================================================
// Main program with grid setup
//===============================================================

void accumulate_result(int &acc, bool good) {
  if(good && acc == 77)
    acc = 0;
  if(!good)
    acc = 1;
}

int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper::instance(argc, argv);

    int result = 77;

    {
      std::cout << "#### Testing YaspGrid Q1 2D" << std::endl;

      // make grid
      Dune::FieldVector<double,2> L(1.0);
      Dune::array<int,2> N(Dune::fill_array<int,2>(1));
      Dune::YaspGrid<2> grid(L,N);
      grid.globalRefine(3);

      // get view
      typedef Dune::YaspGrid<2>::LeafGridView GV;
      const GV& gv=grid.leafGridView();

      // make finite element map
      typedef GV::Grid::ctype DF;
      typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,double,1> FEM;
      FEM fem(gv);

      // solve problem
      bool good = poisson<GV,FEM,Dune::PDELab::ConformingDirichletConstraints>(gv,fem,"batchedtbbgridoperator_yasp_Q1_2d",2);
      accumulate_result(result, good);
    }

    {
      std::cout << "#### Testing YaspGrid Q2 2D" << std::endl;

      // make grid
      Dune::FieldVector<double,2> L(1.0);
      Dune::array<int,2> N(Dune::fill_array<int,2>(1));
      Dune::YaspGrid<2> grid(L,N);
      grid.globalRefine(3);

      // get view
      typedef Dune::YaspGrid<2>::LeafGridView GV;
      const GV& gv=grid.leafGridView();

      // make finite element map
      typedef GV::Grid::ctype DF;
      typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,double,2> FEM;
      FEM fem(gv);

      // solve problem
      bool good = poisson<GV,FEM,Dune::PDELab::ConformingDirichletConstraints>(gv,fem,"batchedtbbgridoperator_yasp_Q2_2d",2);
      accumulate_result(result, good);
    }

    {
      std::cout << "#### Testing YaspGrid Q1 3D" << std::endl;

      // make grid
      Dune::FieldVector<double,3> L(1.0);
      Dune::array<int,3> N(Dune::fill_array<int,3>(1));
      Dune::YaspGrid<3> grid(L,N);
      grid.globalRefine(3);

      // get view
      typedef Dune::YaspGrid<3>::LeafGridView GV;
      const GV& gv=grid.leafGridView();

      // make finite element map
      typedef GV::Grid::ctype DF;
      typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,double,1> FEM;
      FEM fem(gv);

      // solve problem
      bool good = poisson<GV,FEM,Dune::PDELab::ConformingDirichletConstraints>(gv,fem,"batchedtbbgridoperator_yasp_Q1_3d",2);
      accumulate_result(result, good);
    }

    {
      std::cout << "#### Testing YaspGrid Q2 3D" << std::endl;

      // make grid
      Dune::FieldVector<double,3> L(1.0);
      Dune::array<int,3> N(Dune::fill_array<int,3>(1));
      Dune::YaspGrid<3> grid(L,N);
      grid.globalRefine(3);

      // get view
      typedef Dune::YaspGrid<3>::LeafGridView GV;
      const GV& gv=grid.leafGridView();

      // make finite element map
      typedef GV::Grid::ctype DF;
      typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,double,2> FEM;
      FEM fem(gv);

      // solve problem
      bool good = poisson<GV,FEM,Dune::PDELab::ConformingDirichletConstraints>(gv,fem,"batchedtbbgridoperator_yasp_Q2_3d",2);
      accumulate_result(result, good);
    }

    // test passed
    return result;
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