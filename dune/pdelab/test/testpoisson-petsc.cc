// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include<iostream>
#include<vector>
#include<map>
#include<string>
#include<dune/common/mpihelper.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include <dune/common/shared_ptr.hh>
#include<dune/common/static_assert.hh>
#include<dune/grid/yaspgrid.hh>
#include<dune/istl/bvector.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/preconditioners.hh>
#include<dune/istl/io.hh>
#include<dune/istl/paamg/amg.hh>

#include"../finiteelementmap/q12dfem.hh"
#include"../finiteelementmap/q22dfem.hh"
#include"../finiteelementmap/q1fem.hh"
#include"../finiteelementmap/conformingconstraints.hh"
#include"../gridfunctionspace/gridfunctionspace.hh"
#include"../gridfunctionspace/gridfunctionspaceutilities.hh"
#include"../gridfunctionspace/interpolate.hh"
#include"../constraints/constraints.hh"
#include"../common/function.hh"
#include"../common/vtkexport.hh"
#include"../gridoperator/gridoperator.hh"
#include"../backend/petscvectorbackend.hh"
#include"../backend/petscmatrixbackend.hh"
#include"../localoperator/laplacedirichletp12d.hh"
#include"../localoperator/poisson.hh"

// Petsc solvers
#include <petscksp.h>

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
template<typename GV, typename FEM, typename CON, int q>
void poisson (const GV& gv, const FEM& fem, std::string filename)
{
  // constants and types
  typedef typename GV::Grid::ctype DF;
  typedef typename FEM::Traits::FiniteElementType::Traits::
    LocalBasisType::Traits::RangeFieldType R;

  // make function space
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,
    Dune::PDELab::PetscVectorBackend > GFS;
  GFS gfs(gv,fem);

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
  typedef Dune::PDELab::Poisson<FType,ConstraintsParameters,JType,q> LOP;
  LOP lop(f,constraintsparameters,j);

  // make grid operator
  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,
                                     Dune::PDELab::PetscMatrixBackend,
                                     double,double,double,
                                     C,C> GridOperator;
  GridOperator gridoperator(gfs,cg,gfs,cg,lop);

  // make coefficent Vector and initialize it from a function
  typedef typename GridOperator::Traits::Domain DV;
  DV x(gfs,0.0);

  Dune::PDELab::interpolate(g,gfs,x);
  //Dune::PDELab::set_nonconstrained_dofs(cg,0.0,x0);


  // represent operator as a matrix
  typedef typename GridOperator::Traits::Jacobian M;
  M m(gridoperator);
  gridoperator.jacobian(x,m);

  // evaluate residual w.r.t initial guess
  typedef typename GridOperator::Traits::Range RV;
  RV r(gfs,0.0);
  gridoperator.residual(x,r);

  DV z(gfs);

  KSP ksp;
  PETSC_CALL(KSPCreate(PETSC_COMM_SELF,&ksp));
  PETSC_CALL(KSPSetOperators(ksp,m,m,DIFFERENT_NONZERO_PATTERN));
  PETSC_CALL(KSPSolve(ksp,r,z));
  KSPConvergedReason reason;
  PETSC_CALL(KSPGetConvergedReason(ksp,&reason));
  std::cout << "convergence: " << reason << std::endl;
  int its;
  PETSC_CALL(KSPGetIterationNumber(ksp,&its));
  std::cout << "iterations: " << its << std::endl;
  PETSC_CALL(KSPDestroy(&ksp));

  // apply correction
  x -= z;

  // make discrete function object
  typedef Dune::PDELab::DiscreteGridFunction<GFS,DV> DGF;
  DGF dgf(gfs,x);

  // output grid function with VTKWriter
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(dgf,"solution"));
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

    // call after MPIHelper
    PETSC_CALL(PetscInitialize(&argc,&argv,NULL,NULL));

    // YaspGrid Q1 2D test
    {
      // make grid
      Dune::FieldVector<double,2> L(1.0);
      Dune::FieldVector<int,2> N(1);
      Dune::FieldVector<bool,2> B(false);
      Dune::YaspGrid<2> grid(L,N,B,0);
      grid.globalRefine(3);

      // get view
      typedef Dune::YaspGrid<2>::LeafGridView GV;
      const GV& gv=grid.leafView();

      // make finite element map
      typedef GV::Grid::ctype DF;
      typedef Dune::PDELab::Q12DLocalFiniteElementMap<DF,double> FEM;
      FEM fem;

      // solve problem
      poisson<GV,FEM,Dune::PDELab::ConformingDirichletConstraints,2>(gv,fem,"poisson_petsc_yasp_Q1_2d");
    }

    // YaspGrid Q2 2D test
    {
      // make grid
      Dune::FieldVector<double,2> L(1.0);
      Dune::FieldVector<int,2> N(1);
      Dune::FieldVector<bool,2> B(false);
      Dune::YaspGrid<2> grid(L,N,B,0);
      grid.globalRefine(3);

      // get view
      typedef Dune::YaspGrid<2>::LeafGridView GV;
      const GV& gv=grid.leafView();

      // make finite element map
      typedef GV::Grid::ctype DF;
      typedef Dune::PDELab::Q22DLocalFiniteElementMap<DF,double> FEM;
      FEM fem;

      // solve problem
      poisson<GV,FEM,Dune::PDELab::ConformingDirichletConstraints,2>(gv,fem,"poisson_petsc_yasp_Q2_2d");
    }

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
      poisson<GV,FEM,Dune::PDELab::ConformingDirichletConstraints,2>(gv,fem,"poisson_petsc_yasp_Q1_3d");
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

  // call after leaving try-block to make sure all wrappers are destructed
  // before this call
  PETSC_CALL(PetscFinalize());

}
