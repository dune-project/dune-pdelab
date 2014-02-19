// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cmath>
#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/shared_ptr.hh>

#include <dune/geometry/generalvertexorder.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/utility/vertexorderfactory.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/pdelab/backend/istlvectorbackend.hh>
#include <dune/pdelab/backend/istl/bcrsmatrixbackend.hh>
#include <dune/pdelab/backend/istlsolverbackend.hh>
#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/constraints/common/constraintsparameters.hh>
#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/finiteelementmap/q1fem.hh>
#include <dune/pdelab/finiteelementmap/q22dfem.hh>
#include <dune/pdelab/finiteelementmap/pk2dfem.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/localoperator/poisson.hh>

#include <dune/pdelab/gridfunctionspace/vtk.hh>

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
class F :
  public Dune::PDELab::AnalyticGridFunctionBase<
    Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
    F<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,F<GV,RF> > BaseT;

  F (const GV& gv) : BaseT(gv) {}
  inline void evaluateGlobal(const typename Traits::DomainType& x,
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
class G :
  public Dune::PDELab::AnalyticGridFunctionBase<
    Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
    G<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,G<GV,RF> > BaseT;

  G (const GV& gv) : BaseT(gv) {}
  inline void evaluateGlobal(const typename Traits::DomainType& x,
                             typename Traits::RangeType& y) const
  {
    typename Traits::DomainType center;
    for (int i=0; i<GV::dimension; i++) center[i] = 0.5;
    center -= x;
    y = std::exp(-center.two_norm2());
  }
};

// function for defining the flux boundary condition
template<typename GV, typename RF>
class J :
  public Dune::PDELab::AnalyticGridFunctionBase<
    Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
    J<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,J<GV,RF> > BaseT;

  J (const GV& gv) : BaseT(gv) {}
  inline void evaluateGlobal(const typename Traits::DomainType& x,
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
void poisson (const GV& gv, const FEM& fem, std::string filename, int q)
{
  // constants and types
  typedef typename FEM::Traits::FiniteElementType::Traits::Basis::Traits::
    RangeField R;

  // make function space
  typedef Dune::PDELab::GridFunctionSpace<
    GV,FEM,CON,
    Dune::PDELab::ISTLVectorBackend<>
    > GFS;
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
  LOP lop(f,constraintsparameters,j,q);

  // make grid operator
  typedef Dune::PDELab::GridOperator<
    GFS,GFS,LOP,
    Dune::PDELab::ISTLMatrixBackend,
    R,R,R,C,C
    > GridOperator;
  GridOperator gridoperator(gfs,cg,gfs,cg,lop);

  // make coefficent Vector and initialize it from a function
  typedef typename GridOperator::Traits::Domain V;
  V x0(gfs,0.0);

  Dune::PDELab::interpolate(g,gfs,x0);
  Dune::PDELab::set_nonconstrained_dofs(cg,0.0,x0);


  // represent operator as a matrix
  typedef typename GridOperator::Traits::Jacobian M;
  M m(gridoperator);
  m = 0.0;
  gridoperator.jacobian(x0,m);
  // Dune::printmatrix(std::cout,m.base(),"global stiffness matrix","row",9,1);

  // evaluate residual w.r.t initial guess
  V r(gfs);
  r = 0.0;
  gridoperator.residual(x0,r);

  V x(gfs,0.0);

  {
    // make ISTL solver
    using Dune::PDELab::istl::raw;
    using Dune::PDELab::istl::raw_type;

    Dune::MatrixAdapter<
      typename raw_type<M>::type,
      typename raw_type<V>::type,
      typename raw_type<V>::type
      > opa(raw(m));

    Dune::SeqILU0<
      typename raw_type<M>::type,
      typename raw_type<V>::type,
      typename raw_type<V>::type
      > ilu0(raw(m),1.0);

    Dune::CGSolver<typename raw_type<V>::type> solvera(opa,ilu0,1E-10,5000,2);
    Dune::InverseOperatorResult stat;

    // solve the jacobian system
    r *= -1.0; // need -residual
    solvera.apply(raw(x),raw(r),stat);
    x += x0;
  }

  // output grid function with VTKWriter
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
  Dune::PDELab::addSolutionToVTKWriter(vtkwriter,gfs,x);
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

    // YaspGrid Q1 2D test
    {
      // make grid
      Dune::FieldVector<double,2> L(1.0);
      Dune::array<int,2> N(Dune::fill_array<int,2>(1));
      Dune::YaspGrid<2> grid(L,N);
      grid.globalRefine(3);

      // get view
      typedef Dune::YaspGrid<2>::LeafGridView GV;
      const GV& gv=grid.leafGridView();

      // make finite element map
      typedef Dune::PDELab::Q1FiniteElementMap<GV::Codim<0>::Geometry, double>
        FEM;
      FEM fem;

      // solve problem
      poisson<GV,FEM,Dune::PDELab::ConformingDirichletConstraints>
        (gv,fem,"poisson_globalfe_yasp_Q1_2d",2);
    }

    // YaspGrid Q2 2D test
    {
      // make grid
      Dune::FieldVector<double,2> L(1.0);
      Dune::array<int,2> N(Dune::fill_array<int,2>(1));
      Dune::YaspGrid<2> grid(L,N);
      grid.globalRefine(3);

      // get view
      typedef Dune::YaspGrid<2>::LeafGridView GV;
      const GV& gv=grid.leafGridView();

      // make finite element map
      typedef Dune::PDELab::Q22DFiniteElementMap<
        GV::Codim<0>::Geometry, double
        > FEM;
      FEM fem;

      // solve problem
      poisson<GV,FEM,Dune::PDELab::ConformingDirichletConstraints>
        (gv,fem,"poisson_globalfe_yasp_Q2_2d",2);
    }

    // YaspGrid Q1 3D test
    {
      // make grid
      Dune::FieldVector<double,3> L(1.0);
      Dune::array<int,3> N(Dune::fill_array<int,3>(1));
      Dune::YaspGrid<3> grid(L,N);
      grid.globalRefine(3);

      // get view
      typedef Dune::YaspGrid<3>::LeafGridView GV;
      const GV& gv=grid.leafGridView();

      // make finite element map
      typedef Dune::PDELab::Q1FiniteElementMap<GV::Codim<0>::Geometry, double>
        FEM;
      FEM fem;

      // solve problem
      poisson<GV,FEM,Dune::PDELab::ConformingDirichletConstraints>
        (gv,fem,"poisson_globalfe_yasp_Q1_3d",2);
    }

    // UG Pk 2D test
#if HAVE_UG
    {
      // make grid
      Dune::shared_ptr<Dune::UGGrid<2> >
        grid(TriangulatedUnitSquareMaker<Dune::UGGrid<2> >::create());
      grid->globalRefine(4);

      // get view
      typedef Dune::UGGrid<2>::LeafGridView GV;
      const GV& gv=grid->leafGridView();

      // make finite element map
      const int k=3;
      const int q=2*k;
      typedef Dune::VertexOrderByIdFactory<GV::Grid::GlobalIdSet>
        VOFactory;
      VOFactory voFactory(grid->globalIdSet());
      typedef Dune::PDELab::Pk2DFiniteElementMap<
        GV::Codim<0>::Geometry, VOFactory, double, k
        > FEM;
      FEM fem(voFactory);

      // solve problem
      poisson<GV,FEM,Dune::PDELab::ConformingDirichletConstraints>
        (gv,fem,"poisson_globalfe_UG_Pk_2d",q);
    }
#endif

#if HAVE_ALBERTA
    {
      // make grid
      AlbertaUnitSquare grid;
      grid.globalRefine(8);

      // get view
      typedef AlbertaUnitSquare::LeafGridView GV;
      const GV& gv=grid.leafGridView();

      // make finite element map
      const int k=3;
      const int q=2*k;
      typedef Dune::VertexOrderByIdFactory<GV::Grid::GlobalIdSet>
        VOFactory;
      VOFactory voFactory(grid.globalIdSet());
      typedef Dune::PDELab::Pk2DFiniteElementMap<
        GV::Codim<0>::Geometry, VOFactory, double, k
        > FEM;
      FEM fem(voFactory);

      // solve problem
      poisson<GV,FEM,Dune::PDELab::ConformingDirichletConstraints>
        (gv,fem,"poisson_globalfe_Alberta_Pk_2d",q);
    }
#endif

#if HAVE_ALUGRID
    {
      // make grid
      ALUUnitSquare grid;
      grid.globalRefine(4);

      // get view
      typedef ALUUnitSquare::LeafGridView GV;
      const GV& gv=grid.leafGridView();

      // make finite element map
      const int k=3;
      const int q=2*k;
      typedef Dune::VertexOrderByIdFactory<GV::Grid::GlobalIdSet>
        VOFactory;
      VOFactory voFactory(grid.globalIdSet());
      typedef Dune::PDELab::Pk2DFiniteElementMap<
        GV::Codim<0>::Geometry, VOFactory, double, k
        > FEM;
      FEM fem(voFactory);

      // solve problem
      poisson<GV,FEM,Dune::PDELab::ConformingDirichletConstraints>
        (gv,fem,"poisson_globalfe_ALU_Pk_2d",q);
    }
#endif

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
