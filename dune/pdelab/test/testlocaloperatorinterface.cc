// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

/* Test Localoperator Interface

   This test does the following:

   - Create linear and nonlinear LOP from localoperator/interface.hh
     that does nothing but calls implements all localoperator methods.
   - Create jacobian matrix -> tests pattern calls
   - Create several grid operator variants and test the calls to the
     residual(), jacobian(), jacobian_apply() methods.
     The latter is tested for both linear and nonlinear LOP.

   The test has two purposes:
   - See if everything compiles when we call all these methods
   - Use the interface localoperator in order to find bugs if
     interfaces change
 */

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/pdelab/gridoperator/fastdg.hh>
#include <dune/pdelab/gridoperator/onestep.hh>
#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/finiteelementmap/qkdg.hh>
#include <dune/pdelab/localoperator/interface.hh>

int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper::instance(argc, argv);

    // Define parameters
    using Real = double;
    const unsigned int dim = 2;
    const int cells = 8;
    const int degree = 2;

    // Create grid
    Dune::FieldVector<Real,dim> l(1.0);
    Dune::array<int,dim> s;
    std::fill(s.begin(), s.end(), cells);
    std::bitset<dim> p(0);
    int overlap = 0;
    using Grid = Dune::YaspGrid<dim>;
    Grid grid(l,s,p,overlap);

    // Get grid view
    using GV = Grid::LeafGridView;
    GV gv = grid.leafGridView();

    // Make grid function space
    using FEM = Dune::PDELab::QkDGLocalFiniteElementMap<Grid::ctype,Real,degree,dim>;
    FEM fem;
    using CON = Dune::PDELab::NoConstraints;
    const int blocksize = Dune::QkStuff::QkSize<degree,dim>::value;
    using VBE = Dune::PDELab::istl::VectorBackend<Dune::PDELab::istl::Blocking::fixed,blocksize>;
    using GFS = Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE>;
    GFS gfs(gv,fem);
    gfs.name("x_h");

    // Make nonlinear local operator
    using LOP0 = Dune::PDELab::LocalOperatorInterface<false>;
    LOP0 lop0;
    // Make linear local operator
    using LOP1 = Dune::PDELab::LocalOperatorInterface<>;
    LOP1 lop1;

    // Make grid operators
    using MBE = Dune::PDELab::istl::BCRSMatrixBackend<>;
    MBE mbe(9); // number of nonzeroes per row can be cross-checked by patternStatistics().
    using CC = typename GFS::template ConstraintsContainer<Real>::Type;
    CC cc;
    // nonlinear grid operator
    using GO0 = Dune::PDELab::GridOperator<GFS,GFS,LOP0,MBE,Real,Real,Real,CC,CC>;
    GO0 go0(gfs,cc,gfs,cc,lop0,mbe);
    // nonlinear fast DG grid operator
    using FastDGGO = Dune::PDELab::FastDGGridOperator<GFS,GFS,LOP0,MBE,Real,Real,Real,CC,CC>;
    FastDGGO fastdggo(gfs,cc,gfs,cc,lop0,mbe);
    // linear grid operator
    using GO1 = Dune::PDELab::GridOperator<GFS,GFS,LOP1,MBE,Real,Real,Real,CC,CC>;
    GO1 go1(gfs,cc,gfs,cc,lop1,mbe);
    // one-step grid operator
    using OneStepGO = Dune::PDELab::OneStepGridOperator<GO0,GO1>;
    OneStepGO onestepgo(go0,go1);

    // Initialize vectors and matrices for gridoperator calls
    typedef typename GO0::Traits::Domain U;
    U u(gfs,0.0);
    using R = typename GO0::Traits::Range;
    R r(gfs);
    using J = typename GO0::Traits::Jacobian;
    J jac(go0);

    // Call gridoperator methods
    go0.residual(u,r);
    go0.jacobian(u,jac);
    go0.jacobian_apply(u,u,r,Dune::Direction::forward);
    fastdggo.residual(u,r);
    fastdggo.jacobian(u,jac);
    fastdggo.jacobian_apply(u,u,r,Dune::Direction::forward);
    go1.jacobian_apply(u,r); // test deprecated call for backwards compatibility
    go1.jacobian_apply(u,u,r,Dune::Direction::forward);
    onestepgo.jacobian_apply(u,u,r,Dune::Direction::forward);

    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
