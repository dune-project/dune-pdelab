// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

/* Test Localoperator Interface

   This test does the following:

   - Create LOP from localoperator/interface.hh that does nothing but
     calls implements all localoperator methods.
   - Create jacobian matrix -> tests pattern calls
   - Call residual, jacobian, jacobian_apply and
     nonlinear_jacobian_apply methods from the gridoperator

   The test has two purposes:
   - See if everything compiles when we call all these methods
   - Use the interface localoperator in order to find bugs if
     interfaces change
 */

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <iostream>

#include <dune/common/deprecated.hh>
#include <dune/common/exceptions.hh> // We use exceptions
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/pdelab.hh>

class NonLinearLocalOperatorInterface
  : public Dune::PDELab::LocalOperatorInterface
{
public:
  static constexpr bool isLinear = false;
};


template <typename GO, typename GFS>
bool test_linear(const GO& go, const GFS& gfs, bool test_jacobian=true)
{
    // Initialize vectors and matrices for gridoperator calls
    typename GO::Traits::Domain u(gfs,0.0);
    typename GO::Traits::Range r(gfs);
    typename GO::Traits::Jacobian jac(go);

    // Call gridoperator methods
    go.residual(u,r);
    if (test_jacobian)
      go.jacobian(u,jac);
    go.jacobian_apply(u,r);

    // For linear problems these methods should throw errors
    bool jacobian_apply_error = false;
    bool nonlinear_jacobian_apply_error = false;
    try{ go.jacobian_apply(u,u,r); } catch (...) { jacobian_apply_error = true; }
    DUNE_NO_DEPRECATED_BEGIN
    try{ go.nonlinear_jacobian_apply(u,u,r); } catch (...) { nonlinear_jacobian_apply_error = true; }
    DUNE_NO_DEPRECATED_END

    // If both are true we hit exceptions where expected -> true is a sucess
    return (jacobian_apply_error and nonlinear_jacobian_apply_error);
}

template <typename GO, typename GFS>
bool test_nonlinear(const GO& go, const GFS& gfs, bool test_jacobian=true)
{
    // Initialize vectors and matrices for gridoperator calls
    typename GO::Traits::Domain u(gfs,0.0);
    typename GO::Traits::Range r(gfs);
    typename GO::Traits::Jacobian jac(go);

    // Call gridoperator methods
    go.residual(u,r);
    if (test_jacobian)
      go.jacobian(u,jac);
    go.jacobian_apply(u,u,r);

    DUNE_NO_DEPRECATED_BEGIN
    go.nonlinear_jacobian_apply(u,u,r);
    DUNE_NO_DEPRECATED_END

    // For non linear problems this methods should throw errors
    bool jacobian_apply_error = false;
    try{ go.jacobian_apply(u,r); } catch (...) { jacobian_apply_error = true; }

    return jacobian_apply_error;
}


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
    std::array<int,dim> s;
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
    using VBE = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed, blocksize>;
    using GFS = Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE>;
    GFS gfs(gv,fem);
    gfs.name("x_h");

    // Constraints
    using MBE = Dune::PDELab::ISTL::BCRSMatrixBackend<>;
    MBE mbe(9); // number of nonzeroes per row can be cross-checked by patternStatistics().
    using CC = typename GFS::template ConstraintsContainer<Real>::Type;
    CC cc;

    // Make linear local operator
    using LOP = Dune::PDELab::LocalOperatorInterface;
    LOP lop;

    // Make non linear local operator
    using NLLOP = NonLinearLocalOperatorInterface;
    NLLOP nllop;

    //=========================================
    // Test interface for linear local operator
    //=========================================
    using GO = Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,Real,Real,Real,CC,CC>;
    GO go(gfs,cc,gfs,cc,lop,mbe);
    auto linear_test = test_linear(go, gfs);

    //==============================================
    // Test interface for non linear local operators
    //==============================================
    using NLGO = Dune::PDELab::GridOperator<GFS,GFS,NLLOP,MBE,Real,Real,Real,CC,CC>;
    NLGO nlgo(gfs,cc,gfs,cc,nllop,mbe);
    auto nonlinear_test = test_nonlinear(nlgo, gfs);

    //===============================================
    // Test FastDGOperator with linear local operator
    //===============================================
    using FastGO = Dune::PDELab::FastDGGridOperator<GFS,GFS,LOP,MBE,Real,Real,Real,CC,CC>;
    FastGO fast_go(gfs,cc,gfs,cc,lop,mbe);
    auto fast_linear_test = test_linear(fast_go, gfs, false);

    //===================================================
    // Test FastDGOperator with non linear local operator
    //===================================================
    using NLFastGO = Dune::PDELab::FastDGGridOperator<GFS,GFS,NLLOP,MBE,Real,Real,Real,CC,CC>;
    NLFastGO nl_fast_go(gfs,cc,gfs,cc,nllop,mbe);
    auto fast_nonlinear_test = test_nonlinear(nl_fast_go, gfs, false);

    //==================
    // Return testresult
    //==================
    bool testfail = true;
    if (linear_test && nonlinear_test && fast_linear_test && fast_nonlinear_test)
      testfail = false;
    return testfail;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
