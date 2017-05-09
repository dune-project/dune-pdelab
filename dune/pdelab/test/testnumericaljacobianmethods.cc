// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

/* Test numerical methods of localoperators

   A unit or integration test to see if the methods get called.  This
   does not test correctness of their implementation.  We want to test

   - numericaljacobian
   - numericaljacobianapply
   - numericalnonlinearjacobianapply
   - numericalresidual

   Note: We don't need to inherit from
   NumericalNonlinearJacobianApply* since the NumericalJacobianApply*
   method already does.
 */

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>

#include <dune/pdelab/localoperator/defaultimp.hh>
#include <dune/pdelab/localoperator/flags.hh>

#include <dune/grid/yaspgrid.hh>
#include<dune/pdelab/backend/istl.hh>

#include<dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/finiteelementmap/qkdg.hh>

class TestNumericalJacobiansLOP
  : public Dune::PDELab::NumericalJacobianVolume<TestNumericalJacobiansLOP>,
    public Dune::PDELab::NumericalJacobianVolumePostSkeleton<TestNumericalJacobiansLOP>,
    public Dune::PDELab::NumericalJacobianSkeleton<TestNumericalJacobiansLOP>,
    public Dune::PDELab::NumericalJacobianBoundary<TestNumericalJacobiansLOP>,
    public Dune::PDELab::NumericalJacobianApplyVolume<TestNumericalJacobiansLOP>,
    public Dune::PDELab::NumericalJacobianApplyVolumePostSkeleton<TestNumericalJacobiansLOP>,
    public Dune::PDELab::NumericalJacobianApplySkeleton<TestNumericalJacobiansLOP>,
    public Dune::PDELab::NumericalJacobianApplyBoundary<TestNumericalJacobiansLOP>,
    public Dune::PDELab::LocalOperatorDefaultFlags
{
public:
  static const bool doAlphaVolume = true;
  static const bool doAlphaVolumePostSkeleton = true;
  static const bool doAlphaSkeleton = true;
  static const bool doAlphaBoundary = true;

  template<typename EG, typename LFSU, typename X, typename LFSV,
           typename R>
  void alpha_volume
  ( const EG& eg,
    const LFSU& lfsu, const X& x, const LFSV& lfsv,
    R& r) const
  {}

  template<typename EG, typename LFSU, typename X, typename LFSV,
           typename R>
  void alpha_volume_post_skeleton
  ( const EG& eg,
    const LFSU& lfsu, const X& x, const LFSV& lfsv,
    R& r) const
  {}

  template<typename IG, typename LFSU, typename X, typename LFSV,
           typename R>
  void alpha_skeleton
  ( const IG& ig,
    const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
    const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
    R& r_s, R& r_n) const
  {}

  template<typename IG, typename LFSU, typename X, typename LFSV,
           typename R>
  void alpha_boundary
  ( const IG& ig,
    const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
    R& r_s) const
  {}


};


int main(int argc, char** argv)
{
  try{
    // Maybe initialize MPI
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
    using VBE = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed,blocksize>;
    using GFS = Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE>;
    GFS gfs(gv,fem);
    gfs.name("x_h");

    // Make local operator
    using LOP = TestNumericalJacobiansLOP;
    LOP lop;

    // Make grid operator
    using MBE = Dune::PDELab::ISTL::BCRSMatrixBackend<>;
    MBE mbe(9); // number of nonzeroes per row can be cross-checked by patternStatistics().
    using CC = typename GFS::template ConstraintsContainer<Real>::Type;
    CC cc;
    using GO = Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,Real,Real,Real,CC,CC>;
    GO go(gfs,cc,gfs,cc,lop,mbe);

    // Initialize vectors and matrices for gridoperator calls
    typedef typename GO::Traits::Domain U;
    U u(gfs,0.0);
    using R = typename GO::Traits::Range;
    R r(gfs);
    using J = typename GO::Traits::Jacobian;
    J jac(go);

    // Call gridoperator methods
    go.residual(u,r);
    go.jacobian(u,jac);
    go.jacobian_apply(u,r);
    go.nonlinear_jacobian_apply(u,u,r);

    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
