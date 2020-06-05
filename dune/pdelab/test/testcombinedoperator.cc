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

#include <dune/common/exceptions.hh> // We use exceptions
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/pdelab.hh>
#include <dune/pdelab/localoperator/sum.hh>
#include <dune/pdelab/localoperator/weightedsum.hh>

class Scalar :
    public Dune::PDELab::FullVolumePattern,
    public Dune::PDELab::LocalOperatorDefaultFlags,
    public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
public:
    // Pattern assembly flags
    enum { doPatternVolume = true };

    // Residual assembly flags
    enum { doAlphaVolume = true };

    Scalar (double value) : _v(value) {}

    template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
    void alpha_volume(const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
    {
        using size_type = typename LFSU::Traits::SizeType;
        for (size_type i=0; i<lfsu.size(); i++)
            r.accumulate(lfsv,i, _v);
    }

private:
    double _v;
};

// template <typename GFS>
// double test(const GO& go, const GFS& gfs, bool test_jacobian=true)
// {
//     // Initialize vectors and matrices for gridoperator calls
//     typename GO::Traits::Domain u(gfs,0.0);
//     typename GO::Traits::Range r(gfs);
//     typename GO::Traits::Jacobian jac(go);

//     // Call gridoperator methods
//     go.residual(u,r);
//     if (test_jacobian)
//       go.jacobian(u,jac);
//     go.jacobian_apply(u,r);

//     // For linear problems these methods should throw errors
//     bool jacobian_apply_error = false;
//     bool nonlinear_jacobian_apply_error = false;
//     try{ go.jacobian_apply(u,u,r); } catch (...) { jacobian_apply_error = true; }
//     try{ go.nonlinear_jacobian_apply(u,u,r); } catch (...) { nonlinear_jacobian_apply_error = true; }

//     // If both are true we hit exceptions where expected -> true is a sucess
//     return (jacobian_apply_error and nonlinear_jacobian_apply_error);
// }

int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper::instance(argc, argv);

    // Define parameters
    using Real = double;
    const unsigned int dim = 2;
    const int cells = 2;

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
    using FEM = Dune::PDELab::P0LocalFiniteElementMap<double,double,2>;
    FEM fem(Dune::GeometryTypes::quadrilateral);
    using CON = Dune::PDELab::NoConstraints;
    using VBE = Dune::PDELab::ISTL::VectorBackend<>;
    using GFS = Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE>;
    GFS gfs(gv,fem);

    // Matrix Backend
    using MBE = Dune::PDELab::ISTL::BCRSMatrixBackend<>;
    MBE mbe(9); // number of nonzeroes per row can be cross-checked by patternStatistics().

    // test results
    bool success = true;

    // Make sum local operator
    {
        using LOP = Dune::PDELab::InstationarySumLocalOperator<Scalar,Scalar>;
        LOP lop(Scalar(2.),Scalar(3.));
        using GO = Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,double,double,double>;

        GO go(gfs,gfs,lop,mbe);
        typename GO::Traits::Domain x(gfs,0.0);
        go.residual(x,x);
        std::cout << Dune::PDELab::Backend::native(x)[0] << std::endl;
        success &= (Dune::PDELab::Backend::native(x)[0] == 2. + 3.);
    }

    // This block tests a deprecated feature - we suppress the deprecation warning
    DUNE_NO_DEPRECATED_BEGIN
    // Make sum local operator from tuple
    {
        using OPS = std::tuple<Scalar,Scalar>;
        using OPSRefs = std::tuple<Scalar&,Scalar&>;
        using LOP = Dune::PDELab::InstationarySumLocalOperator<OPS>;
        Scalar s1(2.), s2(3.);
        OPSRefs ops(s1,s2);
        LOP lop(ops);

        using GO = Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,double,double,double>;
        GO go(gfs,gfs,lop,mbe);
        typename GO::Traits::Domain x(gfs,0.0);
        go.residual(x,x);
        std::cout << Dune::PDELab::Backend::native(x)[0] << std::endl;
        success &= (Dune::PDELab::Backend::native(x)[0] == 2. + 3.);
    }
    DUNE_NO_DEPRECATED_END

    // Make weighted sum local operator
    {
        using LOP = Dune::PDELab::WeightedSumLocalOperator<double,Scalar,Scalar>;
        LOP lop(Scalar(2.),Scalar(3.),{2.,-1.});

        using GO = Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,double,double,double>;
        GO go(gfs,gfs,lop,mbe);
        typename GO::Traits::Domain x(gfs,0.0);
        go.residual(x,x);
        std::cout << Dune::PDELab::Backend::native(x)[0] << std::endl;
        success &= (Dune::PDELab::Backend::native(x)[0] == 2*2. - 3.);
    }

    // This block tests a deprecated feature - we suppress the deprecation warning
    DUNE_NO_DEPRECATED_BEGIN
    // Make weighted sum local operator from tuple
    {
        using OPS = std::tuple<Scalar,Scalar>;
        using OPSRefs = std::tuple<Scalar&,Scalar&>;
        using Weights = Dune::FieldVector<double, 2>;
        using LOP = Dune::PDELab::WeightedSumLocalOperator<double,OPS>;
        Scalar s1(2.), s2(3.);
        OPSRefs ops(s1,s2);
        Weights w({2.,-1.});
        LOP lop(ops,w);

        using GO = Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,double,double,double>;
        GO go(gfs,gfs,lop,mbe);
        typename GO::Traits::Domain x(gfs,0.0);
        go.residual(x,x);
        std::cout << Dune::PDELab::Backend::native(x)[0] << std::endl;
        success &= (Dune::PDELab::Backend::native(x)[0] == 2*2. - 3.);
    }
    DUNE_NO_DEPRECATED_END

    return (! success);
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
