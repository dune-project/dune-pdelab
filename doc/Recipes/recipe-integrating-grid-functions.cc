#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>

#include <dune/common/filledarray.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/pdelab.hh>

/**
 * \page recipe-integrating-grid-functions Integrating grid functions
 * First, we define some example grid function to integrate.
 * \snippet recipe-integrating-grid-functions.cc Defining an analytic grid function
 *
 * Computing the integral of a grid function is just a single call:
 * \snippet recipe-integrating-grid-functions.cc Compute integral
 *
 * When run in parallel, this integrates over the local interior domain.
 * So, in that case, we still need to sum over all processes for the full integral.
 * \snippet recipe-integrating-grid-functions.cc Sum for parallel case
 *
 * Full example code: @ref recipe-integrating-grid-functions.cc
 * \example recipe-integrating-grid-functions.cc
 * See explanation at @ref recipe-integrating-grid-functions
 */

int main(int argc, char** argv)
{
  try{
    // Initialize Mpi
    Dune::MPIHelper::instance(argc, argv);

    // need a grid in order to test grid functions
    constexpr unsigned int dim = 2;
    Dune::FieldVector<double,dim> L(1.0);
    std::array<int,dim> N(Dune::filledArray<dim,int>(64));

    typedef Dune::YaspGrid<dim> Grid;
    Grid grid(L,N);

    // [Defining an analytic grid function]
    auto analyticFunction = Dune::PDELab::makeGridFunctionFromCallable (grid.leafGridView(), [&](const auto& i, const auto& x){
      return exp(-(x*x));
    });
    //! [Defining an analytic grid function]

    // [Compute integral]
    auto integral = Dune::PDELab::integrateGridFunction(analyticFunction,10);
    //! [Compute integral]

    // [Sum for parallel case]
    integral = grid.leafGridView().comm().sum(integral);
    //! [Sum for parallel case]

    std::cout << "Integral: " << integral << std::endl;

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
