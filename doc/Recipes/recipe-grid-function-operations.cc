#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>

#include <dune/common/filledarray.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/pdelab.hh>

/**
 * \page recipe-grid-function-operations Operations on grid functions
 *
 * A frequent example of pointwise operations on GridFunctions is calculating
 * error norms (e.g. analytical vs numerical solution). Here we take
 * the squared difference of two functions. However, the same approach applies to any
 * operation.
 *
 * First, let's define some grid functions for testing.
 * \snippet recipe-grid-function-operations.cc Defining functions
 *
 * Now, we can either plug them into an existing adapter (itself a GridFunction) computing the difference squared
 * \snippet recipe-grid-function-operations.cc Instantiating DifferenceSquaredAdapter
 *
 * or define it ourselves. Here we use the per-element evaluation interface receiving an element and local
 * coordinates on the reference element.
 * \snippet recipe-grid-function-operations.cc Defining operations
 *
 * We can use the resulting GridFunctions as always, e.g. integrate them:
 * \snippet recipe-grid-function-operations.cc Compute integral
 *
 * Full example code: @ref recipe-grid-function-operations.cc
 * \example recipe-grid-function-operations.cc
 * See explanation at @ref recipe-grid-function-operations
 */


int main(int argc, char** argv)
{
  try{
    // Initialize Mpi
    Dune::MPIHelper::instance(argc, argv);

    // need a grid in order to test grid functions
    constexpr unsigned int dim = 2;
    Dune::FieldVector<double,dim> L(5.0);
    std::array<int,dim> N(Dune::filledArray<dim,int>(64));

    typedef Dune::YaspGrid<dim> Grid;
    Grid grid(L,N);

    // [Defining functions]
    auto analyticFunction1 = Dune::PDELab::makeGridFunctionFromCallable (grid.leafGridView(), [&](const auto& x){
      return exp(-(x*x));
    });
    auto analyticFunction2 = Dune::PDELab::makeGridFunctionFromCallable (grid.leafGridView(), [&](const auto& x){
      return x[0]*3.0;
    });
    //! [Defining functions]


    {
      // [Instantiating DifferenceSquaredAdapter]
      Dune::PDELab::DifferenceSquaredAdapter<decltype(analyticFunction1),decltype(analyticFunction2)> difference(analyticFunction1, analyticFunction2);
      // Note: Newer compiler versions do not require template arguments here!
      //! [Instantiating DifferenceSquaredAdapter]
      auto integral = Dune::PDELab::integrateGridFunction(difference,10);
      std::cout << "Integral: " << integral << std::endl;
    }

    {
      // [Defining operations]
      auto difference = Dune::PDELab::makeGridFunctionFromCallable (grid.leafGridView(), [&](const auto& element, const auto& xlocal){
        Dune::FieldVector<double,1> y1;
        analyticFunction1.evaluate(element, xlocal, y1);
        Dune::FieldVector<double,1> y2;
        analyticFunction2.evaluate(element, xlocal, y2);
        y2 -= y1;
        return y2.two_norm2();
      });
      //! [Defining operations]
      // [Compute integral]
      auto integral = Dune::PDELab::integrateGridFunction(difference,10);
      //! [Compute integral]
      std::cout << "Integral: " << integral << std::endl;
    }


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
