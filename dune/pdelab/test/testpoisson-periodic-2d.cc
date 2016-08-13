// -*- tab-width: 4; indent-tabs-mode: nil -*-

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/pdelab/test/testpoisson-periodic.hh>

int main(int argc, char **argv)
{
  try {

    // initialize MPI, finalize is done automatically on exit
    Dune::MPIHelper::instance(argc,argv);

    // define parameters
    const unsigned int dim = 2;
    typedef double NumberType;

    // Define periodic direction
    std::bitset<dim> periodic (false);
    periodic[0] = true;

    // Set periodic overlap size
    int overlap = 1;

    // Building a non-symmetric grid to make sure coordinate transformations work correctly
    std::array<std::vector<NumberType>,dim> coords;
    for (int i=0; i<2; i++)
    {
      coords[i].resize(9);
      coords[i][0] = 0.0;
      coords[i][1] = 0.25;
      coords[i][2] = 0.375;
      coords[i][3] = 0.4375;
      coords[i][4] = 0.5;
      coords[i][5] = 0.65;
      coords[i][6] = 0.8;
      coords[i][7] = 0.9;
      coords[i][8] = 1.0;
    }

    // make grid
    typedef Dune::YaspGrid<dim, Dune::TensorProductCoordinates<NumberType,dim> > GM;
    GM grid (coords,periodic,overlap);

    // Solve problem
    poisson<GM, NumberType, dim> (grid);

  } catch (Dune::Exception &e) {
    std::cerr << e << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Generic exception!" << std::endl;
    return 2;
  }
  // done
  return 0;
}
