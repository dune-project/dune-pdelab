// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>

#include <dune/common/exceptions.hh>

#include <dune/pdelab/gridoperator/common/localmatrix.hh>

void test()
{
  Dune::PDELab::LocalMatrix<double> m(2, 0);
  for (auto it = m.begin(); it != m.end(); ++it) {
    DUNE_THROW(Dune::Exception, "local matrix iterates over 0 columns");
  }
}

int main(int argc, char** argv)
{
  try {
    test();
    // test passed
    return 0;

  } catch (Dune::Exception& e) {
    std::cerr << "Dune reported error: " << e << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
    return 1;
  }
}
