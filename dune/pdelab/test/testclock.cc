// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include <ostream>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/pdelab/common/clock.hh>

int main(int argc, char** argv) {
  try{
    int result = 0;

    //Maybe initialize Mpi
    Dune::MPIHelper::instance(argc, argv);

    // Check Wall Time
    std::cout << "Testing wall time" << std::endl;
    try {
      std::cout << "  Wall time implemented by: "
                << Dune::PDELab::getWallTimeImp() << std::endl;
      std::cout << "  Wall time resolution: "
                << Dune::PDELab::getWallTimeResolution() << std::endl;
      std::cout << "  Current wall time: " << Dune::PDELab::getWallTime()
                << std::endl;
    }
    catch(Dune::PDELab::ClockError &e) {
      std::cerr << "  Cought ClockError: " << e << std::endl;
      result = 1;
    }

    std::cout << std::endl;

    // Check Process Time
    std::cout << "Testing process time" << std::endl;
    try {
      std::cout << "  Process time implemented by: "
                << Dune::PDELab::getProcessTimeImp() << std::endl;
      std::cout << "  Process time resolution: "
                << Dune::PDELab::getProcessTimeResolution() << std::endl;
      std::cout << "  Current process time: " << Dune::PDELab::getProcessTime()
                << std::endl;
    }
    catch(Dune::PDELab::ClockError &e) {
      std::cerr << "  Cought ClockError: " << e << std::endl;
      result = 1;
    }

    // tests done
    return result;
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
