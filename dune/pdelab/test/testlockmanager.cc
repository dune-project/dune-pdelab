// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

// silence the nullmutex warning
#define SILENCE_NULLMUTEX_WARNING 1

#include <chrono>
#include <cstddef>
#include <iostream>
#include <mutex>
#include <ostream>
#include <thread>
#include <vector>

#include <dune/common/array.hh>
#include <dune/common/classname.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/pdelab/common/nullmutex.hh>
#include <dune/pdelab/common/lockmanager.hh>

template<class GV, class LM>
bool testLM(const GV &gv, LM &lm) {
  Dune::MultipleCodimMultipleGeomTypeMapper<GV, Dune::MCMGElementLayout>
    mapper(gv);
  std::vector<double> data(mapper.size(), 0);
  std::size_t iters = 100;
  auto exercise = [&] {
    auto end = gv.template end<0>();
    for(std::size_t i = 0; i < iters; ++i)
      for(auto it = gv.template begin<0>(); it != end; ++it) {
        auto &lock = lm[*it];
        std::lock_guard<decltype(lock)> guard(lock);
        volatile double value = data[mapper.map(*it)];
        std::this_thread::yield();
        std::this_thread::sleep_for(std::chrono::duration<double>(0.001));
        data[mapper.map(*it)] = value + 1;
      }
  };

  std::size_t nthreads = 10;
  std::thread threads[nthreads];
  for(std::size_t i = 0; i < nthreads; ++i)
    threads[i] = std::thread(exercise);
  for(std::size_t i = 0; i < nthreads; ++i)
    threads[i].join();
  for(auto val : data)
    if(val != nthreads * iters)
      return false;
  return true;
}

void pass(int &result) {
  std::cout << "Result: passed" << std::endl;
  if(result == 77)
    result = 0;
}

void fail(int &result) {
  std::cout << "Result: failed" << std::endl;
  result = 1;
}


int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper::instance(argc, argv);

    // need a grid in order to test grid functions
    Dune::FieldVector<double,2> L(1.0);
    Dune::array<int,2> N(Dune::fill_array<int,2>(1));
    Dune::YaspGrid<2> grid(L,N);
    grid.globalRefine(2);
    auto gv = grid.leafGridView();

    int result = 77;

    {
      Dune::PDELab::GlobalLockManager<std::mutex> lm;
      std::cout << "Checking " << Dune::className(lm) << std::endl;
      if(testLM(gv, lm))
        pass(result);
      else
        fail(result);
    }

    {
      Dune::PDELab::GlobalLockManager<Dune::PDELab::NullMutex> lm;
      std::cout << "Checking " << Dune::className(lm) << std::endl;
      if(testLM(gv, lm))
        std::cout << "Result: passed (but ignored)" << std::endl;
      else
        std::cout << "Result: failed (but ignored)" << std::endl;
    }

    {
      Dune::PDELab::PerElementLockManager<Dune::YaspGrid<2>::LeafGridView,
                                          std::mutex> lm(gv);
      std::cout << "Checking " << Dune::className(lm) << std::endl;
      if(testLM(gv, lm))
        pass(result);
      else
        fail(result);
    }

    // test passed
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
