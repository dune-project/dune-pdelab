// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#if HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include <mutex>
#include <ostream>

#include "nullmutex.hh"

namespace Dune {
  namespace PDELab {

    NullMutex::NullMutex()
    {
      static std::once_flag flag;
      std::call_once(flag, [] {
          std::cerr << "Warning: Dune::PDELab::NullMutex used.  I hope this "
                    << "is just for debugging and benchmarking, any results "
                    << "will probably be bogus." << std::endl;
        });
    }

  } // namespace PDELab
} // namespace Dune
