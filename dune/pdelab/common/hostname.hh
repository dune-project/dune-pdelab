// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_COMMON_HOSTNAME_HH
#define DUNE_PDELAB_COMMON_HOSTNAME_HH

#include <string>

namespace Dune {
  namespace PDELab {

    //! C++ friendly wrapper around POSIX' gethostname()
    /**
     * \note Anything after the first '.' is stripped from the value returned
     *       by the underlying implementation.
     */
    std::string getHostName();

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_COMMON_HOSTNAME_HH
