// -*- tab-width: 4; indent-tabs-mode: nil -*-

#ifndef DOXYGEN

// this whole file is an implementation detail and should never be included directly!

#ifndef _DUNE_PDELAB_HAVE_LOWERCASE_ISTL_NAMESPACE
#define _DUNE_PDELAB_HAVE_LOWERCASE_ISTL_NAMESPACE

/** \brief For backward compatibility -- Do not use this! */
namespace Dune {
  namespace PDELab {
    namespace ISTL {}
    namespace istl {
      using namespace Dune::PDELab::ISTL;
    }
  }
}

#endif // _DUNE_PDELAB_HAVE_LOWERCASE_ISTL_NAMESPACE

#if not defined (DUNE_PDELAB_BACKEND_ISTL_HH) and not defined (_DUNE_PDELAB_SUPPRESS_ISTL_HH_WARNING)

#warning Directly including files from dune/pdelab/backend/istl/ is deprecated in PDELab 2.5.0 and will not be supported \
  after that release, please switch to only including dune/pdelab/backend/istl.hh instead

#endif // not defined (DUNE_PDELAB_BACKEND_ISTL_HH) and not defined (_DUNE_PDELAB_SUPPRESS_ISTL_HH_WARNING)

#endif // DOXYGEN
