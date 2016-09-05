// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_SIMPLE_HH
#define DUNE_PDELAB_BACKEND_SIMPLE_HH

#include <dune/pdelab/backend/simple/descriptors.hh>
#include <dune/pdelab/backend/simple/vector.hh>
#include <dune/pdelab/backend/simple/matrix.hh>
#include <dune/pdelab/backend/simple/sparse.hh>

/** \brief For backward compatibility -- Do not use this! */
namespace Dune {
  namespace PDELab {
    namespace simple {
      using namespace Simple;
    }
  }
}

#endif // DUNE_PDELAB_BACKEND_SIMPLE_HH
