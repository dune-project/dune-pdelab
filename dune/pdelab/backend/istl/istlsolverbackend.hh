// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_ISTL_ISTLSOLVERBACKEND_HH
#define DUNE_PDELAB_BACKEND_ISTL_ISTLSOLVERBACKEND_HH

#include "seqistlsolverbackend.hh"
#include "ovlpistlsolverbackend.hh"
#include "novlpistlsolverbackend.hh"

  /**
   * @brief For better handling istlsolverbackend.hh is now divided into:
   * parallelistlhelper.hh for creation of index sets (mainly for amg)
   * seqistlsolverbackend.hh for sequential solvers
   * ovlpistlsolverbackend.hh for overlapping solvers,operators,...
   * novlpistlsolverbackend.hh with nonoverlapping solvers,operators,...
   */

#endif // DUNE_PDELAB_BACKEND_ISTL_ISTLSOLVERBACKEND_HH
