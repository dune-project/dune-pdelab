// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_ISTL_DUNEFUNCTIONS_HH
#define DUNE_PDELAB_BACKEND_ISTL_DUNEFUNCTIONS_HH

#include <dune/pdelab/backend/istl/vector.hh>

namespace Dune {
  namespace PDELab {
    namespace istl {

      template<typename V>
      struct SimpleVectorBackend
      {};

    // can't have the closing of the namespace inside the #ifndef DOXYGEN block
    } // namespace istl

#ifndef DOXYGEN

    namespace Backend {
      namespace impl {

        template<typename V, typename GFS, typename E>
        struct BackendVectorSelectorHelper<istl::SimpleVectorBackend<V>, GFS, E>
        {
          using type = Dune::PDELab::istl::BlockVector<GFS,V>;
          using Type = type;
        };

      } // namespace impl
    } // namespace Backend

#endif // DOXYGEN

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_BACKEND_ISTL_VECTOR_HH
