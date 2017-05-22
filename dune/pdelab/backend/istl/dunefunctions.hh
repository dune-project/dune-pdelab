// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_ISTL_DUNEFUNCTIONS_HH
#define DUNE_PDELAB_BACKEND_ISTL_DUNEFUNCTIONS_HH

#include <dune/pdelab/backend/istl.hh>

namespace Dune {
  namespace PDELab {
    namespace ISTL {

      template<std::size_t block_size = 1>
      struct SimpleVectorBackend
      {};

    // can't have the closing of the namespace inside the #ifndef DOXYGEN block
    } // namespace istl

#ifndef DOXYGEN

    namespace Backend {
      namespace impl {

        template<std::size_t block_size, typename GFS, typename E>
        struct BackendVectorSelectorHelper<ISTL::SimpleVectorBackend<block_size>, GFS, E>
        {
          using type = Dune::PDELab::ISTL::BlockVector<GFS,Dune::BlockVector<Dune::FieldVector<E,block_size>>>;
          using Type = type;
        };

      } // namespace impl
    } // namespace Backend

#endif // DOXYGEN

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_BACKEND_ISTL_VECTOR_HH
