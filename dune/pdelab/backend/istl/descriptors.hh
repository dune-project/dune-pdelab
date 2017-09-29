// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_ISTL_DESCRIPTORS_HH
#define DUNE_PDELAB_BACKEND_ISTL_DESCRIPTORS_HH

// this is here for backwards compatibility and deprecation warnings, remove after 2.5.0
#include "ensureistlinclude.hh"

#include <dune/pdelab/backend/interface.hh>
#include <dune/pdelab/backend/istl/forwarddeclarations.hh>
#include <dune/pdelab/backend/istl/matrixhelpers.hh>
#include <dune/pdelab/backend/istl/utility.hh>
#include <cstddef>

namespace Dune {
  namespace PDELab {

#ifndef DOXYGEN
      template<typename T>
      constexpr bool deactivate_standard_blocking_for_ordering(const T&)
      {
        return false;
      }
#endif

    namespace ISTL {

      //! The type of blocking employed at this node in the function space tree.
      enum class Blocking
      {
        //! No blocking at this level.
        none,
        //! Creates one macro block for each child space, each block is a BlockVector / BCRS matrix.
        bcrs,
        //! Create fixed size blocks that each group together a fixed number of DOFs from each child space.
        /**
         * Creates a block structure with fixed size child blocks that each group together a fixed number
         * of DOFs from each child space. Typically used for entity-wise blocking of DOFs across child spaces
         * for e.g. velocity in flow problems or concentrations in multi-component transport.
         *
         * \note This type of blocking cannot be nested due to limitations in ISTL.
         */
        fixed,
      };

      //! Tag describing an ISTL BlockVector backend.
      struct vector_backend_tag {};

      template<Blocking blocking = Blocking::none, std::size_t block_size_ = 0>
      struct VectorBackend
      {

        using tag = vector_backend_tag;

        using size_type = std::size_t;

        static const size_type blockSize = block_size_;

        struct Traits
        {
          static const Blocking block_type = blocking;
          static const size_type block_size = block_size_;

          static const bool blocked = blocking != Blocking::none;

          static const size_type max_blocking_depth = blocked ? 1 : 0;
        };

        template<typename GFS>
        bool blocked(const GFS& gfs) const
        {
          if (deactivate_standard_blocking_for_ordering(gfs.orderingTag()))
            return false;
          // We have to make an exception for static blocking and block_size == 1:
          // In that case, the ISTL backends expect the redundant index information
          // at that level to be elided, and keeping it in here will break vector
          // and matrix accesses.
          // To work around that problem, we override the user and just turn off
          // blocking internally.
          return Traits::blocked && (blocking != Blocking::fixed || !GFS::isLeaf || block_size_ > 1);
        }

      };

    }

  } // namespace PDELab
} // namespace Dune


#endif // DUNE_PDELAB_BACKEND_ISTL_DESCRIPTORS_HH
