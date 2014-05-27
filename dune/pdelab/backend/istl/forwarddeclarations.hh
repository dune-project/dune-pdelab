// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_ISTL_FORWARDDECLARATIONS_HH
#define DUNE_PDELAB_BACKEND_ISTL_FORWARDDECLARATIONS_HH

#ifndef DOXYGEN // These forward declarations are of no concern to Doxygen

#include <dune/common/version.hh>

namespace Dune {

  // ********************************************************************************
  // forward declarations of tagged types to avoid including their headers
  // ********************************************************************************

  template<typename F, int n>
  class FieldVector;

  template<typename F, int n, int m>
  class FieldMatrix;

  // DynamicVector grew allocator support some time after the 2.3 release,
  // so we have to adjust the forward declaration accordingly

#if DUNE_VERSION_NEWER(DUNE_COMMON,2,4)

  template<typename F, typename Allocator>
  class DynamicVector;

#else

  template<typename F>
  class DynamicVector;

#endif

  template<typename F>
  class DynamicMatrix;

  template<typename Block, typename Alloc>
  class BlockVector;

  template<typename Block, typename Alloc>
  class BCRSMatrix;

  namespace PDELab {

    template<typename GFS, typename C>
    class ISTLBlockVectorContainer;

    template<typename GFSV, typename GFSU, typename C, typename Stats>
    class ISTLMatrixContainer;

    namespace istl {

      template<typename E, typename VV, typename VU>
      struct build_matrix_type;

    } // namespace istl
  } // namespace PDELab
} // namespace Dune

#endif // DOXYGEN

#endif // DUNE_PDELAB_BACKEND_ISTL_FORWARDDECLARATIONS_HH
