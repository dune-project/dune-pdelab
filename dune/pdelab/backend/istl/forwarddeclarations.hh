// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_ISTL_FORWARDDECLARATIONS_HH
#define DUNE_PDELAB_BACKEND_ISTL_FORWARDDECLARATIONS_HH

#ifndef DOXYGEN // These forward declarations are of no concern to Doxygen

#include <dune/common/memory/domain.hh>

namespace Dune {

  // ********************************************************************************
  // forward declarations of tagged types to avoid including their headers
  // ********************************************************************************

  template<typename F, int n>
  class FieldVector;

  template<typename F, int n, int m>
  class FieldMatrix;

  template<typename F>
  class DynamicVector;

  template<typename F>
  class DynamicMatrix;

  template<typename Block, typename Alloc>
  class BlockVector;

  template<typename Block, typename Alloc>
  class BCRSMatrix;

  namespace ISTL {

    template<typename E, typename Alloc, typename D>
    class Vector;

#if !defined DUNE_ISTL_ELLMATRIX_DECLARED
    template<typename F_, typename A_, typename D_ = typename Memory::allocator_domain<A_>::type>
    class ELLMatrix;
#define DUNE_ISTL_ELLMATRIX_DECLARED 1
#endif

  }

  namespace PDELab {

    template<typename GFS, typename C>
    class ISTLBlockVectorContainer;

    template<typename GFSV, typename GFSU, typename C>
    class ISTLMatrixContainer;

    namespace istl {

      template<typename E, typename VV, typename VU>
      struct build_matrix_type;

      template<typename GFS, typename C>
      class FlatVectorContainer;

      template<typename GFSV, typename GFSU, typename C>
      class FlatELLMatrixContainer;

    } // namespace istl
  } // namespace PDELab
} // namespace Dune

#endif // DOXYGEN

#endif // DUNE_PDELAB_BACKEND_ISTL_FORWARDDECLARATIONS_HH
