// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_ISTL_FORWARDDECLARATIONS_HH
#define DUNE_PDELAB_BACKEND_ISTL_FORWARDDECLARATIONS_HH

#ifndef DOXYGEN // These forward declarations are of no concern to Doxygen

#include <dune/istl/forwarddeclarations.hh>


namespace Dune {

  // ********************************************************************************
  // forward declarations of tagged types to avoid including their headers
  // ********************************************************************************

  namespace PDELab {

    template<typename GFS, typename C>
    class ISTLBlockVectorContainer;

    template<typename GFSV, typename GFSU, typename C, typename Stats>
    class ISTLMatrixContainer;

    namespace istl {

      template<typename E, typename VV, typename VU>
      struct build_matrix_type;

      template<typename GFS, typename C>
      class FlatVectorContainer;

      template<typename GFSV, typename GFSU, typename C>
      class FlatELLMatrixContainer;

      template<typename GFS, typename C>
      class BlockVectorContainer;

      template<typename GFSV, typename GFSU, typename C>
      class BELLMatrixContainer;

    } // namespace istl
  } // namespace PDELab
} // namespace Dune

#endif // DOXYGEN

#endif // DUNE_PDELAB_BACKEND_ISTL_FORWARDDECLARATIONS_HH
