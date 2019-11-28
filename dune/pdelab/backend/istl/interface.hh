// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_ISTL_INTERFACE_HH
#define DUNE_PDELAB_BACKEND_ISTL_INTERFACE_HH

namespace Dune {
  namespace PDELab {
    //! \addtogroup Backend
    //! \ingroup PDELab
    //! \{

    //=========================================================================
    // Base class providing a matrix based apply, a matrix free apply and a
    // method that returns if this solver is supposed to work matrix free. This
    // way we can avoid duplicating the StationaryLinearProblemSolver.
    //=========================================================================

    /** \brief Base class for ISTL backends
     *
     * Provide matrix based and matrix free apply methods.
     */
    class ISTLBackend_Base
    {
    public:
      //! Matrix based apply method
      template<class M, class V, class W>
      void apply(M& A, V& z, W& r, typename Dune::template FieldTraits<typename W::ElementType >::real_type reduction)
      {
        DUNE_THROW(Exception, "This solver does not work matrix based.");
      }

      //! Matrix free apply method
      template<class V, class W>
      void apply(V& z, W& r, typename Dune::template FieldTraits<typename W::ElementType >::real_type reduction)
      {
        DUNE_THROW(Exception, "This solver does not work matrix free.");
      }

      //! Return if this solver works matrix free
      bool matrixFree()
      {
        return false;
      }
    };

    //! \} group Backend
  }
}

#endif
