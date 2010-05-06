// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_DEFAULT_HH
#define DUNE_PDELAB_DEFAULT_HH

#include <dune/pdelab/gridoperatorspace/localmatrix.hh>

namespace Dune {
  namespace PDELab {
    //! \addtogroup LocalOperator
    //! \ingroup PDELab
    //! \{

    //! Default implementation for alpha volume
    /**
     * This class implements alpha_volume() using jacobian_volume().
     *
     * \tparam Imp Type of the class that implements jacobian_volume().  Imp
     *             must be derived from AlphaVoumeFromJacobianVolume<Imp>
     *             (curiously recurring template pattern, CRTP).
     */
    template<typename Imp>
    class AlphaVolumeFromJacobianVolume {
    public:
      //! compute residual on volume
      /**
       * \param eg   ElementGeometry object.
       * \param lfsu LocalFunctionSpace for x (this is the domain space of the
       *             local operator).
       * \param x    Coefficients vector to apply the operator to.
       * \param lfsv LocalFunctionSpace for r (this is the range space of the
       *             local operator).
       * \param r    Resulting residual.  alpha_volume() adds its residual to
       *             this parameter, so if r must be set to zero or similar it
       *             must be done before calling this function.  We access the
       *             first element of the vector to extract the type of the
       *             vectors elements, so the expression \c r[0] must be valid.
       */
      template<typename EG, typename LFSU, typename X,
               typename LFSV, typename R>
      void alpha_volume (const EG& eg,
                         const LFSU& lfsu, const X& x,
                         const LFSV& lfsv, R& r) const
      {
        // Trick to extract the type of the vector elements: call the helper
        // function with an additional dummy argument that has the type of the
        // vectors elements.  Although the argument is not actually used, this
        // allows the helper function to infer its type.  Of course it
        // requires that the expression r[0] is valid even if the value is not
        // actually used. 
        alpha_volume_helper(eg, lfsu, x, lfsv, r, r[0]);
      }

    private:
      template<typename EG, typename LFSU, typename X,
               typename LFSV, typename R, typename E>
      void alpha_volume_helper(const EG& eg,
                               const LFSU& lfsu, const X& x,
                               const LFSV& lfsv, R& r, const E&) const
      {
        LocalMatrix<E> mat(lfsv.size(), lfsu.size(), 0);
        static_cast<const Imp*>(this)->jacobian_volume(eg, lfsu, x, lfsv,
                                                       mat);
        for(unsigned i = 0; i < lfsv.size(); ++i)
          for(unsigned j = 0; j < lfsu.size(); ++j)
            r[i] += mat(i, j) * x[j];
      }
    };

    //! Default implementation for alpha boundary
    /**
     * This class implements alpha_boundary using jacobian boundary.
     *
     * \tparam Imp Type of the class that implements jacobian_boundary().  Imp
     *             must be derived from AlphaBoundaryFromJacobianBoundary<Imp>
     *             (curiously recurring template pattern, CRTP).
     */
    template<typename Imp>
    class AlphaBoundaryFromJacobianBoundary {
    public:
      //! compute residual on boundary
      /**
       * \tparam CV Type of residual.  We use the member typedef \c value_type
       *            to extract the type of the elements (should work with \c
       *            std::vector at least).
       *
       * \param ig   IntersectionGeometry object.
       * \param lfsu LocalFunctionSpace for x (this is the domain space of the
       *             local operator).
       * \param x    Coefficients vector to apply the operator to.
       * \param lfsv LocalFunctionSpace for r (this is the range space of the
       *             local operator).
       * \param r    Resulting residual.  alpha_boundary() adds its residual to
       *             this parameter, so if r must be set to zero or similar it
       *             must be done before calling this function.
       */
      template<typename IG, typename LFSU, typename X,
               typename LFSV, typename CV>
      void alpha_boundary(const IG& ig,
                          const LFSU& lfsu, const X& x,
                          const LFSV& lfsv, CV& r) const
      {
        LocalMatrix<typename CV::value_type> mat(lfsv.size(), lfsu.size(), 0);
        static_cast<const Imp*>(this)->jacobian_boundary(ig, lfsu, x, lfsv,
                                                         mat);
        for(unsigned i = 0; i < lfsv.size(); ++i)
          for(unsigned j = 0; j < lfsu.size(); ++j)
            r[i] += mat(i, j) * x[j];
      }
    };

    //! \} group LocalOperator
  }
}

#endif // DUNE_PDELAB_DEFAULT_HH
