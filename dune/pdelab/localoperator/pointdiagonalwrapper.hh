// -*- tab-width: 2; indent-tabs-mode: nil -*-
// vi: set et ts=2 sw=2 sts=2:
#ifndef DUNE_PDELAB_LOCALOPERATOR_POINTDIAGONALWRAPPER_HH
#define DUNE_PDELAB_LOCALOPERATOR_POINTDIAGONALWRAPPER_HH

namespace Dune {
  namespace PDELab {

    namespace impl {
      // This wraps a matrix accumulation view and only accumulates if diagonal
      // is set to true and if the indices i and j are equal. Can be used to
      // accumulate the point diagonal.
      template <typename AccumulationView>
      class PointDiagonalAccumulationViewWrapper
      {
      public:
        PointDiagonalAccumulationViewWrapper(AccumulationView& view, bool diagonal)
          : _view(view), _diagonal(diagonal)
        {}

        template <typename LFSU, typename I, typename LFSV, typename J, typename Value>
        void accumulate(const LFSU& lfsu, I i, const LFSV& lfsv, J j, Value value)
        {
          if (_diagonal && i == j){
            _view.accumulate(lfsu, i, value);
          }
        }

      private:
        AccumulationView& _view;
        bool _diagonal;
      };

    } // namespace impl


    /** \brief A local operator that accumulates the point diagonal
     *
     * Since the point diagonal has exactly the same size as a residual vector
     * we use the residual methods to implement it (instead of going through
     * the jacobian methods and wasting lots of memory).
     *
     * Note: This operator does skeletons in a one sided fashion. This means
     * that every skeleton-intersection is visited twice for assembling the
     * whole block diagonal. Once from one side and once from the other. This
     * behavior is needed for the implementation of matrix-free block
     * preconditioners.
     *
     * \tparam[in] LocalOperator Type of the local operator that gets wrapped
     */
    template <class LocalOperator>
    class PointDiagonalLocalOperatorWrapper
      : public Dune::PDELab::FullVolumePattern
      ,  public Dune::PDELab::LocalOperatorDefaultFlags
    {
    public:
      // We only want to assemble the point diagonal so we only need volume
      // pattern.
      static constexpr bool doPatternVolume = true;

      // For assembling the point diagonal correctly we might need to call
      // volume, skeleton and boundary
      static constexpr bool doAlphaVolume = LocalOperator::doAlphaVolume;
      static constexpr bool doAlphaSkeleton = LocalOperator::doAlphaSkeleton;
      static constexpr bool doAlphaBoundary = LocalOperator::doAlphaBoundary;

      // If the underlying lop is linear, this one will be linear too
      static constexpr bool isLinear = LocalOperator::isLinear;

      // This wrapper is designed to use two sided assembly. If it was just
      // about assembling the whole diagonal block matrix one sided assembly
      // would be more efficient. For the implementation of matrix-free
      // preconditioner we need to assemble a diagonal block locally for a
      // given cell so we need to assemble in a two sided fashion.
      static constexpr bool doSkeletonTwoSided = true;

      /** \brief Construct new instance of class
       *
       * \param[in] _lop Wrapped local operator instance
       */
      PointDiagonalLocalOperatorWrapper(const LocalOperator& localOperator)
        : _localOperator(localOperator)

      {}

      /** \brief Copy constructor */
      PointDiagonalLocalOperatorWrapper(const PointDiagonalLocalOperatorWrapper& other)
        : _localOperator(other._localOperator)
      {}

      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
        impl::PointDiagonalAccumulationViewWrapper<R> view(r, true);
        _localOperator.jacobian_volume(eg, lfsu, x, lfsv, view);
      }

      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_skeleton (const IG& ig,
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                           R& r_s, R& r_n) const
      {
        impl::PointDiagonalAccumulationViewWrapper<R> view_ss(r_s, true);
        impl::PointDiagonalAccumulationViewWrapper<R> view_other(r_s, false);
        _localOperator.jacobian_skeleton(ig,
                                         lfsu_s, x_s, lfsv_s,
                                         lfsu_n, x_n, lfsv_n,
                                         view_ss, view_other, view_other, view_other);
      }

      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_boundary (const IG& ig,
                         const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                         R& r_s) const
      {
        impl::PointDiagonalAccumulationViewWrapper<R> view(r_s, true);
        _localOperator.jacobian_boundary(ig, lfsu_s, x_s, lfsv_s, view);
      }

    private:
      /** \brief Wrapped original local operator */
      const LocalOperator& _localOperator;
    };

    /** \brief A function for assembling the point diagonal of a single block
     *
     * This assumes that the localoperator that is passed in does assemble in a
     * two sided fashion.
     */
    template <typename LocalOperator, typename EG, typename LFSU, typename X, typename LFSV, typename Y>
    void assembleLocalPointDiagonal(const LocalOperator& localOperator,
                                    const EG& eg,
                                    const LFSU& lfsu,
                                    const X& x,
                                    const LFSV& lfsv,
                                    Y& y)
    {
      // Assemble the volume part
      localOperator.alpha_volume(eg, lfsu, x, lfsv, y);

      // Iterate over the intersections
      auto entitySet = lfsu.gridFunctionSpace().entitySet();
      std::size_t intersectionIndex = 0;
      for (const auto& is : intersections(lfsu.gridFunctionSpace().gridView(), eg.entity()))
      {
        Dune::PDELab::IntersectionGeometry<
          typename std::remove_reference<
            typename std::remove_const<
              decltype(is)
              >::type
            >::type
          > ig(is, intersectionIndex++);
        auto intersectionData = classifyIntersection(entitySet, is);
        auto intersectionType = std::get<0>(intersectionData);

        // Assemble the intersection part
        switch (intersectionType){
        case Dune::PDELab::IntersectionType::skeleton:
          localOperator.alpha_skeleton(ig, lfsu, x, lfsv, lfsu, x, lfsv, y, y);
          break;
        case Dune::PDELab::IntersectionType::boundary:
          localOperator.alpha_boundary(ig, lfsu, x, lfsv, y);
          break;
        default:
          break;
        }
      }
    }

  } // namespace PDELab
} // namespace Dune

#endif
