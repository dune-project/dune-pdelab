// -*- tab-width: 2; indent-tabs-mode: nil -*-
// vi: set et ts=2 sw=2 sts=2:

#ifndef DUNE_PDELAB_LOCALOPERATOR_BLOCKDIAGONALWRAPPER_HH
#define DUNE_PDELAB_LOCALOPERATOR_BLOCKDIAGONALWRAPPER_HH

#include <dune/pdelab/common/intersectiontype.hh>
#include <dune/pdelab/localoperator/jacobianapplyhelper.hh>

namespace Dune {
  namespace PDELab {

    namespace impl{

      // This wraps an accumulation view and only accumulates if diagonal is
      // set to true
      template <typename AccumulationView>
      class BlockDiagonalAccumulationViewWrapper
      {
      public:

        BlockDiagonalAccumulationViewWrapper(AccumulationView& view, bool diagonal)
          : _view(view), _diagonal(diagonal)
        {}

        const auto weight()
        {
          return _view.weight();
        }

        auto data()
        {
          return _view.data();
        }

        const auto data() const
        {
          return _view.data();
        }

        template <typename LFSU, typename I, typename Value>
        void accumulate(const LFSU& lfsu, I i, Value value)
        {
          if (_diagonal){
            _view.accumulate(lfsu, i, value);
          }
        }

        template <typename LFSU, typename I, typename LFSV, typename J, typename Value>
        void accumulate(const LFSU& lfsu, I i, const LFSV& lfsv, J j, Value value)
        {
          if (_diagonal){
            _view.accumulate(lfsu, i, lfsv, j, value);
          }
        }

      private:
        AccumulationView& _view;
        bool _diagonal;
      };
    } // namespace impl


    /** \brief A local operator that accumulates the block diagonal
     *
     * This makes only sense for methods that have a block structure like
     * Discontinuous Galerking methods or Finite Vvolume methods. For those the
     * resulting operator assembles only the diagonal blocks when the jacobian
     * methods are called.
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
    class BlockDiagonalLocalOperatorWrapper
      : public Dune::PDELab::LocalOperatorDefaultFlags
    {
    public:
      // We only want to assemble the diagonal blocks so we only need volume
      // pattern.
      static constexpr bool doPatternVolume = true;

      // For assembling the diagonal block correctly we might need to call
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
       * \param[in] _localOperator Wrapped local operator instance
       */
      BlockDiagonalLocalOperatorWrapper(const LocalOperator& localOperator)
        : _localOperator(localOperator)
      {}

      /** \brief Copy constructor */
      BlockDiagonalLocalOperatorWrapper(const BlockDiagonalLocalOperatorWrapper& other)
        : _localOperator(other._localOperator)
      {}

      // define sparsity pattern of operator representation
      template<typename LFSU, typename LFSV, typename LocalPattern>
      void pattern_volume (const LFSU& lfsu, const LFSV& lfsv,
                           LocalPattern& pattern) const
      {
        _localOperator.pattern_volume(lfsu, lfsv, pattern);
      }

      template<typename EG, typename LFSU, typename X, typename LFSV, typename MAT>
      void jacobian_volume (const EG& eg,
                            const LFSU& lfsu,
                            const X& x,
                            const LFSV& lfsv,
                            MAT& mat) const
      {
        _localOperator.jacobian_volume(eg,lfsu,x,lfsv,mat);
      }

      template<typename IG, typename LFSU, typename X, typename LFSV, typename MAT>
      void jacobian_skeleton (const IG& ig,
                              const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                              const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                              MAT& mat_ss, MAT& mat_sn,
                              MAT& mat_ns, MAT& mat_nn) const
      {
        // In the end this jacobian_skeleton only accumulates the self-self
        // part which corresponds to a diagonal block. The localoperator uses
        // two sided assembly, so we can discard the neighbor-neighbor block
        // since it will be accumulated through another cell.
        impl::BlockDiagonalAccumulationViewWrapper<MAT> view_ss(mat_ss, true);
        impl::BlockDiagonalAccumulationViewWrapper<MAT> view_other(mat_ss, false);
        _localOperator.jacobian_skeleton(ig, lfsu_s, x_s, lfsv_s, lfsu_n, x_n, lfsv_n, view_ss, view_other, view_other, view_other);
      }

      template<typename IG, typename LFSU, typename X, typename LFSV, typename MAT>
      void jacobian_boundary (const IG& ig,
                              const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                              MAT& mat_ss) const
      {
        _localOperator.jacobian_boundary(ig, lfsu_s, x_s, lfsv_s, mat_ss);
      }


      template<typename EG, typename LFSU, typename X, typename LFSV, typename Y>
      void jacobian_apply_volume(const EG& eg, const LFSU& lfsu, const X& z, const LFSV& lfsv, Y& y) const
      {
        Dune::PDELab::impl::jacobianApplyVolume(_localOperator, eg, lfsu, z, lfsv, y);
      }
      template<typename EG, typename LFSU, typename X, typename Z, typename LFSV, typename Y>
      void jacobian_apply_volume(const EG& eg, const LFSU& lfsu, const X& x, const Z& z, const LFSV& lfsv, Y& y) const
      {
        Dune::PDELab::impl::jacobianApplyVolume(_localOperator, eg, lfsu, x, z, lfsv, y);
      }


      template<typename IG, typename LFSU, typename Z, typename LFSV, typename Y>
      void jacobian_apply_skeleton(const IG& ig,
                                    const LFSU& lfsu_s, const Z& z_s, const LFSV& lfsv_s,
                                    const LFSU& lfsu_n, const Z& z_n, const LFSV& lfsv_n,
                                    Y& y_s, Y& y_n) const
      {
        impl::BlockDiagonalAccumulationViewWrapper<Y> view_s(y_s, true);
        impl::BlockDiagonalAccumulationViewWrapper<Y> view_n(y_n, false);
        // Note for the application of only the diagonal block we need to pass
        // a coefficient vector 0 for the test function into the original
        // lop. Otherwise it will also apply one of the off-diagonal blocks.
        Z z_zero(z_n.size(), 0.0);
        Dune::PDELab::impl::jacobianApplySkeleton(_localOperator, ig, lfsu_s, z_s, lfsv_s, lfsu_n, z_zero, lfsv_n, view_s, view_n);
      }
      template<typename IG, typename LFSU, typename X, typename Z, typename LFSV, typename Y>
      void jacobian_apply_skeleton(const IG& ig,
                                   const LFSU& lfsu_s, const X& x_s, const Z& z_s, const LFSV& lfsv_s,
                                   const LFSU& lfsu_n, const X& x_n, const Z& z_n, const LFSV& lfsv_n,
                                   Y& y_s, Y& y_n) const
      {
        impl::BlockDiagonalAccumulationViewWrapper<Y> view_s(y_s, true);
        impl::BlockDiagonalAccumulationViewWrapper<Y> view_n(y_n, false);
        // Note for the application of only the diagonal block we need to pass
        // a coefficient vector 0 for the test function into the original
        // lop. Otherwise it will also apply one of the off-diagonal blocks.
        Z z_zero(z_n.size(), 0.0);
        Dune::PDELab::impl::jacobianApplySkeleton(_localOperator, ig, lfsu_s, x_s, z_s, lfsv_s, lfsu_n, x_n, z_zero, lfsv_n, view_s, view_n);
      }

      template<typename IG, typename LFSU, typename X, typename LFSV, typename Y>
      void jacobian_apply_boundary (const IG& ig,
                                    const LFSU& lfsu_s, const X& z_s, const LFSV& lfsv_s,
                                    Y& y_s) const
      {
        Dune::PDELab::impl::jacobianApplyBoundary(_localOperator, ig, lfsu_s, z_s, lfsv_s, y_s);
      }

      template<typename IG, typename LFSU, typename X, typename Z, typename LFSV, typename Y>
      void jacobian_apply_boundary (const IG& ig,
                                    const LFSU& lfsu_s, const X& x_s, const Z& z_s, const LFSV& lfsv_s,
                                    Y& y_s) const
      {
        Dune::PDELab::impl::jacobianApplyBoundary(_localOperator, ig, lfsu_s, x_s, z_s, lfsv_s, y_s);
      }


    private:
      /** \brief Wrapped original local operator */
      const LocalOperator& _localOperator;
    };


    //! A function for assembling a single diagonal block
    template <typename LocalOperator, typename EG, typename LFSU, typename X, typename LFSV, typename MAT>
    void assembleLocalDiagonalBlock(const LocalOperator& localOperator,
                                    const EG& eg,
                                    const LFSU& lfsu,
                                    const X& x,
                                    const LFSV& lfsv,
                                    MAT& mat)
    {
       // Assemble the volume part
      localOperator.jacobian_volume(eg, lfsu, x, lfsv, mat);

      // Iterate over the intersections
      auto entitySet = lfsu.gridFunctionSpace().entitySet();
      std::size_t intersectionIndex = 0;
      for (const auto& is : intersections(lfsu.gridFunctionSpace().gridView(), eg.entity()))
      {
        Dune::PDELab::IntersectionGeometry<std::decay_t<decltype(is)>> ig(is, intersectionIndex++);
        auto intersectionData = classifyIntersection(entitySet, is);
        auto intersectionType = std::get<0>(intersectionData);

        // Assemble the intersection part
        switch (intersectionType){
        case Dune::PDELab::IntersectionType::skeleton:
          localOperator.jacobian_skeleton(ig, lfsu, x, lfsv, lfsu, x, lfsv, mat, mat, mat, mat);
          break;
        case Dune::PDELab::IntersectionType::boundary:
          localOperator.jacobian_boundary(ig, lfsu, x, lfsv, mat);
          break;
        default:
          break;
        }
      }
    }

    //! A function for applying a single diagonal block
    template<typename LocalOperator, typename EG, typename LFSU, typename X, typename Z, typename LFSV, typename Y>
    void applyLocalDiagonalBlock(const LocalOperator& localOperator,
                                 const EG& eg,
                                 const LFSU& lfsu,
                                 const X& x,
                                 const Z& z,
                                 const LFSV& lfsv,
                                 Y& y)
    {
      // Apply the volume part
      if (LocalOperator::isLinear)
        Dune::PDELab::impl::jacobianApplyVolume(localOperator, eg, lfsu, z, lfsv, y);
      else
        Dune::PDELab::impl::jacobianApplyVolume(localOperator, eg, lfsu, x, z, lfsv, y);

      // Iterate over the intersections
      auto entitySet = lfsu.gridFunctionSpace().entitySet();
      std::size_t intersectionIndex = 0;
      for (const auto& is : intersections(lfsu.gridFunctionSpace().gridView(), eg.entity()))
      {
        Dune::PDELab::IntersectionGeometry<std::decay_t<decltype(is)>> ig(is, intersectionIndex++);
        auto intersectionData = classifyIntersection(entitySet, is);
        auto intersectionType = std::get<0>(intersectionData);

        // Assemble the intersection part
        switch (intersectionType){
        case Dune::PDELab::IntersectionType::skeleton:
          if (LocalOperator::isLinear)
            Dune::PDELab::impl::jacobianApplySkeleton(localOperator, ig, lfsu, z, lfsv, lfsu, z, lfsv, y, y);
          else
            Dune::PDELab::impl::jacobianApplySkeleton(localOperator, ig, lfsu, x, z, lfsv, lfsu, x, z, lfsv, y, y);
          break;
        case Dune::PDELab::IntersectionType::boundary:
          if (LocalOperator::isLinear)
            Dune::PDELab::impl::jacobianApplyBoundary(localOperator, ig, lfsu, z, lfsv, y);
          else
            Dune::PDELab::impl::jacobianApplyBoundary(localOperator, ig, lfsu, x, z, lfsv, y);
          break;
        default:
          break;
        }
      }
    }

  } // namespace PDELab
} // namespace Dune

#endif
