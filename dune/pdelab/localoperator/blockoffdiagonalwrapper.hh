// -*- tab-width: 2; indent-tabs-mode: nil -*-
// vi: set et ts=2 sw=2 sts=2:
#ifndef DUNE_PDELAB_LOCALOPERATOR_BLOCKOFFDIAGONALWRAPPER_HH
#define DUNE_PDELAB_LOCALOPERATOR_BLOCKOFFDIAGONALWRAPPER_HH

#include <dune/pdelab/common/intersectiontype.hh>
#include <dune/pdelab/localoperator/blockdiagonalwrapper.hh>
namespace Dune {
  namespace PDELab {

    namespace impl {

      // This can be used to get a vector view that returns a zero coefficient.
      template <typename View>
      class ZeroViewWrapper
      {
      public:
        using Container = typename View::Container;
        using ElementType = typename View::value_type;
        using SizeType = typename View::size_type;

        // If zero is set to false this class will forward all container
        // accesses to the vector view that is passed as first argument. This
        // means it will basically behave the same way as the view itself.
        //
        // If you set zero to false it will return 0.0 when it is asked for a
        // coefficient.
        //
        // Use case: The methods in the localoperator interface sometimes get
        // multiple coefficient views that need to have the same type (e.g. x_s
        // and x_n for the ansatz function in skeleton terms). This can be used
        // to 'null' one of those vectors without actually changing any values
        // in memory.
        ZeroViewWrapper(const View& view, bool zero)
          : _view(view), _zero(zero), _zeroCoefficient(0.0)
        {}

        template <typename LFS>
        const ElementType& operator()(const LFS& lfs, SizeType i) const
        {
          if (_zero)
            return _zeroCoefficient;
          else
            return _view.container()(lfs, i);
        }

      private:
        const View& _view;
        bool _zero;
        ElementType _zeroCoefficient;
      };

      // Interfaces look different in the fast-DG case
      template <typename Container, typename LocalFunctionSpaceCache>
      class ZeroViewWrapper<AliasedVectorView<Container, LocalFunctionSpaceCache>>
      {
      public:
        using View = ConstAliasedVectorView<Container, LocalFunctionSpaceCache>;
        using ElementType = typename View::ElementType;
        using SizeType = typename View::size_type;

        ZeroViewWrapper(const View& view, bool zero)
          : _view(view), _zero(zero), _zeroCoefficient(0.0)
        {}

        template <typename LFS>
        const ElementType& operator()(const LFS& lfs, SizeType i) const
        {
          if (_zero)
            return _zeroCoefficient;
          else
            return _view(lfs, i);
        }

        const ElementType* data() const
        {
          // Note: There is no principal problem in implementing this. Create
          // an array of ElementType that has the correct size and contains
          // only zeros. This was not implemted since there was no way of
          // testing the implementation. Better to have a clear error message
          // than a delicate implementation bug.
          DUNE_THROW(Dune::Exception, "So far the ZeroViewWrapper does not support fast DG local operators using the data() method to access coefficients. .");
        }

      private:
        const View& _view;
        bool _zero;
        ElementType _zeroCoefficient;
      };


    } // namespace impl

    /** \brief A local operator that accumulates the off block diagonal
     *
     * This makes only sense for methods that have a block structure like
     * Discontinuous Galerking methods or Finite Vvolume methods. For those the
     * resulting operator assembles only the off diagonal blocks when the
     * jacobian methods are called.
     *
     * Note: This operator does skeletons in a one sided fashion. This means
     * that every skeleton-intersection is visited twice for assembling the
     * whole block diagonal. Once from one side and once from the other. This
     * behavior is needed for the implementation of matrix-free block
     * preconditioners.
     *
     * \tparam[in] LocalOperator Type of the local operator that gets wrapped
     */
    template <typename LocalOperator>
    class BlockOffDiagonalLocalOperatorWrapper
      : public Dune::PDELab::LocalOperatorDefaultFlags
    {
    public:
      // We only want to assemble the off-diagonal blocks so we only need
      // skeleton pattern
      static constexpr bool doPatternSkeleton = true;

      // For assembling the off-diagonal blocks we only need alpha skeleton
      static constexpr bool doAlphaSkeleton = LocalOperator::doAlphaSkeleton;

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
      BlockOffDiagonalLocalOperatorWrapper(const LocalOperator& localOperator)
        : _localOperator(localOperator)
      {}

      /** \brief Copy constructor */
      BlockOffDiagonalLocalOperatorWrapper(const BlockOffDiagonalLocalOperatorWrapper& other)
        : _localOperator(other._localOperator)
      {}

      // define sparsity pattern connecting self and neighbor dofs
      template<typename LFSU, typename LFSV, typename LocalPattern>
      void pattern_skeleton (const LFSU& lfsu_s, const LFSV& lfsv_s, const LFSU& lfsu_n, const LFSV& lfsv_n,
                             LocalPattern& pattern_sn,
                             LocalPattern& pattern_ns) const
      {
        _localOperator.pattern_skeleton (lfsu_s, lfsv_s, lfsu_n, lfsv_n, pattern_sn, pattern_ns);
      }

      template<typename IG, typename LFSU, typename X, typename LFSV, typename MAT>
      void jacobian_skeleton (const IG& ig,
                              const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                              const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                              MAT& mat_ss, MAT& mat_sn,
                              MAT& mat_ns, MAT& mat_nn) const
      {
        // Since we do two sided assembly (visiting each intersection twice) we
        // have options. We could assemble either mat_sn or mat_ns. Both lead
        // to the same solution. In order to be consistent with the choice for
        // jacobian_apply_skeleton we will assemble m_ns.
        impl::BlockDiagonalAccumulationViewWrapper<MAT> view_ns(mat_ns, true);
        impl::BlockDiagonalAccumulationViewWrapper<MAT> view_other(mat_ns, false);
        _localOperator.jacobian_skeleton(ig, lfsu_s, x_s, lfsv_s, lfsu_n, x_n, lfsv_n, view_other, view_other, view_ns, view_other);
      }

      template<typename IG, typename LFSU, typename Z, typename LFSV, typename Y>
      void jacobian_apply_skeleton(const IG& ig,
                                   const LFSU& lfsu_s, const Z& z_s, const LFSV& lfsv_s,
                                   const LFSU& lfsu_n, const Z& z_n, const LFSV& lfsv_n,
                                   Y& y_s, Y& y_n) const
      {
        // This is more tricky. For a full Jacobian apply you would do:
        //
        // y_s = A_ss z_s + A_ns z_n
        // y_n = A_sn z_s + A_nn z_n
        //
        // Instead we want:
        //
        // (1) y_s = A_ns z_n
        // (2) y_n = A_sn z_s
        //
        // Since we do two sided assembly we only need to assemble one of these
        // two. For the matrix free preconditioner we need to apply the
        // jacobian locally for one block so we need to implement equation (1)
        // to get all contributions of other cell on the current one.

        // Set input coefficients z_s to zero
        impl::ZeroViewWrapper<Z> z_zero(z_s, true);
        impl::ZeroViewWrapper<Z> z_neigh(z_n, false);

        // Only accumulate in y_s
        impl::BlockDiagonalAccumulationViewWrapper<Y> view_s_on(y_s, true);
        impl::BlockDiagonalAccumulationViewWrapper<Y> view_n_off(y_n, false);

        // Apply Jacobian
        Dune::PDELab::impl::jacobianApplySkeleton(_localOperator, ig, lfsu_s, z_zero, lfsv_s, lfsu_n, z_neigh, lfsv_n, view_s_on, view_n_off);
      }

      template<typename IG, typename LFSU, typename X, typename Z, typename LFSV, typename Y>
      void jacobian_apply_skeleton(const IG& ig,
                                   const LFSU& lfsu_s, const X& x_s, const Z& z_s, const LFSV& lfsv_s,
                                   const LFSU& lfsu_n, const X& x_n, const Z& z_n, const LFSV& lfsv_n,
                                    Y& y_s, Y& y_n) const
      {
        // This is more tricky. For a full Jacobian apply you would do:
        //
        // y_s = A_ss z_s + A_ns z_n
        // y_n = A_sn z_s + A_nn z_n
        //
        // Instead we want:
        //
        // (1) y_s = A_ns z_n
        // (2) y_n = A_sn z_s
        //
        // Since we do two sided assembly we only need to assemble one of these
        // two. For the matrix free preconditioner we need to apply the
        // jacobian locally for one block so we need to implement equation (1)
        // to get all contributions of other cell on the current one.

        // Set input coefficients z_s to zero
        impl::ZeroViewWrapper<Z> z_zero(z_s, true);
        impl::ZeroViewWrapper<Z> z_neigh(z_n, false);

        // Only accumulate in y_s
        impl::BlockDiagonalAccumulationViewWrapper<Y> view_s_on(y_s, true);
        impl::BlockDiagonalAccumulationViewWrapper<Y> view_n_off(y_n, false);

        // Apply Jacobian
        Dune::PDELab::impl::jacobianApplySkeleton(_localOperator, ig, lfsu_s, x_s, z_zero, lfsv_s, lfsu_n, x_n, z_neigh, lfsv_n, view_s_on, view_n_off);
      }

    private:
      /** \brief Wrapped original local operator */
      const LocalOperator& _localOperator;
    };

  } // namespace PDELab
} // namespace Dune

#endif
