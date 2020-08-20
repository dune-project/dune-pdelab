// -*- tab-width: 2; indent-tabs-mode: nil -*-
// vi: set et ts=2 sw=2 sts=2:
#ifndef DUNE_PDELAB_LOCALOPERATOR_BLOCKOFFDIAGONALWRAPPER_HH
#define DUNE_PDELAB_LOCALOPERATOR_BLOCKOFFDIAGONALWRAPPER_HH

#include <dune/pdelab/common/intersectiontype.hh>
#include <dune/pdelab/localoperator/blockdiagonalwrapper.hh>
namespace Dune {
  namespace PDELab {

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
    template <class LocalOperator>
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

        impl::BlockDiagonalAccumulationViewWrapper<Y> view_s_on(y_s, true);
        impl::BlockDiagonalAccumulationViewWrapper<Y> view_n_off(y_n, false);
        Z z_s_zero(z_s.size(), 0.0);
        Dune::PDELab::impl::jacobianApplySkeleton(_localOperator, ig, lfsu_s, z_s_zero, lfsv_s, lfsu_n, z_n, lfsv_n, view_s_on, view_n_off);
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

        impl::BlockDiagonalAccumulationViewWrapper<Y> view_s_on(y_s, true);
        impl::BlockDiagonalAccumulationViewWrapper<Y> view_n_off(y_n, false);
        Z z_s_zero(z_s.size(), 0.0);
        Dune::PDELab::impl::jacobianApplySkeleton(_localOperator, ig, lfsu_s, x_s, z_s_zero, lfsv_s, lfsu_n, x_n, z_n, lfsv_n, view_s_on, view_n_off);
      }

    private:
      /** \brief Wrapped original local operator */
      const LocalOperator& _localOperator;
    };

  } // namespace PDELab
} // namespace Dune

#endif
