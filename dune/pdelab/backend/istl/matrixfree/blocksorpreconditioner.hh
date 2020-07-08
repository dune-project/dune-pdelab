// -*- tab-width: 2; indent-tabs-mode: nil -*-
// vi: set et ts=2 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_ISTL_MATRIXFREE_BLOCKSORPRECONDITIONER_HH
#define DUNE_PDELAB_BACKEND_ISTL_MATRIXFREE_BLOCKSORPRECONDITIONER_HH

namespace Dune {
  namespace PDELab {

    /** \brief Turn a matrix-free Jacobi-type local preconditioner to a SOR one
     *
     * This preconditioner assumes that we have a discretization that leads to
     * a block structure, e.g. DG methods. It can be used to do a matrix-free
     * block-SOR preconditioner step.
     *
     * Given a linear system of equations Ax=b, a preconditioner step solves
     * Wv=d for a given vector d and an approximation $W \approx A$.
     *
     * Using the block decomposition $A=D+L+U$ block-SOR can be implemented in
     * the following way (approximately solving Av=d).
     *
     * for element T_i, i=1,...,m do:
     *   (1) a_i = d_i - \sum_{j<i} A_{ij} v_j^{(k)} - \sum_{j>i} A_{ij} v_j^{(k-1)}
     *   (2) Solve D_i b_i = a_i
     *   (3) Update v_i^{(k)} = (1-\omega) v_i^{(k-1)} + \omega b_i
     *
     * Here d_i, a_i, ... denote a local vector from the global block vector
     * and A_{ij} is the block (i,j) of the global block matrix A.
     *
     * See the artice "Matrix-free multigrid block-preconditioners for higher
     * order discontinuous Galerkin discretisations" by P. Bastian, E. Mueller,
     * S. Muething and M. Piatkowski.
     *
     * \tparam JacobianLOP Type of the Jacobi preconditioner local operator
     * \tparam BlockOffDiagonalLOP Type of the local operator for assembling the block off diagonal
     * \tparam GridFunctionSpace The type of grid function space.
     */
    template<typename JacobianLOP, typename BlockOffDiagonalLOP, typename GridFunctionSpace>
    class BlockSORPreconditionerLocalOperator
      : public Dune::PDELab::FullVolumePattern
      , public Dune::PDELab::LocalOperatorDefaultFlags
    {
      using value_type = typename GridFunctionSpace::Traits::GridView::Grid::ctype;
      using LocalVector = Dune::PDELab::LocalVector<value_type>;
      using W = typename LocalVector::WeightedAccumulationView;

    public :
      // A SOR preconditioner includes non-diagonal blocks so we need volume
      // and skeleton methods. We do two-sided assembly!
      static constexpr bool doPatternVolume = true;
      static constexpr bool doPatternSkeleton = true;
      static constexpr bool doPatternVolumePostSkeleton = true;
      static constexpr bool doAlphaVolume = true;
      static constexpr bool doAlphaSkeleton = true;
      static constexpr bool doAlphaVolumePostSkeleton = true;
      static constexpr bool doSkeletonTwoSided = true;
      static constexpr bool isLinear = JacobianLOP::isLinear;

      BlockSORPreconditionerLocalOperator(JacobianLOP& jacobianlop,
                                          BlockOffDiagonalLOP& blockOffDiagonalLOP,
                                          const GridFunctionSpace& gridFunctionSpace,
                                          const double omega=1.0)
        : _jacobianlop(jacobianlop)
        , _blockOffDiagonalLOP(blockOffDiagonalLOP)
        , _omega(omega)
        , _a_i(gridFunctionSpace.ordering().maxLocalSize())
        , _b_i(gridFunctionSpace.ordering().maxLocalSize())
      {}

      /** We use a Jacobi preconditioner that requires a setup. The setup will
       * be done in the alpha-volume method and later be used during the apply
       * methods.
       */
      bool requireSetup()
      {
        return _jacobianlop.requireSetup();
      }
      void setRequireSetup(bool v)
      {
        _jacobianlop.setRequireSetup(v);
      }

      //! Prepare underlying diagonal block preconditioner
      template<typename EG, typename LFSU, typename X, typename LFSV, typename Y>
      void alpha_volume(const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, Y& y) const
      {
        _jacobianlop.alpha_volume(eg,lfsu,x,lfsv,y);
      }

      //! Provide this method, but it actually does not nothing
      template<typename IG, typename LFSU, typename X, typename LFSV, typename Y>
      void alpha_skeleton(const IG& ig,
                          const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                          const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                          Y& y_s, Y& y_n) const
      {}

      //! Provide this method, but it actually does nothing
      template<typename EG, typename LFSU, typename X, typename LFSV, typename Y>
      void alpha_volume_post_skeleton(const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, Y& y) const
      {}

      //! Linear operator application, volume terms
      template<typename EG, typename LFSU, typename X, typename LFSV, typename Y>
      void jacobian_apply_volume(const EG& eg,
                                 const LFSU& lfsu, const X& x,
                                 const LFSV& lfsv, Y& y) const
      {
        //====================================
        // How does this preconditioner work
        //====================================
        //
        // When jacobian_apply is called on the grid operator the assembly
        // machine will iterate over the grid. Since we have set twosided to
        // true the order of methods will be:
        //
        // jacobian_apply_volume
        // jacobian_apply_skeleton for each intersection
        // jacobian_apply_volume_post_skeleton
        //
        // In the volume part we set _a_i to zero. During the skeleton part we
        // apply the off-diagonal blocks to the old and new solution, step (1)
        // of the algorithm description at the top of this class. In the
        // post_skeleton section we will do steps (2) and (3) of the algorithm.
        _a_i = 0.0;
      }

      //! linearized operator application, volume terms
      template<typename EG, typename LFSU, typename X, typename LFSV, typename Y>
      void jacobian_apply_volume(const EG& eg, const LFSU& lfsu, const X& x, const X& z, const LFSV& lfsv, Y& y) const
      {
        _a_i = 0.0;
      }

      //! Gather off-block-diagonals in Gauss-Seidel process of linear operator
      template<typename IG, typename LFSU, typename Z, typename LFSV, typename Y>
      void jacobian_apply_skeleton(const IG& ig,
                                   const LFSU& lfsu_s, const Z& z_s, const LFSV& lfsv_s,
                                   const LFSU& lfsu_n, const Z& z_n, const LFSV& lfsv_n,
                                   Y& y_s, Y& y_n) const
      {
        // This step is bit tricky.
        //
        // In the apply method of the GridOperatorPreconditioner class
        // jacobian_apply(d, v) is called on the preconditioner grid
        // operator.
        //
        // z_s and z_n will be the local entries of the block vector d.
        //
        // y_s will always be the solution of the previous step, e.i
        // $v_i^{(k-1)} in the above algorithm. This will be zero since that is
        // the value you get from ISTL for v^{(0)} and we do only one
        // preconditioner step.
        //
        // For cells you have not yet visited y_n will also be the old solution
        // $v_j^{(k-1)}. If you have already visited the neigbor cell it will
        // instead be the new solution $v_j^{(k)}$.
        //
        // For this reason we pass the local vectors of y_s and y_n as
        // coefficients to the block off-diagonal local operator. The result
        // will be accumulated in _a_i. Following the above algorithm we will
        // have
        //
        // _a_i = \sum_{j<i} A_ij v_j^{(k)} + \sum_{j>i} A_ij v_j^{(k-1)}
        //
        // after visiting all the intersections.
        W _a_i_view(_a_i, y_s.weight());
        if (ig.inside().partitionType() == Dune::InteriorEntity){
          // Note: Since the block off diagonal local operator works two sided
          // we can pass _a_i_view two times to the jacobian_apply_skeleton but
          // will only accumulate once.

          // TODO: Only works for FastDGGridOperator (y_* being an AliasedVectorView)
          _blockOffDiagonalLOP.jacobian_apply_skeleton(ig, lfsu_s, y_s, lfsv_s, lfsu_n, y_n, lfsv_s, _a_i_view, _a_i_view);
        }
      }

      //! Gather off-block-diagonals in Gauss-Seidel process of linearized operator
      template<typename IG, typename LFSU, typename X, typename Z, typename LFSV, typename Y>
      void jacobian_apply_skeleton(const IG& ig,
                                   const LFSU& lfsu_s, const X& x_s, const Z& z_s, const LFSV& lfsv_s,
                                   const LFSU& lfsu_n, const X& x_n, const Z& z_n, const LFSV& lfsv_n,
                                   Y& y_s, Y& y_n) const
      {
        W _a_i_view(_a_i, y_s.weight());
        if (ig.inside().partitionType() == Dune::InteriorEntity)
          // TODO: Only works for FastDGGridOperator (y_* being an AliasedVectorView)
          _blockOffDiagonalLOP.jacobian_apply_skeleton(ig, lfsu_s, x_s, y_s, lfsv_s, lfsu_n, x_n, y_n, lfsv_s, _a_i_view, _a_i_view);
      }

      //! Apply preconditioner after skeleton terms, linear version
      template<typename EG, typename LFSU, typename X, typename LFSV, typename Y>
      void jacobian_apply_volume_post_skeleton(const EG& eg,
                                               const LFSU& lfsu, const X& x,
                                               const LFSV& lfsv, Y& y) const
      {
        // Subtracting x: _a_i -= x_e. Looking at the algorithm above this
        // finishes step (1) since the coefficient x is the vector d
        // above. After this we have r_tmp = - a_i
        std::transform(_a_i.data(),
                       _a_i.data() + _a_i.size(),
                       x.data(),
                       _a_i.data(),
                       std::minus<value_type>{});

        // Solve the block diagonal system. This is step (2) of the
        // algorithm. After this we have _b_i = -b_i.
        _b_i = 0.0;
        W _b_i_view(_b_i, y.weight());
        _jacobianlop.jacobian_apply_volume(eg, lfsu, _a_i, lfsv, _b_i_view);


        // Do the Update step (3) of the algorithm
        std::transform(y.data(),
                       y.data() + y.size(),
                       _b_i.data(),
                       y.data(),
                       [=](const double &a,
                           const double& b)
                       {
                         return (1-_omega)*a - _omega*b;
                       });
      }

      //! apply preconditioner after skeleton terms, linearized version
      template<typename EG, typename LFSU, typename X, typename LFSV, typename Y>
      void jacobian_apply_volume_post_skeleton(const EG& eg, const LFSU& lfsu, const X& x, const X& z, const LFSV& lfsv, Y& y) const
      {
        // static checks
        //----------------------
        // static_assert(std::is_same<decltype(x),decltype(_a_i)>::value,"Both types have to be the same for nonlinear Jacobian apply");
        // static_assert(not std::is_same<decltype(x),decltype(_a_i)>::value,"Both types have to be the same for nonlinear Jacobian apply");
        //
        // dynamic crashes
        //----------------------
        // x.bar(); // type ConstAliasedVectorView
        // _a_i.bar(); // type LocalVector<double, >

        // Subtract x_e: _a_i -= x_e
        std::transform(_a_i.data(),_a_i.data()+_a_i.size(),
                       z.data(),
                       _a_i.data(),
                       std::minus<value_type>{});
        // Divide by block diagonal, _b_i = D_{ee}^{-1} _a_i
        _b_i = 0.0;
        W _b_i_view(_b_i,y.weight());
        _jacobianlop.jacobian_apply_volume(eg, lfsu, x, _a_i, lfsv, _b_i_view);
        // Update vector r_e -= omega*_b_i
        // <=> r_e = (1-\omega)r_e + \omega D_{ee}^{-1}(x_e - \sum_{e'}A_{e,e'}r_{e'})
        std::transform(y.data(),
                       y.data()+y.size(),
                       _b_i.data(),
                       y.data(),
                       [=](const double& a,
                           const double& b)
                       {
                         return (1-_omega)*a - _omega*b;
                       });
      }

    private :

      JacobianLOP& _jacobianlop;
      BlockOffDiagonalLOP& _blockOffDiagonalLOP;
      const double _omega;
      mutable LocalVector _a_i;
      mutable LocalVector _b_i;
    }; // end class BlockSORPreconditionerLocalOperator

  } // namespace PDELab
} // namespace Dune

#endif
