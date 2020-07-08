// -*- tab-width: 2; indent-tabs-mode: nil -*-
// vi: set et ts=2 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_ISTL_MATRIXFREE_ASSEMBLEDBLOCKJACOBIPREONDITIONER_HH
#define DUNE_PDELAB_BACKEND_ISTL_MATRIXFREE_ASSEMBLEDBLOCKJACOBIPREONDITIONER_HH

#if HAVE_EIGEN

#include<dune/grid/common/scsgmapper.hh>

#include <dune/pdelab/localoperator/blockdiagonalwrapper.hh>

#include<Eigen/Dense>


namespace Dune {
  namespace PDELab {

    namespace impl {
      // Below we will assemble a block of the matrix as Eigen matrix. This
      // accumulation view makes sure that we apply the weights.
      template<typename C>
      class WeightedEigenMatrixAccumulationView
      {
      public :
        typedef C Container;

        typedef typename Container::Scalar value_type;
        using size_type = typename Container::StorageIndex;
        using weight_type = value_type;

        const weight_type& weight()
        {
          return _weight;
        }

        auto data()
        {
          return _container.data();
        }

        const auto data() const
        {
          return _container.data();
        }

        auto& container()
        {
          return _container;
        }

        template <typename LFSV, typename LFSU>
        void accumulate(const LFSV& lfsv, size_type i, const LFSU& lfsu, size_type j, value_type value)
        {
          _container(i,j) += _weight * value;
        }

        WeightedEigenMatrixAccumulationView(Container& container, weight_type weight)
          : _container(container)
          , _weight(weight)
        {}

      private :
        Container& _container;
        weight_type _weight;
      };
    }


    /**\brief Local operator that can be used to create a partially matrix-free Jacobi preconditioner
     *
     * Partially matrix-free means the following: A matrix-based block Jacobi
     * preconditioner would need to invert the diagonal block matrix D that
     * comes from the block decomposition A=L+D+U in the lower left blocks,
     * diagonal blocks and upper right blocks.
     *
     * Instead of assembling the diagonal block, inverting it and applying it
     * to a vector in the Krylow solver (i.e. D^{-1}v) we never assemble the
     * matrix D or D^{-1} and do the application in a matrix free way.
     *
     * Since the inverse of a block diagonal matrix is just the block matrix
     * of the inverse of the blocks we need a way to apply the inverse of a
     * block (we denote this with D_i). This solver will assemble the diagonal
     * block and compute the inverse with LU decomposition.
     *
     * In this sense the resulting preconditioner is partially matrix-free: We
     * never assemble the complete matrix D but we will assemble all the
     * diagonal blocks D_i. For a complete matrix-free implementation have a
     * look at the class IterativeBlockJacobiPreconditionerLocalOperator.
     *
     * On the technical side this operator will calculate the LU decomposition
     * of the diagonal blocks in the alpha_volume method. Later on the inverse
     * can be applied through the jacobian_volume methods.
     *
     * For examples see dune-pdelab/dune/pdelab/test/matrixfree/.
     */
    template<class BlockDiagonalLocalOperator, class GridFunctionSpace, class DomainField>
    class AssembledBlockJacobiPreconditionerLocalOperator
      : public Dune::PDELab::FullVolumePattern
      , public Dune::PDELab::LocalOperatorDefaultFlags
    {
      // Extract some types
      using GridView = typename GridFunctionSpace::Traits::GridViewType;
      using Grid = typename GridView::Traits::Grid;
      using EntityType = typename Grid::template Codim<0>::Entity;
      using value_type = DomainField;
      using DiagonalBlock = ::Eigen::Matrix<value_type, ::Eigen::Dynamic, ::Eigen::Dynamic, ::Eigen::RowMajor>;
      using PermutationIndices = ::Eigen::Matrix<int, ::Eigen::Dynamic, 1>;
      using TupleType = std::tuple<DiagonalBlock, PermutationIndices>;

    public :

      // Since this class implements something like D^{-1}v for a diagonal
      // block matrix D we will only have volume methods. The underlying local
      // operator that describes the block diagonal might of course have
      // skeleton or boundary parts.
      static constexpr bool doPatternVolume = true;
      static constexpr bool doAlphaVolume = true;
      static constexpr bool isLinear = BlockDiagonalLocalOperator::isLinear;

      /**\brief Constructor
       *
       * \param blockDiagonalLocalOperator A local operator that can be used to
       *   make a local assembly of a diagonal block through the
       *   <assembleLocalDiagonalBlock> function. You can create such a local
       *   operator by wrapping your local operator into a
       *   <BlockDiagonalLocalOperatorWrapper>
       * \param gridFunctionSpace A grid function space
       * \param verbose Controls the amount of output
       */
      AssembledBlockJacobiPreconditionerLocalOperator(const BlockDiagonalLocalOperator& blockDiagonalLocalOperator,
                                                      const GridFunctionSpace& gridFunctionSpace,
                                                      const bool verbose=0)
        : _blockDiagonalLocalOperator(blockDiagonalLocalOperator)
        , _gridFunctionSpace(gridFunctionSpace)
        , _mapper(gridFunctionSpace.gridView().grid())
        , _precCache(_mapper.size())
        , _verbose(verbose)
        , _requireSetup(true)
      {
      }


      // Before we can call the jacobian_apply methods we need to compute the
      // LU decomposition of the diagonal block. We use the bool _requireSetup
      // to make sure that we have computed this decomposition before we try to
      // use it.
      bool requireSetup()
      {
        return _requireSetup;
      }
      void setRequireSetup(bool requireSetup)
      {
        _requireSetup = requireSetup;
      }


      /** \brief Prepare partially matrix-free preconditioner
       *
       * This assembles the local block diagonal and computes a
       * LU-decomposition with partial pivoting.
       */
      template<typename EG, typename LFSU, typename X, typename LFSV, typename Y>
      void alpha_volume(const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, Y& y) const
      {
        // Matrix for storing the LU decomposition of the block diagonal
        const int size = lfsu.size();
        assert(lfsv.size() == size);

        // Assemble local block diagonal
        DiagonalBlock diagonalBlock = DiagonalBlock::Constant(size, size, 0.0);
        impl::WeightedEigenMatrixAccumulationView<DiagonalBlock> viewDiagonalBlock(diagonalBlock, y.weight());
        assembleLocalDiagonalBlock(_blockDiagonalLocalOperator, eg, lfsu, x, lfsv, viewDiagonalBlock);

        // Compute the LU-factorization directly inside the matrix
        ::Eigen::PartialPivLU<::Eigen::Ref<DiagonalBlock>> lu(diagonalBlock);
        auto inversePermutationIndices = lu.permutationP().indices();

        // We need the inverse of the permutation matrix that Eigen gives
        PermutationIndices permutationIndices(size);
        for(int i=0; i<inversePermutationIndices.size(); ++i)
          permutationIndices[inversePermutationIndices[i]] = i;

        // Fill up correct entry of preconditioner cache
        std::size_t cache_idx = _mapper.index(eg.entity());
        _precCache[cache_idx] = std::make_tuple(diagonalBlock, permutationIndices);
      }

      template<typename EG, typename LFSU, typename X, typename Z, typename LFSV, typename Y>
      void jacobian_apply_volume(const EG& eg, const LFSU& lfsu, const X& x, const Z& z, const LFSV& lfsv, Y& y) const
      {
        // We don't need the linearization point as the preconditioner is
        // already set up
        jacobian_apply_volume(eg,lfsu,z,lfsv,y);
      }


      /**\brief Apply partially matrix-free preconditioner
       *
       * This computes D_i^{-1}z using the LU-decomposition from the
       * preparation step.
       */
      template<typename EG, typename LFSU, typename Z, typename LFSV, typename Y>
      void jacobian_apply_volume(const EG& eg, const LFSU& lfsu, const Z& z, const LFSV& lfsv, Y& y) const
      {
        assert(not _requireSetup);

        // Get correct LU decomposition of the block diagonal from the cache
        std::size_t cache_idx = _mapper.index(eg.entity());
        const auto& cacheEntry = _precCache[cache_idx];
        auto diagonalBlockLU = std::get<0>(cacheEntry);
        auto permutationIndices = std::get<1>(cacheEntry);

        // Apply the preconditioner
        const int size = lfsu.size();
        const auto& z_container = accessBaseContainer(z);
        std::vector<value_type> ztmp_container(size);
        auto& y_container = accessBaseContainer(y);

        // Forward solve with lower triangle
        for(int i=0; i<size; ++i){
          value_type rhs(z_container[permutationIndices[i]]);
          for(int j=0; j<i; ++j)
            rhs -= diagonalBlockLU(i, j) * ztmp_container[j];
          ztmp_container[i] = rhs; // L has ones on the diagonal
        }
        // Backward solve with upper triangular
        for(int i=size-1; i>=0; i--){
          value_type rhs(ztmp_container[i]);
          for(int j=size-1; j>i; j--)
            rhs -= diagonalBlockLU(i, j) * y_container[j];
          y_container[i] = rhs/diagonalBlockLU(i,i);
        }
      }

    private :
      BlockDiagonalLocalOperator _blockDiagonalLocalOperator;
      const GridFunctionSpace& _gridFunctionSpace;
      typename Dune::LeafSingleCodimSingleGeomTypeMapper<Grid, 0> _mapper;
      mutable std::vector<TupleType> _precCache;
      const int _verbose;
      bool _requireSetup;
    };

  } // namespace PDELab
} // namespace Dune

#endif // HAVE_EIGEN

#endif
