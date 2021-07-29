// -*- tab-width: 2; indent-tabs-mode: nil -*-
// vi: set et ts=2 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_ISTL_MATRIXFREE_ITERATIVEBLOCKJACOBIPRECONDITIONER_HH
#define DUNE_PDELAB_BACKEND_ISTL_MATRIXFREE_ITERATIVEBLOCKJACOBIPRECONDITIONER_HH

#include<dune/grid/common/scsgmapper.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/operators.hh>

#include <dune/pdelab/localoperator/pointdiagonalwrapper.hh>
#include <dune/pdelab/localoperator/blockdiagonalwrapper.hh>

#include<dune/pdelab/backend/istl/matrixfree/solverstatistics.hh>
#include<dune/pdelab/backend/istl/matrixfree/assembledblockjacobipreconditioner.hh>

namespace Dune {
  namespace PDELab {

    namespace impl{

      /* \brief Point Jacobi preconditioner for local diagonal block systems
       *
       * This class will do the Jacobi preconditioning for a diagonal
       * block. The values of the inverse of the diagonal need to be provided
       * during construction.
       */
      template<typename X>
      class LocalPointJacobiPreconditioner : public Dune::Preconditioner<X,X>
      {
      public :
        typedef X domain_type;
        typedef X range_type;
        typedef typename X::BaseContainer InvDiagonal;
        typedef typename X::value_type value_type;

        Dune::SolverCategory::Category category() const override
        {
          return Dune::SolverCategory::sequential;
        }

        /** \brief Constructor
         *
         * \param[in] invDiagonal Vector holding the inverse of the matrix
         * \param[in] diagonalWeight Single value for scaling
         * \param[in] precondition Bool that decides if the preconditioner should be applied or not
         */
        LocalPointJacobiPreconditioner(const InvDiagonal& invDiagonal,
                                       const value_type diagonalWeight,
                                       const bool precondition=true)
          : _precondition(precondition)
          , _invDiagonal(invDiagonal)
          , _diagonalWeight(diagonalWeight)
        {}

        void pre (domain_type& x, range_type& b) override {}

        void apply (domain_type& v, const range_type& d) override
        {
          if(_precondition){
            std::transform(d.data(),
                           d.data()+d.size(),
                           _invDiagonal.begin(),
                           v.data(),
                           [=](const value_type& a, const value_type& b){
                             return _diagonalWeight*a*b;
                           });
          }
          else{
            std::copy(d.data(),d.data()+d.size(),v.data());
          }
        }

        void post (domain_type& x) override {}

      private :
        const bool _precondition;
        const InvDiagonal& _invDiagonal;
        const value_type _diagonalWeight;
      };


      /** \brief Create ISTL operator that applies a local block diagonal
       *
       * This class is similar to the OnTheFlyOperator in the sense that a
       * local operator is used to create an ISTL operator can be used to solve
       * a linear system of equations in a matrix free way.
       *
       * Where the OnTheFlyOperator takes a grid operator and solves a global
       * system, this class only solves the local block diagonal system.
       */
      template<typename BlockDiagonalLocalOperator, typename W, typename XView, typename EG, typename LFSU, typename LFSV>
      struct BlockDiagonalOperator
        : public Dune::LinearOperator<W,W>
      {
        using Base = Dune::LinearOperator<W,W>;
        using field_type = typename Base::field_type;
        using weight_type = typename W::WeightedAccumulationView::weight_type;
        static constexpr bool isLinear = BlockDiagonalLocalOperator::isLinear;

        Dune::SolverCategory::Category category() const override
        {
          return Dune::SolverCategory::sequential;
        }

        BlockDiagonalOperator(const BlockDiagonalLocalOperator& blockDiagonalLocalOperator,
                              const EG& eg,
                              const LFSU& lfsu,
                              const LFSV& lfsv)
          : _blockDiagonalLocalOperator(blockDiagonalLocalOperator)
          , _eg(eg)
          , _lfsu(lfsu)
          , _lfsv(lfsv)
          , _tmp(lfsu.size(),0.0)
          , _u(nullptr)
        {}

        void apply(const W& z, W& y) const override
        {
          y = 0.0;
          typename W::WeightedAccumulationView y_view(y, _weight);
          if (isLinear)
            applyLocalDiagonalBlock(_blockDiagonalLocalOperator, _eg, _lfsu, z, z, _lfsv, y_view);
          else
            applyLocalDiagonalBlock(_blockDiagonalLocalOperator, _eg, _lfsu, *_u, z, _lfsv, y_view);
        }

        void applyscaleadd(field_type alpha, const W& z, W& y) const override
        {
          apply(z, _tmp);
          y.axpy(alpha, _tmp);
        }

        void setLinearizationPoint(const XView& u)
        {
          _u = &u;
        }

        void setWeight(const weight_type& weight)
        {
          _weight = weight;
        }

      private :
        const BlockDiagonalLocalOperator& _blockDiagonalLocalOperator;
        const EG& _eg;
        const LFSU& _lfsu;
        const LFSV& _lfsv;
        mutable W _tmp;
        const XView* _u;
        weight_type _weight;
      };

      template<typename C>
      class WeightedPointDiagonalAccumulationView
      {
      public :
        using Container = C;

        using value_type = typename Container::value_type;
        using size_type = typename Container::size_type;
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

        template<typename LFSV>
        void accumulate(const LFSV& lfsv, size_type i, value_type value)
        {
          _container.data()[lfsv.localIndex(i)] += _weight * value;
        }

        WeightedPointDiagonalAccumulationView(Container& container, weight_type weight)
          : _container(container)
          , _weight(weight)
        {}

      private :
        Container& _container;
        weight_type _weight;
      };

    } // namespace impl


    /** \brief Options for IterativeBlockJacobiPreconditionerLocalOperator
     *
     * Controls the options of the iterative solver for the individual blocks.
     */
    struct BlockSolverOptions
    {
      /** \brief Constructor
       *
       * \param[in] resreduction Residual reduction
       * \param[in] maxiter Maximal number of iterations
       * \param[in] precondition Precondition each block with point-Jacobi?
       * \param[in] weight Weight for point-Jacobi preconditioner
       * \param[in] verbose Verbosity level
       */
      BlockSolverOptions(const double resreduction=1e-5,
                         const size_t maxiter=100,
                         const bool precondition=true,
                         const double weight=1.0,
                         const int verbose=0)
        : _resreduction(resreduction)
        , _maxiter(maxiter)
        , _precondition(precondition)
        , _weight(weight)
        , _verbose(verbose)
      {}
      /** \brief Residual reduction, i.e. solver accuracy */
      double _resreduction;
      /** \brief Maximal number of iterations */
      size_t _maxiter;
      /** \brief Precondition with point-Jacobi? */
      bool _precondition;
      /** \brief Weight for point-jacobi */
      double _weight;
      /** \brief verbosity level */
      int _verbose;
    }; // end struct BlockSolverOptions


    /**\brief Local operator that can be used to create a fully matrix-free Jacobi preconditioner
     *
     * Similar to the partial matrix-free class
     * AssembledBlockJacobiPreconditionerLocalOperator this implements a local
     * operator that can be used to implement a matrix-free Jacobi
     * preconditioner. In contrast to the other class this implementation will
     * be fully matrix-free.
     *
     * A matrix-free Jacobi preconditioner needs to be able to apply the
     * inverse of the block diagonal D to a vector. In the partial matrix-free
     * class mentioned above this was done by assembling the diagonal blocks
     * and inverting this matrix. In order to get a fully matrix-free version
     * we instead use a Krylow solver on the diagonal block. This can once
     * again be done in a matrix-free way but we will need a preconditioner for
     * iterative solver. For this purpose we use another Jacobi preconditioner
     * that operates on a single block. This means we need to assemble the
     * point diagonal of this block and apply the inverse.
     *
     * For examples see dune-pdelab/dune/pdelab/test/matrixfree/.
     *
     * \tparam BlockDiagonalLocalOperator A lop for local application of diagonal blocks
     * \tparam PointDiagonalLocalOperator A lop for local assembly of point diagonal
     * \tparam GridFunctionSpace A grid function space
     * \tparam DomainField Domain field
     * \tparam IterativeSolver Solver that will be used to 'invert' the diagonal blocks
     */
    template<typename BlockDiagonalLocalOperator,
             typename PointDiagonalLocalOperator,
             typename GridFunctionSpace,
             typename DomainField,
             template<typename> class IterativeSolver>
    class IterativeBlockJacobiPreconditionerLocalOperator
      : public Dune::PDELab::FullVolumePattern
      , public Dune::PDELab::LocalOperatorDefaultFlags
    {
      // Extract some types
      using GridView = typename GridFunctionSpace::Traits::GridViewType;
      using Grid = typename GridView::Traits::Grid;
      using EntityType = typename Grid::template Codim<0>::Entity;
      using value_type = DomainField;
      using LocalVector = Dune::PDELab::LocalVector<value_type>;
      using InvDiagonal = typename LocalVector::BaseContainer;

    public:
      // Since this class implements something like D^{-1}v for a diagonal
      // block matrix D we will only have volume methods. The underlying local
      // operator that describes the block diagonal might of course have
      // skeleton or boundary parts.
      static constexpr bool doPatternVolume = true;
      static constexpr bool doAlphaVolume = true;
      static constexpr bool isLinear = BlockDiagonalLocalOperator::isLinear;

      /**\brief Constructor
       *
       * \param blockDiagonalLocalOperator A local operator that can be used to (locally)
       *   apply a diagonal block through the <applyLocalDiagonalBlock> function. You can
       *   create such a local operator by wrapping your local operator into a
       *   <BlockDiagonalLocalOperatorWrapper>
       * \param pointDiagonalLocalOperator A local operator that can be used to (locally) assemble
       *   the point diagonal of a diagonal block. You can create such a local operator by
       *   wrapping your local operator into a <PointDiagonalLocalOperatorWrapper>
       * \param gridFunctionSpace A grid function space
       * \param solverStat A class for export solver statistics
       * \param solveroptions A class for configuring your solver
       * \param verbose Controls the amount of output
       */
      IterativeBlockJacobiPreconditionerLocalOperator(const BlockDiagonalLocalOperator& blockDiagonalLocalOperator,
                                                      const PointDiagonalLocalOperator& pointDiagonalLocalOperator,
                                                      const GridFunctionSpace& gridFunctionSpace,
                                                      SolverStatistics<int>& solverStatiscits,
                                                      BlockSolverOptions solveroptions,
                                                      const bool verbose=0)
        : _blockDiagonalLocalOperator(blockDiagonalLocalOperator)
        , _pointDiagonalLocalOperator(pointDiagonalLocalOperator)
        , _gridFunctionSpace(gridFunctionSpace)
        , _solverStatistics(solverStatiscits)
        , _mapper(gridFunctionSpace.gridView())
        , _invDiagonalCache(_mapper.size())
        , _solveroptions(solveroptions)
        , _verbose(verbose)
        , _requireSetup(true)
      {}

      // Before we can call the jacobian_apply methods we need to assemble the
      // point diagonal of the diagonal block. This will be used as a preconditioner
      // for the iterative matrix free solver on the diagonal block.
      bool requireSetup()
      {
        return _requireSetup;
      }
      void setRequireSetup(bool v)
      {
        _requireSetup = v;
      }

      /** \brief Prepare fully matrix-free preconditioner
       *
       * This assembles the point diagonal of the diagonal block and stores its
       * inverse. During the apply step this will be used as a preconditioner
       * for the solver on the local block.
       */
      template<typename EG, typename LFSU, typename X, typename LFSV, typename Y>
      void alpha_volume(const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, Y& y) const
      {
        const int size = lfsu.size();
        assert(lfsv.size() == size);

        // Assemble point diagonal
        std::size_t cache_idx = _mapper.index(eg.entity());
        _invDiagonalCache[cache_idx].resize(lfsu.size());
        using TmpView = impl::WeightedPointDiagonalAccumulationView<InvDiagonal>;
        TmpView view(_invDiagonalCache[cache_idx], y.weight());
        assembleLocalPointDiagonal(_pointDiagonalLocalOperator, eg, lfsu, x, lfsv, view);

        // Invert this and store it (will later be used as a preconditioner)
        for(auto& val : _invDiagonalCache[cache_idx]){
          val = 1. / val;
        }
      }


      //! Apply fully matrix-free preconditioner - linear case
      template<typename EG, typename LFSU, typename Z, typename LFSV, typename Y>
      void jacobian_apply_volume(const EG& eg, const LFSU& lfsu, const Z& z, const LFSV& lfsv, Y& y) const
      {
        assert(not _requireSetup);

        // Get correct point diagonal from the cache
        std::size_t cache_idx = _mapper.index(eg.entity());
        const auto& invDiagonal = _invDiagonalCache[cache_idx];

        // Create iterative solver with point Jacobi preconditioner for solving
        // the local block diagonal system
        impl::LocalPointJacobiPreconditioner<LocalVector> pointJacobi(invDiagonal,
                                                                      _solveroptions._weight,
                                                                      _solveroptions._precondition);
        impl::BlockDiagonalOperator<BlockDiagonalLocalOperator, LocalVector, Z, EG, LFSU, LFSV>
          op(_blockDiagonalLocalOperator, eg, lfsu, lfsv);
        op.setWeight(y.weight());
        IterativeSolver<LocalVector> solver(op,
                                            pointJacobi,
                                            _solveroptions._resreduction,
                                            _solveroptions._maxiter,
                                            _solveroptions._verbose);
        Dune::InverseOperatorResult stat;

        // Set up local right hand side
        LocalVector z_tmp(lfsu.size(), 0.0);
        std::copy(z.data(), z.data()+z.size(), z_tmp.data());

        // Solve local blocks iteratively
        LocalVector y_tmp(lfsv.size(), 0.0);
        solver.apply(y_tmp, z_tmp, stat);
        _solverStatistics.append(stat.iterations);
        std::transform(y.data(),
                       y.data()+y.size(),
                       y_tmp.data(),
                       y.data(),
                       std::plus<value_type>{});
      } // end jacobian_apply_volume


      /**\brief Apply fully matrix-free preconditioner - nonlinear case
       *
       * Compared to the linear case this needs to set the correct
       * linearization point in the BlockDiagonalOperator
       */
      template<typename EG, typename LFSU, typename X, typename Z, typename LFSV, typename Y>
      void jacobian_apply_volume(const EG& eg, const LFSU& lfsu, const X& x, const Z& z, const LFSV& lfsv, Y& y) const
      {
        assert(not _requireSetup);

        // Get correct point diagonal from the cache
        std::size_t cache_idx = _mapper.index(eg.entity());
        const auto& invDiagonal = _invDiagonalCache[cache_idx];

        // Create iterative solver with point Jacobi preconditioner for solving
        // the local block diagonal system
        impl::LocalPointJacobiPreconditioner<LocalVector> pointJacobi(invDiagonal,
                                                                      _solveroptions._weight,
                                                                      _solveroptions._precondition);
        impl::BlockDiagonalOperator<BlockDiagonalLocalOperator, LocalVector, X, EG, LFSU, LFSV>
          op(_blockDiagonalLocalOperator, eg, lfsu, lfsv);
        op.setLinearizationPoint(x);
        op.setWeight(y.weight());
        IterativeSolver<LocalVector> solver(op,
                                            pointJacobi,
                                            _solveroptions._resreduction,
                                            _solveroptions._maxiter,
                                            _solveroptions._verbose);
        Dune::InverseOperatorResult stat;

        // Set up local right hand side
        LocalVector z_tmp(lfsu.size(), 0.0);
        std::copy(z.data(), z.data()+z.size(), z_tmp.data());

        // Solve local blocks iteratively
        LocalVector y_tmp(lfsv.size(), 0.0);
        solver.apply(y_tmp, z_tmp, stat);
        _solverStatistics.append(stat.iterations);
        std::transform(y.data(),
                       y.data()+y.size(),
                       y_tmp.data(),
                       y.data(),
                       std::plus<value_type>{});
      } // end jacobian_apply_volume

    private :
      BlockDiagonalLocalOperator _blockDiagonalLocalOperator;
      PointDiagonalLocalOperator _pointDiagonalLocalOperator;
      const GridFunctionSpace& _gridFunctionSpace;
      SolverStatistics<int>& _solverStatistics;
      typename Dune::SingleCodimSingleGeomTypeMapper<GridView, 0> _mapper;
      mutable std::vector<InvDiagonal> _invDiagonalCache;
      mutable BlockSolverOptions _solveroptions;
      const int _verbose;
      bool _requireSetup;
    }; // end class IterativeBlockJacobiPreconditionerLocalOperator

  } // namespace PDELab
} // namespace Dune
#endif
