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

      /** \brief Point Jacobi preconditioner for local diagonal block systems
       *
       * The values of the inverse of the diagonal need to be provided during
       * construction.
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
      template<typename BlockDiagonalLOP, typename W, typename XView, typename EG, typename LFSU, typename LFSV>
      struct BlockDiagonalOperator
        : public Dune::LinearOperator<W,W>
      {
        using Base = Dune::LinearOperator<W,W>;
        using field_type = typename Base::field_type;
        using weight_type = typename W::WeightedAccumulationView::weight_type;
        static constexpr bool isLinear = BlockDiagonalLOP::isLinear;

        Dune::SolverCategory::Category category() const override
        {
          return Dune::SolverCategory::sequential;
        }

        BlockDiagonalOperator(const BlockDiagonalLOP& blockDiagonalLop,
                              const EG& eg,
                              const LFSU& lfsu,
                              const LFSV& lfsv)
          : _blockDiagonalLop(blockDiagonalLop)
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
            applyLocalDiagonalBlock(_blockDiagonalLop, _eg, _lfsu, z, z, _lfsv, y_view);
          else
            applyLocalDiagonalBlock(_blockDiagonalLop, _eg, _lfsu, *_u, z, _lfsv, y_view);
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
        const BlockDiagonalLOP& _blockDiagonalLop;
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

        // Note: This whole method is only needed since the generated point
        // diagonal localoperator use matrix based accumulation calls for
        // now. This is something that needs to be fixed in the nonfast case of
        // code generation
        template <typename LFSV, typename LFSU>
        void accumulate(const LFSV& lfsv, size_type i, const LFSU& lfsu, size_type j, value_type value)
        {
          assert (i==j);

          // palpo TODO: The if statement is only a workaround. The nonfast
          // generated local operator need to be adjusted to drop the
          // nondiagonal entries...
          if (i == j)
            _container.data()[lfsv.localIndex(i)] += _weight * value;
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


    /** \struct BlockSolverOptions
     * \brief Options for IterativeBlockJacobiPreconditionerLocalOperator
     *
     * Controls the options of the iterative solver for the individual
     * blocks.
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


    template<typename BlockDiagonalLOP,
             typename PointDiagonalLOP,
             typename GridFunctionSpace,
             typename DomainField,
             template<typename> class IterativeSolver>
    class IterativeBlockJacobiPreconditionerLocalOperator
      : public Dune::PDELab::FullVolumePattern
      , public Dune::PDELab::LocalOperatorDefaultFlags
    {
      using GridView = typename GridFunctionSpace::Traits::GridViewType;
      using Grid = typename GridView::Traits::Grid;
      using EntityType = typename Grid::template Codim<0>::Entity;
      using value_type = DomainField;

      using LocalVector = Dune::PDELab::LocalVector<value_type>;
      using InvDiagonal = typename LocalVector::BaseContainer;

    public:
      static constexpr bool doPatternVolume = true;
      static constexpr bool doAlphaVolume = true;
      static constexpr bool isLinear = BlockDiagonalLOP::isLinear;

      IterativeBlockJacobiPreconditionerLocalOperator(const BlockDiagonalLOP& blockDiagonalLop,
                                                      const PointDiagonalLOP& pointDiagonalLop,
                                                      const GridFunctionSpace& gridFunctionSpace,
                                                      SolverStatistics<int>& solver_stat,
                                                      BlockSolverOptions solveroptions,
                                                      const bool verbose=0)
        : _blockDiagonalLOP(blockDiagonalLop)
        , _pointDiagonalLop(pointDiagonalLop)
        , _gridFunctionSpace(gridFunctionSpace)
        , _solverStat(solver_stat)
        , _mapper(gridFunctionSpace.gridView().grid())
        , _invDiagonalCache(_mapper.size())
        , _solveroptions(solveroptions)
        , _verbose(verbose)
        , _requireSetup(true)
      {}

      bool requireSetup()
      {
        return _requireSetup;
      }

      void setRequireSetup(bool v)
      {
        _requireSetup = v;
      }

      /** \brief prepare tensor product preconditioner

          - call the point diagonal operator
          - store inverses of the result
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
        assembleLocalPointDiagonal(_pointDiagonalLop, eg, lfsu, x, lfsv, view);

        // Invert this once now to be able to do multiplications lateron
        for(auto& val : _invDiagonalCache[cache_idx]){
          val = 1. / val;
        }
      }


      // Jacobian apply for linear problems
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
        impl::BlockDiagonalOperator<BlockDiagonalLOP, LocalVector, Z, EG, LFSU, LFSV> op(_blockDiagonalLOP,
                                                                                         eg,
                                                                                         lfsu,
                                                                                         lfsv);
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
        _solverStat.append(stat.iterations);
        std::transform(y.data(),
                       y.data()+y.size(),
                       y_tmp.data(),
                       y.data(),
                       std::plus<value_type>{});
      } // end jacobian_apply_volume


      // Jacobian apply for nonlinear problems
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
        impl::BlockDiagonalOperator<BlockDiagonalLOP, LocalVector, X, EG, LFSU, LFSV> op(_blockDiagonalLOP,
                                                                                         eg,
                                                                                         lfsu,
                                                                                         lfsv);
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
        _solverStat.append(stat.iterations);
        std::transform(y.data(),
                       y.data()+y.size(),
                       y_tmp.data(),
                       y.data(),
                       std::plus<value_type>{});
      } // end jacobian_apply_volume

    private :
      BlockDiagonalLOP _blockDiagonalLOP;
      PointDiagonalLOP _pointDiagonalLop;
      const GridFunctionSpace& _gridFunctionSpace;
      SolverStatistics<int>& _solverStat;
      typename Dune::LeafSingleCodimSingleGeomTypeMapper<Grid, 0> _mapper;
      mutable std::vector<InvDiagonal> _invDiagonalCache;
      mutable BlockSolverOptions _solveroptions;
      const int _verbose;
      bool _requireSetup;
    }; // end class IterativeBlockJacobiPreconditionerLocalOperator

  } // namespace PDELab
} // namespace Dune
#endif
