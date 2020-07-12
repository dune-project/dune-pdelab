// -*- tab-width: 2; indent-tabs-mode: nil -*-
// vi: set et ts=2 sw=2 sts=2:

#ifndef DUNE_PDELAB_BACKEND_ISTL_MATRIXFREE_BACKENDS_HH
#define DUNE_PDELAB_BACKEND_ISTL_MATRIXFREE_BACKENDS_HH

#include<dune/grid/common/rangegenerators.hh>
#include<dune/istl/preconditioner.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/solvers.hh>

#include<dune/pdelab/backend/istl/matrixfree/gridoperatorpreconditioner.hh>
#include<dune/pdelab/gridoperator/fastdg.hh>
#include<dune/pdelab/backend/istl/matrixfree/blocksorpreconditioner.hh>

namespace Dune {
  namespace PDELab {

    namespace impl{

      // Can be used to check if a local operator is a
      // BlockSORPreconditionerLocalOperator
      template <typename>
      struct isBlockSORPreconditionerLocalOperator : std::false_type {};

      template <typename JacobianLOP, typename BlockOffDiagonalLOP, typename GridFunctionSpace>
      struct isBlockSORPreconditionerLocalOperator<BlockSORPreconditionerLocalOperator<JacobianLOP,
                                                                                       BlockOffDiagonalLOP,
                                                                                       GridFunctionSpace>>
        : std::true_type {};

      // Can be used to check if a grid operator is a FastDGGridOperator
      template <typename>
      struct isFastDGGridOperator : std::false_type {};

      template<typename GFSU, typename GFSV, typename LOP,
               typename MB, typename DF, typename RF, typename JF,
               typename CU, typename CV>
      struct isFastDGGridOperator<FastDGGridOperator<GFSU, GFSV, LOP, MB, DF, RF, JF, CU, CV >>
        : std::true_type {};

    }

    /**\brief Sequential matrix-free solver backend
     *
     * This can be used with a combination of {CGSolver, BiCGSTABSolver,
     * MINRESSolver} as a solver and grid operator build from one of the
     * following local operators:
     *
     * - AssembledBlockJacobiPreconditionerLocalOperator for partial matrix-free Jacobi
     * - IterativeBlockJacobiPreconditionerLocalOperator for fully matrix-free Jacobi
     * - BlockSORPreconditionerLocalOperator for matrix-free SOR
     *
     * Note: If you use BlockSORPreconditionerLocalOperator you need to use
     * FastDGGridOperator!

     *  \tparam GO     Grid operator implementing the matrix-free operator application
     *  \tparam PrecGO Grid operator implementing matrix-free preconditioning
     */
    template<class GO, class PrecGO,
             template<class> class Solver>
    class ISTLBackend_SEQ_MatrixFree_Base
      : public Dune::PDELab::SequentialNorm
      , public Dune::PDELab::LinearResultStorage
    {
      using V = typename GO::Traits::Domain;
      using W = typename GO::Traits::Range;
      using Operator = Dune::PDELab::OnTheFlyOperator<V,W,GO>;
      using SeqPrec = GridOperatorPreconditioner<PrecGO>;

    public :
      ISTLBackend_SEQ_MatrixFree_Base(const GO& go, const PrecGO& precgo,
                                      unsigned maxiter=5000, int verbose=1)
        : opa_(go), seqprec_(precgo)
        , maxiter_(maxiter), verbose_(verbose)
      {
        // First lets brake that down: The case we want to avoid is having SOR
        // with regular grid operator. We check for SOR and not fast grid
        // operator and assert that this is not the case.
        //
        // For more information why this is not working have a look at the
        // documentation of the class BlockSORPreconditionerLocalOperator
        static_assert(not(impl::isBlockSORPreconditionerLocalOperator<
                          typename PrecGO::Traits::LocalAssembler::LocalOperator>::value
                          and not impl::isFastDGGridOperator<PrecGO>::value),
                      "If you use the BlockSORPreconditioner you need to use FastDGGridOperator!");
      }

      void apply(V& z, W& r, typename V::ElementType reduction)
      {
        Solver<V> solver(opa_,seqprec_,reduction,maxiter_,verbose_);
        Dune::InverseOperatorResult stat;
        solver.apply(z, r, stat);
        res.converged  = stat.converged;
        res.iterations = stat.iterations;
        res.elapsed    = stat.elapsed;
        res.reduction  = stat.reduction;
        res.conv_rate  = stat.conv_rate;
      }

      //! Set position of jacobian.
      //! Must be called before apply().
      void setLinearizationPoint(const V& u)
      {
        opa_.setLinearizationPoint(u);
        seqprec_.setLinearizationPoint(u);
      }

    private :
      Operator opa_;
      SeqPrec seqprec_;
      unsigned maxiter_;
      int verbose_;
    }; // end class ISTLBackend_SEQ_MatrixFree_Base


    // A matrix based SOR backend for comparing the matrix-free version
    class ISTLBackend_SEQ_BCGS_SOR
      : public ISTLBackend_SEQ_Base<Dune::SeqSOR, Dune::BiCGSTABSolver>
    {
    public:
      /*! \brief make a linear solver object

        \param[in] maxiter_ maximum number of iterations to do
        \param[in] verbose_ print messages if true
      */
      explicit ISTLBackend_SEQ_BCGS_SOR (unsigned maxiter_=5000, int verbose_=1)
        : ISTLBackend_SEQ_Base<Dune::SeqSOR, Dune::BiCGSTABSolver>(maxiter_, verbose_)
      {}
    };

  } // namespace PDELab
} // namespace Dune
#endif
