// -*- tab-width: 2; indent-tabs-mode: nil -*-
// vi: set et ts=2 sw=2 sts=2:

#ifndef DUNE_PDELAB_BACKEND_ISTL_MATRIXFREE_BACKENDS_HH
#define DUNE_PDELAB_BACKEND_ISTL_MATRIXFREE_BACKENDS_HH

#include<dune/grid/common/rangegenerators.hh>
#include<dune/istl/preconditioner.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/solvers.hh>

#include<dune/pdelab/backend/istl/matrixfree/gridoperatorpreconditioner.hh>

namespace Dune {
  namespace PDELab {

    /** \brief Sequential matrix-free solver backend for the combinations {CGSolver,BiCGSTABSolver,MINRESSolver} x {SeqMatrixFreeBlockSOR}.

        \tparam GO     Grid operator implementing the matrix-free operator application, plugged in SeqMatrixFreeBlockSOR.
        \tparam PrecGO Grid operator implementing matrix-free preconditioning. Plugged into SeqMatrixFreeBlockSOR.

        The grid operator \p PrecGO can take as a \p LOP
        the local operator wrapper IterativeInverseBlockDiagonalLocalOperatorWrapper
        to build a matrix-free undamped block Jacobi only usable with zero input
        or the local operator wrapper BlockDiagonalSORLocalOperatorWrapper to build a matrix-free block SOR preconditioner.
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
      {}

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
