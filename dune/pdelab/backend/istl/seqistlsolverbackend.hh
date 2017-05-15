// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_ISTL_SEQISTLSOLVERBACKEND_HH
#define DUNE_PDELAB_BACKEND_ISTL_SEQISTLSOLVERBACKEND_HH

// this is here for backwards compatibility and deprecation warnings, remove after 2.5.0
#include "ensureistlinclude.hh"

#include <dune/common/deprecated.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/solvercategory.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/scalarproducts.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/istl/paamg/pinfo.hh>
#include <dune/istl/io.hh>
#include <dune/istl/superlu.hh>
#include <dune/istl/umfpack.hh>

#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include <dune/pdelab/backend/solver.hh>
#include <dune/pdelab/backend/istl/vector.hh>
#include <dune/pdelab/backend/istl/bcrsmatrix.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup Backend
    //! \ingroup PDELab
    //! \{

    template<typename X, typename Y, typename GOS>
    class OnTheFlyOperator : public Dune::LinearOperator<X,Y>
    {
    public:
      typedef X domain_type;
      typedef Y range_type;
      typedef typename X::field_type field_type;

      OnTheFlyOperator (GOS& gos_)
        : gos(gos_)
      {}

      virtual void apply (const X& x, Y& y) const
      {
        y = 0.0;
        gos.jacobian_apply(x,y);
      }

      virtual void applyscaleadd (field_type alpha, const X& x, Y& y) const
      {
        Y temp(y);
        temp = 0.0;
        gos.jacobian_apply(x,temp);
        y.axpy(alpha,temp);
      }

      SolverCategory::Category category() const override
      {
        return SolverCategory::sequential;
      }

    private:
      GOS& gos;
    };

    //==============================================================================
    // Here we add some standard linear solvers conforming to the linear solver
    // interface required to solve linear and nonlinear problems.
    //==============================================================================

    template<template<class,class,class,int> class Preconditioner,
             template<class> class Solver>
    class ISTLBackend_SEQ_Base
      : public SequentialNorm, public LinearResultStorage
    {
    public:
      /*! \brief make a linear solver object

        \param[in] maxiter_ maximum number of iterations to do
        \param[in] verbose_ print messages if true
      */
      explicit ISTLBackend_SEQ_Base(unsigned maxiter_=5000, int verbose_=1)
        : maxiter(maxiter_), verbose(verbose_)
      {}



      /*! \brief solve the given linear system

        \param[in] A the given matrix
        \param[out] z the solution vector to be computed
        \param[in] r right hand side
        \param[in] reduction to be achieved
      */
      template<class M, class V, class W>
      void apply(M& A, V& z, W& r, typename Dune::template FieldTraits<typename W::ElementType >::real_type reduction)
      {
        using Backend::Native;
        using Backend::native;

        Dune::MatrixAdapter<Native<M>,
                            Native<V>,
                            Native<W>> opa(native(A));
        Preconditioner<Native<M>,
                       Native<V>,
                       Native<W>,
                       1> prec(native(A), 3, 1.0);
        Solver<Native<V>> solver(opa, prec, reduction, maxiter, verbose);
        Dune::InverseOperatorResult stat;
        solver.apply(native(z), native(r), stat);
        res.converged  = stat.converged;
        res.iterations = stat.iterations;
        res.elapsed    = stat.elapsed;
        res.reduction  = stat.reduction;
        res.conv_rate  = stat.conv_rate;
      }

    private:
      unsigned maxiter;
      int verbose;
    };

    template<template<typename> class Solver>
    class ISTLBackend_SEQ_ILU0
      :  public SequentialNorm, public LinearResultStorage
    {
    public:
      /*! \brief make a linear solver object

        \param[in] maxiter_ maximum number of iterations to do
        \param[in] verbose_ print messages if true
      */
      explicit ISTLBackend_SEQ_ILU0 (unsigned maxiter_=5000, int verbose_=1)
        : maxiter(maxiter_), verbose(verbose_)
       {}
      /*! \brief solve the given linear system

        \param[in] A the given matrix
        \param[out] z the solution vector to be computed
        \param[in] r right hand side
        \param[in] reduction to be achieved
      */
      template<class M, class V, class W>
      void apply(M& A, V& z, W& r, typename Dune::template FieldTraits<typename W::ElementType >::real_type reduction)
      {
        using Backend::Native;
        using Backend::native;
        Dune::MatrixAdapter<Native<M>,
                            Native<V>,
                            Native<W>> opa(native(A));
        Dune::SeqILU0<Native<M>,
                      Native<V>,
                      Native<W>
                      > ilu0(native(A), 1.0);
        Solver<Native<V>> solver(opa, ilu0, reduction, maxiter, verbose);
        Dune::InverseOperatorResult stat;
        solver.apply(native(z), native(r), stat);
        res.converged  = stat.converged;
        res.iterations = stat.iterations;
        res.elapsed    = stat.elapsed;
        res.reduction  = stat.reduction;
        res.conv_rate  = stat.conv_rate;
       }
    private:
      unsigned maxiter;
      int verbose;
    };

    template<template<typename> class Solver>
    class ISTLBackend_SEQ_ILUn
      :  public SequentialNorm, public LinearResultStorage
    {
    public:
      /*! \brief make a linear solver object
        \param[in] n The number of levels to be used.
        \param[in] w The relaxation factor.
        \param[in] maxiter_ maximum number of iterations to do
        \param[in] verbose_ print messages if true
      */
      ISTLBackend_SEQ_ILUn (int n, double w, unsigned maxiter_=5000, int verbose_=1)
        : n_(n), w_(w), maxiter(maxiter_), verbose(verbose_)
       {}
      /*! \brief solve the given linear system

        \param[in] A the given matrix
        \param[out] z the solution vector to be computed
        \param[in] r right hand side
        \param[in] reduction to be achieved
      */
      template<class M, class V, class W>
      void apply(M& A, V& z, W& r, typename Dune::template FieldTraits<typename W::ElementType >::real_type reduction)
      {
        using Backend::Native;
        using Backend::native;
        Dune::MatrixAdapter<Native<M>,
                            Native<V>,
                            Native<W>
                            > opa(native(A));
        Dune::SeqILUn<Native<M>,
                      Native<V>,
                      Native<W>
                      > ilun(native(A), n_, w_);
        Solver<Native<V>> solver(opa, ilun, reduction, maxiter, verbose);
        Dune::InverseOperatorResult stat;
        solver.apply(native(z), native(r), stat);
        res.converged  = stat.converged;
        res.iterations = stat.iterations;
        res.elapsed    = stat.elapsed;
        res.reduction  = stat.reduction;
        res.conv_rate  = stat.conv_rate;
       }
    private:
      int n_;
      double w_;

      unsigned maxiter;
      int verbose;
    };

    //! \addtogroup PDELab_seqsolvers Sequential Solvers
    //! \{

    /**
     * @brief Backend for sequential loop solver with Jacobi preconditioner.
     */
    class ISTLBackend_SEQ_LOOP_Jac
      : public ISTLBackend_SEQ_Base<Dune::SeqJac, Dune::LoopSolver>
    {
    public:
      /*! \brief make a linear solver object
        \param[in] maxiter_ maximum number of iterations to do
        \param[in] verbose_ print messages if true
      */
      explicit ISTLBackend_SEQ_LOOP_Jac (unsigned maxiter_=5000, int verbose_=1)
        : ISTLBackend_SEQ_Base<Dune::SeqJac, Dune::LoopSolver>(maxiter_, verbose_)
      {}
    };

    /**
     * @brief Backend for sequential BiCGSTAB solver with Jacobi preconditioner.
     */
    class ISTLBackend_SEQ_BCGS_Jac
      : public ISTLBackend_SEQ_Base<Dune::SeqJac, Dune::BiCGSTABSolver>
    {
    public:
      /*! \brief make a linear solver object
        \param[in] maxiter_ maximum number of iterations to do
        \param[in] verbose_ print messages if true
      */
      explicit ISTLBackend_SEQ_BCGS_Jac (unsigned maxiter_=5000, int verbose_=1)
        : ISTLBackend_SEQ_Base<Dune::SeqJac, Dune::BiCGSTABSolver>(maxiter_, verbose_)
      {}
    };

    /**
     * @brief Backend for sequential BiCGSTAB solver with SSOR preconditioner.
     */
    class ISTLBackend_SEQ_BCGS_SSOR
      : public ISTLBackend_SEQ_Base<Dune::SeqSSOR, Dune::BiCGSTABSolver>
    {
    public:
      /*! \brief make a linear solver object

        \param[in] maxiter_ maximum number of iterations to do
        \param[in] verbose_ print messages if true
      */
      explicit ISTLBackend_SEQ_BCGS_SSOR (unsigned maxiter_=5000, int verbose_=1)
        : ISTLBackend_SEQ_Base<Dune::SeqSSOR, Dune::BiCGSTABSolver>(maxiter_, verbose_)
      {}
    };

     /**
     * @brief Backend for sequential BiCGSTAB solver with ILU0 preconditioner.
     */
    class ISTLBackend_SEQ_BCGS_ILU0
      : public ISTLBackend_SEQ_ILU0<Dune::BiCGSTABSolver>
    {
    public:
      /*! \brief make a linear solver object

        \param[in] maxiter_ maximum number of iterations to do
        \param[in] verbose_ print messages if true
      */
      explicit ISTLBackend_SEQ_BCGS_ILU0 (unsigned maxiter_=5000, int verbose_=1)
        : ISTLBackend_SEQ_ILU0<Dune::BiCGSTABSolver>(maxiter_, verbose_)
      {}
    };

    /**
     * @brief Backend for sequential conjugate gradient solver with ILU0 preconditioner.
     */
    class ISTLBackend_SEQ_CG_ILU0
      : public ISTLBackend_SEQ_ILU0<Dune::CGSolver>
    {
    public:
      /*! \brief make a linear solver object

        \param[in] maxiter_ maximum number of iterations to do
        \param[in] verbose_ print messages if true
      */
      explicit ISTLBackend_SEQ_CG_ILU0 (unsigned maxiter_=5000, int verbose_=1)
        : ISTLBackend_SEQ_ILU0<Dune::CGSolver>(maxiter_, verbose_)
      {}
    };

    //! \brief Sequential BiCGStab solver with ILU0 preconditioner
    class ISTLBackend_SEQ_BCGS_ILUn
      : public ISTLBackend_SEQ_ILUn<Dune::BiCGSTABSolver>
    {
    public:
      /*! \brief make a linear solver object


        \param[in] n_ The number of levels to be used.
        \param[in] w_ The relaxation factor.
        \param[in] maxiter_ maximum number of iterations to do
        \param[in] verbose_ print messages if true
      */
      explicit ISTLBackend_SEQ_BCGS_ILUn (int n_, double w_=1.0, unsigned maxiter_=5000, int verbose_=1)
        : ISTLBackend_SEQ_ILUn<Dune::BiCGSTABSolver>(n_, w_, maxiter_, verbose_)
      {}
    };

    //! \brief Sequential congute gradient solver with ILU0 preconditioner
    class ISTLBackend_SEQ_CG_ILUn
      : public ISTLBackend_SEQ_ILUn<Dune::CGSolver>
    {
    public:
      /*! \brief make a linear solver object


        \param[in] n_ The number of levels to be used.
        \param[in] w_ The relaxation factor.
        \param[in] maxiter_ maximum number of iterations to do
        \param[in] verbose_ print messages if true
      */
      explicit ISTLBackend_SEQ_CG_ILUn (int n_, double w_=1.0, unsigned maxiter_=5000, int verbose_=1)
        : ISTLBackend_SEQ_ILUn<Dune::CGSolver>(n_, w_, maxiter_, verbose_)
      {}
    };

    /**
     * @brief Backend for sequential conjugate gradient solver with SSOR preconditioner.
     */
    class ISTLBackend_SEQ_CG_SSOR
      : public ISTLBackend_SEQ_Base<Dune::SeqSSOR, Dune::CGSolver>
    {
    public:
      /*! \brief make a linear solver object

        \param[in] maxiter_ maximum number of iterations to do
        \param[in] verbose_ print messages if true
      */
      explicit ISTLBackend_SEQ_CG_SSOR (unsigned maxiter_=5000, int verbose_=1)
        : ISTLBackend_SEQ_Base<Dune::SeqSSOR, Dune::CGSolver>(maxiter_, verbose_)
      {}
    };

    /**
     * @brief Backend using a MINRes solver preconditioned by SSOR.
     */
    class ISTLBackend_SEQ_MINRES_SSOR
      : public ISTLBackend_SEQ_Base<Dune::SeqSSOR, Dune::MINRESSolver>
    {
    public:
      /*! \brief make a linear solver object

        \param[in] maxiter_ maximum number of iterations to do
        \param[in] verbose_ print messages if true
      */
      explicit ISTLBackend_SEQ_MINRES_SSOR (unsigned maxiter_=5000, int verbose_=1)
        : ISTLBackend_SEQ_Base<Dune::SeqSSOR, Dune::MINRESSolver>(maxiter_, verbose_)
      {}
    };

    /**
     * @brief Backend for conjugate gradient solver with Jacobi preconditioner.
     */
    class ISTLBackend_SEQ_CG_Jac
      : public ISTLBackend_SEQ_Base<Dune::SeqJac, Dune::CGSolver>
    {
    public:
      /*! \brief make a linear solver object
        \param[in] maxiter_ maximum number of iterations to do
        \param[in] verbose_ print messages if true
      */
      explicit ISTLBackend_SEQ_CG_Jac (unsigned maxiter_=5000, int verbose_=1)
        : ISTLBackend_SEQ_Base<Dune::SeqJac, Dune::CGSolver>(maxiter_, verbose_)
      {}
    };

#if HAVE_SUPERLU || DOXYGEN
    /**
     * @brief Solver backend using SuperLU as a direct solver.
     */
    class ISTLBackend_SEQ_SuperLU
      : public SequentialNorm, public LinearResultStorage
    {
    public:
      /*! \brief make a linear solver object

        \param[in] verbose_ print messages if true
      */
      explicit ISTLBackend_SEQ_SuperLU (int verbose_=1)
        : verbose(verbose_)
      {}


      /*! \brief make a linear solver object

        \param[in] maxiter Maximum number of allowed steps (ignored)
        \param[in] verbose_ print messages if true
      */
      ISTLBackend_SEQ_SuperLU (int maxiter, int verbose_)
        : verbose(verbose_)
      {}

      /*! \brief solve the given linear system

        \param[in] A the given matrix
        \param[out] z the solution vector to be computed
        \param[in] r right hand side
        \param[in] reduction to be achieved
      */
      template<class M, class V, class W>
      void apply(M& A, V& z, W& r, typename Dune::template FieldTraits<typename W::ElementType >::real_type reduction)
      {
        using Backend::Native;
        using Backend::native;
        using ISTLM = Native<M>;
        Dune::SuperLU<ISTLM> solver(native(A), verbose);
        Dune::InverseOperatorResult stat;
        solver.apply(native(z), native(r), stat);
        res.converged  = stat.converged;
        res.iterations = stat.iterations;
        res.elapsed    = stat.elapsed;
        res.reduction  = stat.reduction;
        res.conv_rate  = stat.conv_rate;
      }

    private:
      int verbose;
    };
#endif // HAVE_SUPERLU || DOXYGEN

#if HAVE_SUITESPARSE_UMFPACK || DOXYGEN
    /**
     * @brief Solver backend using UMFPack as a direct solver.
     */
    class ISTLBackend_SEQ_UMFPack
      : public SequentialNorm, public LinearResultStorage
    {
    public:
      /*! \brief make a linear solver object

        \param[in] verbose_ print messages if true
      */
      explicit ISTLBackend_SEQ_UMFPack (int verbose_=1)
        : verbose(verbose_)
      {}


      /*! \brief make a linear solver object

        \param[in] maxiter Maximum number of allowed steps (ignored)
        \param[in] verbose_ print messages if true
      */
      ISTLBackend_SEQ_UMFPack (int maxiter, int verbose_)
        : verbose(verbose_)
      {}

      /*! \brief solve the given linear system

        \param[in] A the given matrix
        \param[out] z the solution vector to be computed
        \param[in] r right hand side
        \param[in] reduction to be achieved
      */
      template<class M, class V, class W>
      void apply(M& A, V& z, W& r, typename Dune::template FieldTraits<typename W::ElementType >::real_type reduction)
      {
        using Backend::native;
        using ISTLM = Backend::Native<M>;
        Dune::UMFPack<ISTLM> solver(native(A), verbose);
        Dune::InverseOperatorResult stat;
        solver.apply(native(z), native(r), stat);
        res.converged  = stat.converged;
        res.iterations = stat.iterations;
        res.elapsed    = stat.elapsed;
        res.reduction  = stat.reduction;
        res.conv_rate  = stat.conv_rate;
      }

    private:
      int verbose;
    };
#endif // HAVE_SUITESPARSE_UMFPACK || DOXYGEN

    //! Solver to be used for explicit time-steppers with (block-)diagonal mass matrix
    class ISTLBackend_SEQ_ExplicitDiagonal
      : public SequentialNorm, public LinearResultStorage
    {
    public:
      /*! \brief make a linear solver object
      */
      ISTLBackend_SEQ_ExplicitDiagonal ()
      {}

      /*! \brief solve the given linear system

        \param[in] A the given matrix
        \param[out] z the solution vector to be computed
        \param[in] r right hand side
        \param[in] reduction to be achieved
      */
      template<class M, class V, class W>
      void apply(M& A, V& z, W& r, typename Dune::template FieldTraits<typename W::ElementType >::real_type reduction)
      {
        using Backend::Native;
        Dune::SeqJac<Native<M>,
                     Native<V>,
                     Native<W>
                     > jac(Backend::native(A),1,1.0);
        jac.pre(z,r);
        jac.apply(z,r);
        jac.post(z);
        res.converged  = true;
        res.iterations = 1;
        res.elapsed    = 0.0;
        res.reduction  = reduction;
        res.conv_rate  = reduction; // pow(reduction,1.0/1)
      }
    };

    //! \} Sequential Solvers

    /**
     * @brief Class providing some statistics of the AMG solver.
     *
     */
    struct ISTLAMGStatistics
    {
      /**
       * @brief The needed for computing the parallel information and
       * for adapting the linear system.
       */
      double tprepare;
      /** @brief the number of levels in the AMG hierarchy. */
      int levels;
      /** @brief The time spent in solving the system (without building the hierarchy. */
      double tsolve;
      /** @brief The time needed for building the AMG hierarchy (coarsening). */
      double tsetup;
      /** @brief The number of iterations performed until convergence was reached. */
      int iterations;
      /** @brief True if a direct solver was used on the coarset level. */
      bool directCoarseLevelSolver;
    };

    template<class GO, template<class,class,class,int> class Preconditioner, template<class> class Solver,
              bool skipBlocksizeCheck = false>
    class ISTLBackend_SEQ_AMG : public LinearResultStorage
    {
      typedef typename GO::Traits::TrialGridFunctionSpace GFS;
      typedef typename GO::Traits::Jacobian M;
      typedef Backend::Native<M> MatrixType;
      typedef typename GO::Traits::Domain V;
      typedef Backend::Native<V> VectorType;
      typedef Preconditioner<MatrixType,VectorType,VectorType,1> Smoother;
      typedef Dune::MatrixAdapter<MatrixType,VectorType,VectorType> Operator;
      typedef typename Dune::Amg::SmootherTraits<Smoother>::Arguments SmootherArgs;
      typedef Dune::Amg::AMG<Operator,VectorType,Smoother> AMG;
      typedef Dune::Amg::Parameters Parameters;

    public:
      ISTLBackend_SEQ_AMG(unsigned maxiter_=5000, int verbose_=1,
                          bool reuse_=false, bool usesuperlu_=true)
        : maxiter(maxiter_), params(15,2000), verbose(verbose_),
          reuse(reuse_), firstapply(true), usesuperlu(usesuperlu_)
      {
        params.setDefaultValuesIsotropic(GFS::Traits::GridViewType::Traits::Grid::dimension);
        params.setDebugLevel(verbose_);
#if !HAVE_SUPERLU
        if (usesuperlu == true)
          {
            std::cout << "WARNING: You are using AMG without SuperLU!"
                      << " Please consider installing SuperLU,"
                      << " or set the usesuperlu flag to false"
                      << " to suppress this warning." << std::endl;
          }
#endif
      }

       /*! \brief set AMG parameters

        \param[in] params_ a parameter object of Type Dune::Amg::Parameters
      */
      void setparams(Parameters params_)
      {
        params = params_;
      }

      //! Set whether the AMG should be reused again during call to apply().
      void setReuse(bool reuse_)
      {
        reuse = reuse_;
      }

      //! Return whether the AMG is reused during call to apply()
      bool getReuse() const
      {
        return reuse;
      }

      /*! \brief compute global norm of a vector

        \param[in] v the given vector
      */
      typename V::ElementType norm (const V& v) const
      {
        return Backend::native(v).two_norm();
      }

      /*! \brief solve the given linear system

        \param[in] A the given matrix
        \param[out] z the solution vector to be computed
        \param[in] r right hand side
        \param[in] reduction to be achieved
      */
      void apply(M& A, V& z, V& r, typename Dune::template FieldTraits<typename V::ElementType >::real_type reduction)
      {
        Timer watch;
        MatrixType& mat = Backend::native(A);
        typedef Dune::Amg::CoarsenCriterion<Dune::Amg::SymmetricCriterion<MatrixType,
          Dune::Amg::FirstDiagonal> > Criterion;
        SmootherArgs smootherArgs;
        smootherArgs.iterations = 1;
        smootherArgs.relaxationFactor = 1;

        Criterion criterion(params);
        //only construct a new AMG if the matrix changes
        if (reuse==false || firstapply==true){
          oop.reset(new Operator(mat));
          amg.reset(new AMG(*oop, criterion, smootherArgs));
          firstapply = false;
          stats.tsetup = watch.elapsed();
          stats.levels = amg->maxlevels();
          stats.directCoarseLevelSolver=amg->usesDirectCoarseLevelSolver();
        }
        watch.reset();
        Dune::InverseOperatorResult stat;

        Solver<VectorType> solver(*oop,*amg,reduction,maxiter,verbose);
        solver.apply(Backend::native(z),Backend::native(r),stat);
        stats.tsolve= watch.elapsed();
        res.converged  = stat.converged;
        res.iterations = stat.iterations;
        res.elapsed    = stat.elapsed;
        res.reduction  = stat.reduction;
        res.conv_rate  = stat.conv_rate;
      }


      /**
       * @brief Get statistics of the AMG solver (no of levels, timings).
       * @return statistis of the AMG solver.
       */
      const ISTLAMGStatistics& statistics() const
      {
        return stats;
      }

    private:
      unsigned maxiter;
      Parameters params;
      int verbose;
      bool reuse;
      bool firstapply;
      bool usesuperlu;
      std::shared_ptr<Operator> oop;
      std::shared_ptr<AMG> amg;
      ISTLAMGStatistics stats;
    };

    //! \addtogroup PDELab_seqsolvers Sequential Solvers
    //! \{

    /**
     * @brief Sequential conjugate gradient solver preconditioned with AMG smoothed by SSOR
     * @tparam GO The type of the grid operator
     * (or the fakeGOTraits class for the old grid operator space).
     */
    template<class GO>
    class ISTLBackend_SEQ_CG_AMG_SSOR
      : public ISTLBackend_SEQ_AMG<GO, Dune::SeqSSOR, Dune::CGSolver>
    {

    public:
      /**
       * @brief Constructor
       * @param maxiter_ The maximum number of iterations allowed.
       * @param verbose_ The verbosity level to use.
       * @param reuse_ Set true, if the Matrix to be used is always identical
       * (AMG aggregation is then only performed once).
       * @param usesuperlu_ Set false, to suppress the no SuperLU warning
       */
      ISTLBackend_SEQ_CG_AMG_SSOR(unsigned maxiter_=5000, int verbose_=1,
                                  bool reuse_=false, bool usesuperlu_=true)
        : ISTLBackend_SEQ_AMG<GO, Dune::SeqSSOR, Dune::CGSolver>
          (maxiter_, verbose_, reuse_, usesuperlu_)
      {}
    };

    /**
     * @brief Sequential BiCGStab solver preconditioned with AMG smoothed by SSOR
     * @tparam GO The type of the grid operator
     * (or the fakeGOTraits class for the old grid operator space).
     */
    template<class GO>
    class ISTLBackend_SEQ_BCGS_AMG_SSOR
      : public ISTLBackend_SEQ_AMG<GO, Dune::SeqSSOR, Dune::BiCGSTABSolver>
    {

    public:
      /**
       * @brief Constructor
       * @param maxiter_ The maximum number of iterations allowed.
       * @param verbose_ The verbosity level to use.
       * @param reuse_ Set true, if the Matrix to be used is always identical
       * (AMG aggregation is then only performed once).
       * @param usesuperlu_ Set false, to suppress the no SuperLU warning
       */
      ISTLBackend_SEQ_BCGS_AMG_SSOR(unsigned maxiter_=5000, int verbose_=1,
                                    bool reuse_=false, bool usesuperlu_=true)
        : ISTLBackend_SEQ_AMG<GO, Dune::SeqSSOR, Dune::BiCGSTABSolver>
          (maxiter_, verbose_, reuse_, usesuperlu_)
      {}
    };

    /**
     * @brief Sequential BiCGSTAB solver preconditioned with AMG smoothed by SOR
     * @tparam GO The type of the grid operator
     * (or the fakeGOTraits class for the old grid operator space).
     */
    template<class GO>
    class ISTLBackend_SEQ_BCGS_AMG_SOR
      : public ISTLBackend_SEQ_AMG<GO, Dune::SeqSOR, Dune::BiCGSTABSolver>
    {

    public:
      /**
       * @brief Constructor
       * @param maxiter_ The maximum number of iterations allowed.
       * @param verbose_ The verbosity level to use.
       * @param reuse_ Set true, if the Matrix to be used is always identical
       * (AMG aggregation is then only performed once).
       * @param usesuperlu_ Set false, to suppress the no SuperLU warning
       */
      ISTLBackend_SEQ_BCGS_AMG_SOR(unsigned maxiter_=5000, int verbose_=1,
                                   bool reuse_=false, bool usesuperlu_=true)
        : ISTLBackend_SEQ_AMG<GO, Dune::SeqSOR, Dune::BiCGSTABSolver>
          (maxiter_, verbose_, reuse_, usesuperlu_)
      {}
    };

    /**
     * @brief Sequential Loop solver preconditioned with AMG smoothed by SSOR
     * @tparam GO The type of the grid operator
     * (or the fakeGOTraits class for the old grid operator space).
     */
    template<class GO>
    class ISTLBackend_SEQ_LS_AMG_SSOR
      : public ISTLBackend_SEQ_AMG<GO, Dune::SeqSSOR, Dune::LoopSolver>
    {

    public:
      /**
       * @brief Constructor
       * @param maxiter_ The maximum number of iterations allowed.
       * @param verbose_ The verbosity level to use.
       * @param reuse_ Set true, if the Matrix to be used is always identical
       * (AMG aggregation is then only performed once).
       * @param usesuperlu_ Set false, to suppress the no SuperLU warning
       */
      ISTLBackend_SEQ_LS_AMG_SSOR(unsigned maxiter_=5000, int verbose_=1,
                                  bool reuse_=false, bool usesuperlu_=true)
        : ISTLBackend_SEQ_AMG<GO, Dune::SeqSSOR, Dune::LoopSolver>
          (maxiter_, verbose_, reuse_, usesuperlu_)
      {}
    };

    /**
     * @brief Sequential Loop solver preconditioned with AMG smoothed by SOR
     * @tparam GO The type of the grid operator
     * (or the fakeGOTraits class for the old grid operator space).
     */
    template<class GO>
    class ISTLBackend_SEQ_LS_AMG_SOR
      : public ISTLBackend_SEQ_AMG<GO, Dune::SeqSOR, Dune::LoopSolver>
    {

    public:
      /**
       * @brief Constructor
       * @param maxiter_ The maximum number of iterations allowed.
       * @param verbose_ The verbosity level to use.
       * @param reuse_ Set true, if the Matrix to be used is always identical
       * (AMG aggregation is then only performed once).
       * @param usesuperlu_ Set false, to suppress the no SuperLU warning
       */
      ISTLBackend_SEQ_LS_AMG_SOR(unsigned maxiter_=5000, int verbose_=1,
                                 bool reuse_=false, bool usesuperlu_=true)
        : ISTLBackend_SEQ_AMG<GO, Dune::SeqSOR, Dune::LoopSolver>
          (maxiter_, verbose_, reuse_, usesuperlu_)
      {}
    };

    /** \brief Linear solver backend for Restarted GMRes
        preconditioned with ILU(0)

    */

    class ISTLBackend_SEQ_GMRES_ILU0
      : public SequentialNorm, public LinearResultStorage
    {
    public :

      /** \brief make linear solver object

          \param[in] restart_ number of iterations when GMRes has to be restarted
          \param[in] maxiter_ maximum number of iterations to do
          \param[in] verbose_ print messages if true
      */
      explicit ISTLBackend_SEQ_GMRES_ILU0(int restart_ = 200, int maxiter_ = 5000, int verbose_ = 1)
        : restart(restart_), maxiter(maxiter_), verbose(verbose_)
      {}

      /** \brief solve the given linear system

          \param[in] A the given matrix
          \param[out] z the solution vector to be computed
          \param[in] r right hand side
          \param[in] reduction to be achieved
      */
      template<class M, class V, class W>
      void apply(M& A, V& z, W& r, typename Dune::template FieldTraits<typename W::ElementType>::real_type reduction)
      {
        using Backend::Native;
        using Backend::native;
        Dune::MatrixAdapter<
          Native<M>,
          Native<V>,
          Native<W>
          > opa(native(A));
        Dune::SeqILU0<
          Native<M>,
          Native<V>,
          Native<W>
          > ilu0(native(A), 1.0);
        Dune::RestartedGMResSolver<Native<V>> solver(opa,ilu0,reduction,restart,maxiter,verbose);
        Dune::InverseOperatorResult stat;
        solver.apply(native(z), native(r), stat);
        res.converged  = stat.converged;
        res.iterations = stat.iterations;
        res.elapsed    = stat.elapsed;
        res.reduction  = stat.reduction;
        res.conv_rate  = stat.conv_rate;
      }

    private :
      int restart, maxiter, verbose;
    };

    //! \} group Sequential Solvers
    //! \} group Backend

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_BACKEND_ISTL_SEQISTLSOLVERBACKEND_HH
