// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_SEQISTLSOLVERBACKEND_HH
#define DUNE_SEQISTLSOLVERBACKEND_HH

#include <dune/common/deprecated.hh>
#include <dune/common/mpihelper.hh>

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

#include "../constraints/constraints.hh"
#include "../gridfunctionspace/genericdatahandle.hh"
#include "../newton/newton.hh"
#include "istlvectorbackend.hh"
#include "parallelistlhelper.hh"

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

      enum {category=Dune::SolverCategory::sequential};

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

    private:
      GOS& gos;
    };

    //==============================================================================
    // Here we add some standard linear solvers conforming to the linear solver
    // interface required to solve linear and nonlinear problems.
    //==============================================================================

    struct SequentialNorm
    {/*! \brief compute global norm of a vector

        \param[in] v the given vector
      */
      template<class V>
      typename V::ElementType norm(const V& v) const
      {
        return v.base().two_norm();
      }
    };

    class LinearResultStorage
    {
    public:
      /*! \brief Return access to result data */
      const Dune::PDELab::LinearSolverResult<double>& result() const
      {
        return res;
      }

    protected:
      Dune::PDELab::LinearSolverResult<double> res;
    };
    

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
      void apply(M& A, V& z, W& r, typename W::ElementType reduction)
      {
        Dune::MatrixAdapter<typename M::BaseT, 
                            typename V::BaseT, 
                            typename W::BaseT> opa(A);
        Preconditioner<typename M::BaseT, 
                            typename V::BaseT, 
                       typename W::BaseT,1> ssor(A, 3, 1.0);
        Solver<typename V::BaseT> solver(opa, ssor, reduction, maxiter, verbose);
        Dune::InverseOperatorResult stat;
        solver.apply(z, r, stat);
        res.converged  = stat.converged;
        res.iterations = stat.iterations;
        res.elapsed    = stat.elapsed;
        res.reduction  = stat.reduction;
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
      void apply(M& A, V& z, W& r, typename W::ElementType reduction)
      {
        Dune::MatrixAdapter<typename M::BaseT, 
                            typename V::BaseT, 
                            typename W::BaseT> opa(A);
        Dune::SeqILU0<typename M::BaseT, 
                      typename V::BaseT, 
                      typename W::BaseT> ilu0(A, 1.0);
        Solver<typename V::BaseT> solver(opa, ilu0, reduction, maxiter, verbose);
        Dune::InverseOperatorResult stat;
        solver.apply(z, r, stat);
        res.converged  = stat.converged;
        res.iterations = stat.iterations;
        res.elapsed    = stat.elapsed;
        res.reduction  = stat.reduction;
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
      explicit ISTLBackend_SEQ_ILUn (int n, double w, unsigned maxiter_=5000, int verbose_=1)
        : n_(n), w_(w), maxiter(maxiter_), verbose(verbose_)
       {}
      /*! \brief solve the given linear system

        \param[in] A the given matrix
        \param[out] z the solution vector to be computed
        \param[in] r right hand side
        \param[in] reduction to be achieved
      */
      template<class M, class V, class W>
      void apply(M& A, V& z, W& r, typename W::ElementType reduction)
      {
        Dune::MatrixAdapter<M,V,W> opa(A);
        Dune::SeqILUn<typename M::BaseT,V,W> ilu0(A.base(), n_, w_);
        Solver<V> solver(opa, ilu0, reduction, maxiter, verbose);
        Dune::InverseOperatorResult stat;
        solver.apply(z, r, stat);
        res.converged  = stat.converged;
        res.iterations = stat.iterations;
        res.elapsed    = stat.elapsed;
        res.reduction  = stat.reduction;
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
     * @brief Backend for sequential BiCGSTAB solver with Jacobi preconditioner.
     */
    class ISTLBackend_SEQ_BCGS_Jac
      : public ISTLBackend_SEQ_Base<Dune::SeqJac, Dune::BiCGSTABSolver>
    {
    public:
      /*! \brief make a linear solver object
        \param[in] maxiter maximum number of iterations to do
        \param[in] verbose print messages if true
      */
      explicit ISTLBackend_SEQ_BCGS_Jac (unsigned maxiter_=5000, bool verbose_=true)
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

    //! Sequential BiCGStab solver with ILU0 preconditioner
    class ISTLBackend_SEQ_BCGS_ILUn
      : public ISTLBackend_SEQ_ILUn<Dune::BiCGSTABSolver>
    {
    public:
      /*! \brief make a linear solver object

        
        \param[in] n The number of levels to be used.
        \param[in] w The relaxation factor.
        \param[in] maxiter_ maximum number of iterations to do
        \param[in] verbose_ print messages if true
      */
      explicit ISTLBackend_SEQ_BCGS_ILUn (int n, double w=1.0, unsigned maxiter_=5000, int verbose_=1)
        : ISTLBackend_SEQ_ILUn<Dune::BiCGSTABSolver>(n, w, maxiter_, verbose_)
      {}
    }; 

    //! Sequential congute gradient solver with ILU0 preconditioner
    class ISTLBackend_SEQ_CG_ILUn
      : public ISTLBackend_SEQ_ILUn<Dune::CGSolver>
    {
    public:
      /*! \brief make a linear solver object

        
        \param[in] n The number of levels to be used.
        \param[in] w The relaxation factor.
        \param[in] maxiter_ maximum number of iterations to do
        \param[in] verbose_ print messages if true
      */
      explicit ISTLBackend_SEQ_CG_ILUn (int n, double w=1.0, unsigned maxiter_=5000, int verbose_=1)
        : ISTLBackend_SEQ_ILUn<Dune::CGSolver>(n, w, maxiter_, verbose_)
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
     * @brief Backend for conjugate gradient solver with Jacobi preconditioner.
     */
    class ISTLBackend_SEQ_CG_Jac
      : public ISTLBackend_SEQ_Base<Dune::SeqJac, Dune::CGSolver>
    {
    public:
      /*! \brief make a linear solver object
        \param[in] maxiter maximum number of iterations to do
        \param[in] verbose print messages if true
      */
      explicit ISTLBackend_SEQ_CG_Jac (unsigned maxiter_=5000, bool verbose_=true)
        : ISTLBackend_SEQ_Base<Dune::SeqJac, Dune::CGSolver>(maxiter_, verbose_)
      {}
    };

#if HAVE_SUPERLU
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

      /*! \brief solve the given linear system

        \param[in] A the given matrix
        \param[out] z the solution vector to be computed
        \param[in] r right hand side
        \param[in] reduction to be achieved
      */
      template<class M, class V, class W>
      void apply(M& A, V& z, W& r, typename W::ElementType reduction)
      {
        typedef typename M::BaseT ISTLM;
        Dune::SuperLU<ISTLM> solver(A, verbose);
        Dune::InverseOperatorResult stat;
        solver.apply(z, r, stat);
        res.converged  = stat.converged;
        res.iterations = stat.iterations;
        res.elapsed    = stat.elapsed;
        res.reduction  = stat.reduction;
      }

    private:
      int verbose;
    };
#endif // HAVE_SUPERLU

    //! Solver to be used for explicit time-steppers with (block-)diagonal mass matrix
    class ISTLBackend_SEQ_ExplicitDiagonal
      : public SequentialNorm, public LinearResultStorage
    {
    public:
      /*! \brief make a linear solver object
      */
      explicit ISTLBackend_SEQ_ExplicitDiagonal ()
      {}

      /*! \brief solve the given linear system

        \param[in] A the given matrix
        \param[out] z the solution vector to be computed
        \param[in] r right hand side
        \param[in] reduction to be achieved
      */
      template<class M, class V, class W>
      void apply(M& A, V& z, W& r, typename W::ElementType reduction)
      {
        Dune::SeqJac<typename M::BaseT, 
                     typename V::BaseT, 
                     typename W::BaseT> jac(A,1,1.0);
        jac.pre(z,r);
        jac.apply(z,r);
        jac.post(z);
        res.converged  = true;
        res.iterations = 1;
        res.elapsed    = 0.0;
        res.reduction  = reduction;
      }
    };

    //! \} Sequential Solvers

    template<class GFS, template<class,class,class,int> class SMI, template<class> class SOI>
    class ISTLBackend_SEQ_AMG
    {

    public:
      ISTLBackend_SEQ_AMG(int smoothsteps=2,
                          unsigned maxiter_=5000, int verbose_=1)
        : maxiter(maxiter_), steps(smoothsteps), verbose(verbose_)
      {}


      /*! \brief compute global norm of a vector

        \param[in] v the given vector
      */
      template<class V>
      typename V::ElementType norm (const V& v) const
      {
        return v.base().two_norm();
      }

      /*! \brief solve the given linear system

        \param[in] A the given matrix
        \param[out] z the solution vector to be computed
        \param[in] r right hand side
        \param[in] reduction to be achieved
      */
      template<class M, class V>
      void apply(M& A, V& z, V& r, typename V::ElementType reduction)
      {
        typedef typename M::BaseT MatrixType;
        typedef typename BlockProcessor<GFS>::template AMGVectorTypeSelector<V>::Type
          VectorType;
        MatrixType& mat=A.base();
        typedef Dune::Amg::CoarsenCriterion<Dune::Amg::SymmetricCriterion<MatrixType,
          Dune::Amg::FirstDiagonal> >
          Criterion;
        typedef SMI<MatrixType,VectorType,VectorType,1> Smoother;
        typedef typename Dune::Amg::SmootherTraits<Smoother>::Arguments
          SmootherArgs;
        typedef Dune::MatrixAdapter<MatrixType,VectorType,VectorType> Operator;
        typedef Dune::Amg::AMG<Operator,VectorType,Smoother> AMG;
        SmootherArgs smootherArgs;
        smootherArgs.iterations = 1;
        smootherArgs.relaxationFactor = 1;

        Criterion criterion(15,2000);
        criterion.setDebugLevel(verbose?3:0);
        Operator oop(mat);
        AMG amg=AMG(oop, criterion, smootherArgs, 1, steps, steps);

        Dune::InverseOperatorResult stat;

        SOI<VectorType> solver(oop,amg,reduction,maxiter,verbose);
        solver.apply(BlockProcessor<GFS>::getVector(z),BlockProcessor<GFS>::getVector(r),stat);
        res.converged  = stat.converged;
        res.iterations = stat.iterations;
        res.elapsed    = stat.elapsed;
        res.reduction  = stat.reduction;
      }

      /*! \brief Return access to result data */
      const Dune::PDELab::LinearSolverResult<double>& result() const
      {
        return res;
      }

    private:
      Dune::PDELab::LinearSolverResult<double> res;
      unsigned maxiter;
      int steps;
      int verbose;
    };

    //! \addtogroup PDELab_seqsolvers Sequential Solvers
    //! \{

    /**
     * @brief Sequential conjugate gradient solver preconditioned with AMG smoothed by SSOR
     * @tparam GFS The type of the grid functions space.
     */
    template<class GFS>
    class ISTLBackend_SEQ_CG_AMG_SSOR
      : public ISTLBackend_SEQ_AMG<GFS, Dune::SeqSSOR, Dune::CGSolver>
    {

    public:
      /**
       * @brief Constructor
       * @param smoothsteps_ The number of steps to use for both pre and post smoothing.
       * @param maxiter_ The maximum number of iterations allowed.
       * @param verbose_ The verbosity level to use.
       */
      ISTLBackend_SEQ_CG_AMG_SSOR(int smoothsteps_=2,
                                  unsigned maxiter_=5000, int verbose_=1)
        : ISTLBackend_SEQ_AMG<GFS, Dune::SeqSSOR, Dune::CGSolver>(smoothsteps_, maxiter_,verbose_)
      {}
    };

    /**
     * @brief Sequential BiCGStab solver preconditioned with AMG smoothed by SSOR
     * @tparam GFS The type of the grid functions space.
     */
    template<class GFS>
    class ISTLBackend_SEQ_BCGS_AMG_SSOR
      : public ISTLBackend_SEQ_AMG<GFS, Dune::SeqSSOR, Dune::BiCGSTABSolver>
    {

    public:
      /**
       * @brief Constructor
       * @param smoothsteps_ The number of steps to use for both pre and post smoothing.
       * @param maxiter_ The maximum number of iterations allowed.
       * @param verbose_ The verbosity level to use.
       */
      ISTLBackend_SEQ_BCGS_AMG_SSOR(int smoothsteps_=2,
                                    unsigned maxiter_=5000, int verbose_=1)
        : ISTLBackend_SEQ_AMG<GFS, Dune::SeqSSOR, Dune::BiCGSTABSolver>(smoothsteps_, maxiter_, verbose_)
      {}
    };
    
    /**
     * @brief Sequential BiCGSTAB solver preconditioned with AMG smoothed by SOR     * @tparam GFS The type of the grid functions space.
     */
    template<class GFS>
    class ISTLBackend_SEQ_BCGS_AMG_SOR
      : public ISTLBackend_SEQ_AMG<GFS, Dune::SeqSOR, Dune::BiCGSTABSolver>
    {

    public:
      /**
       * @brief Constructor
       * @param smoothsteps_ The number of steps to use for both pre and post smoothing.
       * @param maxiter_ The maximum number of iterations allowed.
       * @param verbose_ The verbosity level to use.
       */
      ISTLBackend_SEQ_BCGS_AMG_SOR(int smoothsteps_=2,
                                    unsigned maxiter_=5000, int verbose_=1)
        : ISTLBackend_SEQ_AMG<GFS, Dune::SeqSOR, Dune::BiCGSTABSolver>(smoothsteps_, maxiter_, verbose_)
      {}
    };

    /**
     * @brief Sequential Loop solver preconditioned with AMG smoothed by SSOR
     * @tparam GFS The type of the grid functions space.
     */
    template<class GFS>
    class ISTLBackend_SEQ_LS_AMG_SSOR
      : public ISTLBackend_SEQ_AMG<GFS, Dune::SeqSSOR, Dune::LoopSolver>
    {

    public:
      /**
       * @brief Constructor
       * @param smoothsteps_ The number of steps to use for both pre and post smoothing.
       * @param maxiter_ The maximum number of iterations allowed.
       * @param verbose_ The verbosity level to use.
       */
      ISTLBackend_SEQ_LS_AMG_SSOR(int smoothsteps_=2,
                                    unsigned maxiter_=5000, int verbose_=1)
        : ISTLBackend_SEQ_AMG<GFS, Dune::SeqSSOR, Dune::LoopSolver>(smoothsteps_, maxiter_, verbose_)
      {}
    };

    /**
     * @brief Sequential Loop solver preconditioned with AMG smoothed by SOR
     * @tparam GFS The type of the grid functions space.
     */
    template<class GFS>
    class ISTLBackend_SEQ_LS_AMG_SOR
      : public ISTLBackend_SEQ_AMG<GFS, Dune::SeqSOR, Dune::LoopSolver>
    {

    public:
      /**
       * @brief Constructor
       * @param smoothsteps_ The number of steps to use for both pre and post smoothing.
       * @param maxiter_ The maximum number of iterations allowed.
       * @param verbose_ The verbosity level to use.
       */
      ISTLBackend_SEQ_LS_AMG_SOR(int smoothsteps_=2,
                                    unsigned maxiter_=5000, int verbose_=1)
        : ISTLBackend_SEQ_AMG<GFS, Dune::SeqSOR, Dune::LoopSolver>(smoothsteps_, maxiter_, verbose_)
      {}
    };

    //! \} group Sequential Solvers
    //! \} group Backend

  } // namespace PDELab
} // namespace Dune

#endif
