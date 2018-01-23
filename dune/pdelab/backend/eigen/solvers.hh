// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_EIGEN_SOLVERS_HH
#define DUNE_PDELAB_BACKEND_EIGEN_SOLVERS_HH

#include <dune/common/timer.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/pdelab/constraints/common/constraints.hh>

#include "../solver.hh"

#if HAVE_EIGEN

#include <Eigen/Eigen>
#include <Eigen/Sparse>

namespace Dune {
  namespace PDELab {

    //==============================================================================
    // Here we add some standard linear solvers conforming to the linear solver
    // interface required to solve linear and nonlinear problems.
    //==============================================================================

  template<class PreconditionerImp>
    class EigenBackend_BiCGSTAB_Base
      : public LinearResultStorage
    {
    public:
      /*! \brief make a linear solver object

        \param[in] maxiter_ maximum number of iterations to do
        \param[in] verbose_ print messages if true
      */
      explicit EigenBackend_BiCGSTAB_Base(unsigned maxiter_=5000)
        : maxiter(maxiter_)
      {}

      /*! \brief solve the given linear system

        \param[in] A the given matrix
        \param[out] z the solution vector to be computed
        \param[in] r right hand side
        \param[in] reduction to be achieved
      */
      template<class M, class V, class W>
      void apply(M& A, V& z, W& r, typename W::field_type reduction)
      {
        using Backend::Native;
        using Backend::native;
        // PreconditionerImp preconditioner;
        using Mat = Native<M>;
        ::Eigen::BiCGSTAB<Mat, PreconditionerImp> solver;
        solver.setMaxIterations(maxiter);
        solver.setTolerance(reduction);
        Dune::Timer watch;
        watch.reset();
        solver.compute(native(A));
        native(z) = solver.solve(native(r));
        double elapsed = watch.elapsed();

        res.converged  = solver.info() == ::Eigen::ComputationInfo::Success;
        res.iterations = solver.iterations();
        res.elapsed    = elapsed;
        res.reduction  = solver.error();
        res.conv_rate  = 0;
      }

    public:
      template<class V>
      typename Dune::template FieldTraits<typename V::ElementType >::real_type norm(const V& v) const
      {
        return Backend::native(v).norm();
      }

    private:
      unsigned maxiter;
      int verbose;
    };

    class EigenBackend_BiCGSTAB_IILU
      : public EigenBackend_BiCGSTAB_Base<::Eigen::IncompleteLUT<double> >
    {
    public:
        explicit EigenBackend_BiCGSTAB_IILU(unsigned maxiter_=5000)
          : EigenBackend_BiCGSTAB_Base(maxiter_)
        {
        }
    };

    class EigenBackend_BiCGSTAB_Diagonal
      : public EigenBackend_BiCGSTAB_Base<::Eigen::DiagonalPreconditioner<double> >
    {
    public:
        explicit EigenBackend_BiCGSTAB_Diagonal(unsigned maxiter_=5000)
          : EigenBackend_BiCGSTAB_Base(maxiter_)
        {}
    };

    template< class Preconditioner, int UpLo >
      class EigenBackend_CG_Base
      : public SequentialNorm, public LinearResultStorage
    {
    public:
      /*! \brief make a linear solver object

        \param[in] maxiter_ maximum number of iterations to do
        \param[in] verbose_ print messages if true
        */
      explicit EigenBackend_CG_Base(unsigned maxiter_=5000)
        : maxiter(maxiter_)
      {}

      /*! \brief solve the given linear system

        \param[in] A the given matrix
        \param[out] z the solution vector to be computed
        \param[in] r right hand side
        \param[in] reduction to be achieved
        */
      template<class M, class V, class W>
      void apply(M& A, V& z, W& r, typename W::field_type reduction)
      {
        using Backend::native;
        using Mat = Backend::Native<M>;
        ::Eigen::ConjugateGradient<Mat, UpLo, Preconditioner> solver;
        solver.setMaxIterations(maxiter);
        solver.setTolerance(reduction);
        Dune::Timer watch;
        watch.reset();
        solver.compute(native(A));
        native(z) = solver.solve(native(r));
        double elapsed = watch.elapsed();


        res.converged  = solver.info() == ::Eigen::ComputationInfo::Success;
        res.iterations = solver.iterations();
        res.elapsed    = elapsed;
        res.reduction  = solver.error();
        res.conv_rate  = 0;
      }

    private:
      unsigned maxiter;
      int verbose;
    };


    class EigenBackend_CG_IILU_Up
      : public EigenBackend_CG_Base<::Eigen::IncompleteLUT<double>, ::Eigen::Upper >
    {
    public:
        explicit EigenBackend_CG_IILU_Up(unsigned maxiter_=5000)
          : EigenBackend_CG_Base(maxiter_)
        {}
    };

    class EigenBackend_CG_Diagonal_Up
      : public EigenBackend_CG_Base<::Eigen::DiagonalPreconditioner<double>, ::Eigen::Upper >
    {
    public:
        explicit EigenBackend_CG_Diagonal_Up(unsigned maxiter_=5000)
          : EigenBackend_CG_Base(maxiter_)
        {}
    };

    class EigenBackend_CG_IILU_Lo
      : public EigenBackend_CG_Base<::Eigen::IncompleteLUT<double>, ::Eigen::Lower >
    {
    public:
        explicit EigenBackend_CG_IILU_Lo(unsigned maxiter_=5000)
          : EigenBackend_CG_Base(maxiter_)
        {}
    };

    class EigenBackend_CG_Diagonal_Lo
      : public EigenBackend_CG_Base<::Eigen::DiagonalPreconditioner<double>, ::Eigen::Lower >
    {
    public:
        explicit EigenBackend_CG_Diagonal_Lo(unsigned maxiter_=5000)
          : EigenBackend_CG_Base(maxiter_)
        {}
    };

#if EIGEN_VERSION_AT_LEAST(3,2,2)
    template<template<class, int, class> class Solver, int UpLo>
#else
    template<template<class, int> class Solver, int UpLo>
#endif
      class EigenBackend_SPD_Base
      : public SequentialNorm, public LinearResultStorage
    {
    public:
      /*! \brief make a linear solver object

        \param[in] maxiter_ maximum number of iterations to do
        \param[in] verbose_ print messages if true
        */
      explicit EigenBackend_SPD_Base()
      {}

      /*! \brief solve the given linear system

        \param[in] A the given matrix
        \param[out] z the solution vector to be computed
        \param[in] r right hand side
        \param[in] reduction to be achieved
        */
      template<class M, class V, class W>
      void apply(M& A, V& z, W& r, typename W::field_type reduction)
      {
        using Backend::native;
        using Mat = Backend::Native<M>;
#if EIGEN_VERSION_AT_LEAST(3,2,2)
        // use the approximate minimum degree algorithm for the ordering in
        // the solver. This reproduces the default ordering for the Cholesky
        // type solvers.
        Solver<Mat, UpLo, ::Eigen::AMDOrdering<typename Mat::Index> > solver;
#else
        Solver<Mat, UpLo> solver;
#endif
        Dune::Timer watch;
        watch.reset();
        solver.compute(native(A));
        native(z) = solver.solve(native(r));
        double elapsed = watch.elapsed();

        res.converged  = solver.info() == ::Eigen::ComputationInfo::Success;
        res.iterations = solver.iterations();
        res.elapsed    = elapsed;
        res.reduction  = solver.error();
        res.conv_rate  = 0;
      }
    private:
    };

    class EigenBackend_SimplicialCholesky_Up
      : public EigenBackend_SPD_Base<::Eigen::SimplicialCholesky, ::Eigen::Upper >
    {
    public:
        explicit EigenBackend_SimplicialCholesky_Up()
        {}
    };

    class EigenBackend_SimplicialCholesky_Lo
      : public EigenBackend_SPD_Base<::Eigen::SimplicialCholesky, ::Eigen::Lower >
    {
    public:
        explicit EigenBackend_SimplicialCholesky_Lo()
        {}
    };

/*    class EigenBackend_SimplicialLDLt_Up
 *      : public EigenBackend_SPD_Base<::Eigen::SimplicialLDLt, ::Eigen::Upper >
 *    {
 *    public:
 *        explicit EigenBackend_SimplicialLDLt_Up()
 *        {}
 *    };

 *    class EigenBackend_SimplicialLDLt_Lo
 *      : public EigenBackend_SPD_Base<::Eigen::SimplicialLDLt, ::Eigen::Lower >
 *    {
 *    public:
 *        explicit EigenBackend_SimplicialLDLt_Lo()
 *        {}
 *    };

 *    class EigenBackend_SimplicialLLt_Up
 *      : public EigenBackend_SPD_Base<::Eigen::SimplicialLLt, ::Eigen::Upper >
 *    {
 *    public:
 *        explicit EigenBackend_SimplicialLLt_Up()
 *        {}
 *    };

 *    class EigenBackend_SimplicialLLt_Lo
 *      : public EigenBackend_SPD_Base<::Eigen::SimplicialLLt, ::Eigen::Lower >
 *    {
 *    public:
 *        explicit EigenBackend_SimplicialLLt_Lo()
 *        {}
 *    };

 *    class EigenBackend_Cholmod_Up
 *      : public EigenBackend_SPD_Base<::Eigen::CholmodDecomposition, ::Eigen::Upper >
 *    {
 *    public:
 *        explicit EigenBackend_Cholmod_Up()
 *        {}
 *    };

 *    class EigenBackend_Cholmod_Lo
 *      : public EigenBackend_SPD_Base<::Eigen::CholmodDecomposition, ::Eigen::Lower >
 *    {
 *    public:
 *        explicit EigenBackend_Cholmod_Lo()
 *        {}
 *    };*/

    class EigenBackend_LeastSquares
      : public SequentialNorm, public LinearResultStorage
    {
    private:
      const static unsigned int defaultFlag = ::Eigen::ComputeThinU | ::Eigen::ComputeThinV;
    public:
      /*! \brief make a linear solver object

        \param[in] maxiter_ maximum number of iterations to do
        \param[in] verbose_ print messages if true
        */
      explicit EigenBackend_LeastSquares(unsigned int flags = defaultFlag)
        : flags_(flags)
      {}

      /*! \brief solve the given linear system

        \param[in] A the given matrix
        \param[out] z the solution vector to be computed
        \param[in] r right hand side
        \param[in] reduction to be achieved
        */
      template<class M, class V, class W>
      void apply(M& A, V& z, W& r, typename W::field_type reduction)
      {
        Dune::Timer watch;
        watch.reset();
        using Backend::native;
        using Mat = Backend::Native<M>;
        ::Eigen::JacobiSVD<Mat, ::Eigen::ColPivHouseholderQRPreconditioner> solver(A, flags_);
        native(z) = solver.solve(native(r));
        double elapsed = watch.elapsed();

        res.converged  = solver.info() == ::Eigen::ComputationInfo::Success;
        res.iterations = solver.iterations();
        res.elapsed    = elapsed;
        res.reduction  = solver.error();
        res.conv_rate  = 0;
      }
    private:
      unsigned int flags_;
    };

  } // namespace PDELab
} // namespace Dune

#endif // HAVE_EIGEN

#endif // DUNE_PDELAB_BACKEND_EIGEN_SOLVERS_HH
