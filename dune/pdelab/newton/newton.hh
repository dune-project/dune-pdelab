// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_NEWTON_NEWTON_HH
#define DUNE_PDELAB_NEWTON_NEWTON_HH

#include <iostream>
#include <iomanip>
#include <cmath>
#include <memory>

#include <type_traits>

#include <math.h>

#include <dune/common/debugstream.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/ios_state.hh>
#include <dune/common/timer.hh>
#include <dune/common/parametertree.hh>

#include <dune/pdelab/backend/solver.hh>

namespace Dune
{
  namespace PDELab
  {
  namespace Impl
  {
  // Some SFinae magic to execute setReuse(bool) on a backend
  template<typename T1, typename = void>
  struct HasSetReuse
  : std::false_type
  {};

  template<typename T>
  struct HasSetReuse<T, decltype(std::declval<T>().setReuse(true), void())>
  : std::true_type
  {};

  template<typename T>
  inline void setLinearSystemReuse(T& solver_backend, bool reuse, std::true_type)
  {
    if (!solver_backend.getReuse() && reuse)
      dwarn << "WARNING: Newton needed to override your choice to reuse the linear system in order to work!" << std::endl;
    solver_backend.setReuse(reuse);
  }

  template<typename T>
  inline void setLinearSystemReuse(T&, bool, std::false_type)
  {}

  template<typename T>
  inline void setLinearSystemReuse(T& solver_backend, bool reuse)
  {
    setLinearSystemReuse(solver_backend, reuse, HasSetReuse<T>());
  }
  }

    // Exception classes used in NewtonSolver
    class NewtonError : public Exception {};
    class NewtonDefectError : public NewtonError {};
    class NewtonLinearSolverError : public NewtonError {};
    class NewtonLineSearchError : public NewtonError {};
    class NewtonNotConverged : public NewtonError {};

    // Status information of Newton's method
    template<class RFType>
    struct NewtonResult : LinearSolverResult<RFType>
    {
      RFType first_defect;       // the first defect
      RFType defect;             // the final defect
      double assembler_time;     // Cumulative time for matrix assembly
      double linear_solver_time; // Cumulative time for linear solver
      int linear_solver_iterations; // Total number of linear iterations

      NewtonResult() :
        first_defect(0.0), defect(0.0), assembler_time(0.0), linear_solver_time(0.0),
        linear_solver_iterations(0) {}
    };

    template<class GOS, class TrlV, class TstV>
    class NewtonBase
    {
      typedef GOS GridOperator;
      typedef TrlV TrialVector;
      typedef TstV TestVector;

      typedef typename TestVector::ElementType RFType;
      typedef typename GOS::Traits::Jacobian Matrix;


    public:
      // export result type
      typedef NewtonResult<RFType> Result;

      void setVerbosityLevel(unsigned int verbosity_level)
      {
        if (gridoperator_.trialGridFunctionSpace().gridView().comm().rank()>0)
          verbosity_level_ = 0;
        else
          verbosity_level_ = verbosity_level;
      }

      //! Set whether the jacobian matrix should be kept across calls to apply().
      void setKeepMatrix(bool b)
      {
        keep_matrix_ = b;
      }

      //! Return whether the jacobian matrix is kept across calls to apply().
      bool keepMatrix() const
      {
        return keep_matrix_;
      }

      //! Discard the stored Jacobian matrix.
      void discardMatrix()
      {
        if(A_)
          A_.reset();
      }

    protected:
      const GridOperator& gridoperator_;
      TrialVector *u_;
      std::shared_ptr<TrialVector> z_;
      std::shared_ptr<TestVector> r_;
      std::shared_ptr<Matrix> A_;
      Result res_;
      unsigned int verbosity_level_;
      RFType prev_defect_;
      RFType linear_reduction_;
      bool reassembled_;
      RFType reduction_;
      RFType abs_limit_;
      bool keep_matrix_;

      NewtonBase(const GridOperator& go, TrialVector& u)
        : gridoperator_(go)
        , u_(&u)
        , verbosity_level_(1)
        , keep_matrix_(true)
      {
        if (gridoperator_.trialGridFunctionSpace().gridView().comm().rank()>0)
          verbosity_level_ = 0;
      }

      NewtonBase(const GridOperator& go)
        : gridoperator_(go)
        , u_(0)
        , verbosity_level_(1)
        , keep_matrix_(true)
      {
        if (gridoperator_.trialGridFunctionSpace().gridView().comm().rank()>0)
          verbosity_level_ = 0;
      }

      virtual ~NewtonBase() { }

      virtual bool terminate() = 0;
      virtual void prepare_step(Matrix& A, TestVector& r) = 0;
      virtual void line_search(TrialVector& z, TestVector& r) = 0;
      virtual void defect(TestVector& r) = 0;
    }; // end class NewtonBase

    template<class GOS, class S, class TrlV, class TstV>
    class NewtonSolver : public virtual NewtonBase<GOS,TrlV,TstV>
    {
      typedef S Solver;
      typedef GOS GridOperator;
      typedef TrlV TrialVector;
      typedef TstV TestVector;

      typedef typename TestVector::ElementType RFType;
      typedef typename GOS::Traits::Jacobian Matrix;

    public:
      typedef NewtonResult<RFType> Result;

      NewtonSolver(const GridOperator& go, TrialVector& u_, Solver& solver)
        : NewtonBase<GOS,TrlV,TstV>(go,u_)
        , solver_(solver)
        , result_valid_(false)
      {}

      NewtonSolver(const GridOperator& go, Solver& solver)
        : NewtonBase<GOS,TrlV,TstV>(go)
        , solver_(solver)
        , result_valid_(false)
      {}

      void apply();

      void apply(TrialVector& u_);

      const Result& result() const
      {
        if (!result_valid_)
          DUNE_THROW(NewtonError,
                     "NewtonSolver::result() called before NewtonSolver::solve()");
        return this->res_;
      }

    protected:
      virtual void defect(TestVector& r)
      {
        r = 0.0;
        this->gridoperator_.residual(*this->u_, r);
        this->res_.defect = this->solver_.norm(r);
        if (!std::isfinite(this->res_.defect))
          DUNE_THROW(NewtonDefectError,
                     "NewtonSolver::defect(): Non-linear defect is NaN or Inf");
      }


    private:
      void linearSolve(Matrix& A, TrialVector& z, TestVector& r) const
      {
        if (this->verbosity_level_ >= 4)
          std::cout << "      Solving linear system..." << std::endl;
        z = 0.0;
        // If possible tell solver backend to reuse linear system when we did not reassemble.
        Impl::setLinearSystemReuse(this->solver_, this->reassembled_);
        this->solver_.apply(A, z, r, this->linear_reduction_);

        ios_base_all_saver restorer(std::cout); // store old ios flags

        if (!this->solver_.result().converged)
          DUNE_THROW(NewtonLinearSolverError,
                     "NewtonSolver::linearSolve(): Linear solver did not converge "
                     "in " << this->solver_.result().iterations << " iterations");
        if (this->verbosity_level_ >= 4)
          std::cout << "          linear solver iterations:     "
                    << std::setw(12) << solver_.result().iterations << std::endl
                    << "          linear defect reduction:      "
                    << std::setw(12) << std::setprecision(4) << std::scientific
                    << solver_.result().reduction << std::endl;
      }

      Solver& solver_;
      bool result_valid_;
    }; // end class NewtonSolver

    template<class GOS, class S, class TrlV, class TstV>
    void NewtonSolver<GOS,S,TrlV,TstV>::apply(TrialVector& u)
    {
      this->u_ = &u;
      apply();
    }

    template<class GOS, class S, class TrlV, class TstV>
    void NewtonSolver<GOS,S,TrlV,TstV>::apply()
    {
      this->res_.iterations = 0;
      this->res_.converged = false;
      this->res_.reduction = 1.0;
      this->res_.conv_rate = 1.0;
      this->res_.elapsed = 0.0;
      this->res_.assembler_time = 0.0;
      this->res_.linear_solver_time = 0.0;
      this->res_.linear_solver_iterations = 0;
      result_valid_ = true;
      Timer timer;

      try
        {
          if(!this->r_) {
            // std::cout << "=== Setting up residual vector ..." << std::endl;
            this->r_ = std::make_shared<TestVector>(this->gridoperator_.testGridFunctionSpace());
          }
          // residual calculation in member function "defect":
          //--------------------------------------------------
          // - set residual vector to zero
          // - calculate new residual
          // - store norm of residual in "this->res_.defect"
          this->defect(*this->r_);
          this->res_.first_defect = this->res_.defect;
          this->prev_defect_ = this->res_.defect;

          if (this->verbosity_level_ >= 2)
            {
              // store old ios flags
              ios_base_all_saver restorer(std::cout);
              std::cout << "  Initial defect: "
                        << std::setw(12) << std::setprecision(4) << std::scientific
                        << this->res_.defect << std::endl;
            }

          if(!this->A_) {
            // std::cout << "==== Setting up jacobian matrix ... " << std::endl;
            this->A_ = std::make_shared<Matrix>(this->gridoperator_);
          }
          if(!this->z_) {
            // std::cout << "==== Setting up correction vector ... " << std::endl;
            this->z_ = std::make_shared<TrialVector>(this->gridoperator_.trialGridFunctionSpace());
          }

          while (!this->terminate())
            {
              if (this->verbosity_level_ >= 3)
                std::cout << "  Newton iteration " << this->res_.iterations
                          << " --------------------------------" << std::endl;

              Timer assembler_timer;
              try
                {
                  // jacobian calculation in member function "prepare_step"
                  //-------------------------------------------------------
                  // - if above reassemble threshold
                  //   - set jacobian to zero
                  //   - calculate new jacobian
                  // - set linear reduction
                  this->prepare_step(*this->A_,*this->r_);
                }
              catch (...)
                {
                  this->res_.assembler_time += assembler_timer.elapsed();
                  throw;
                }
              double assembler_time = assembler_timer.elapsed();
              this->res_.assembler_time += assembler_time;
              if (this->verbosity_level_ >= 3)
                std::cout << "      matrix assembly time:             "
                          << std::setw(12) << std::setprecision(4) << std::scientific
                          << assembler_time << std::endl;

              Timer linear_solver_timer;
              try
                {
                  // solution of linear system in member function "linearSolve"
                  //-----------------------------------------------------------
                  // - set initial guess for correction z to zero
                  // - call linear solver
                  this->linearSolve(*this->A_, *this->z_, *this->r_);
                }
              catch (...)
                {
                  this->res_.linear_solver_time += linear_solver_timer.elapsed();
                  this->res_.linear_solver_iterations += this->solver_.result().iterations;
                  throw;
                }
              double linear_solver_time = linear_solver_timer.elapsed();
              this->res_.linear_solver_time += linear_solver_time;
              this->res_.linear_solver_iterations += this->solver_.result().iterations;

              try
                {
                  // line search with correction z
                  // the undamped version is also integrated in here
                  this->line_search(*this->z_, *this->r_);
                }
              catch (NewtonLineSearchError)
                {
                  if (this->reassembled_)
                    throw;
                  if (this->verbosity_level_ >= 3)
                    std::cout << "      line search failed - trying again with reassembled matrix" << std::endl;
                  continue;
                }

              this->res_.reduction = this->res_.defect/this->res_.first_defect;
              this->res_.iterations++;
              this->res_.conv_rate = std::pow(this->res_.reduction, 1.0/this->res_.iterations);

              // store old ios flags
              ios_base_all_saver restorer(std::cout);

              if (this->verbosity_level_ >= 3)
                std::cout << "      linear solver time:               "
                          << std::setw(12) << std::setprecision(4) << std::scientific
                          << linear_solver_time << std::endl
                          << "      defect reduction (this iteration):"
                          << std::setw(12) << std::setprecision(4) << std::scientific
                          << this->res_.defect/this->prev_defect_ << std::endl
                          << "      defect reduction (total):         "
                          << std::setw(12) << std::setprecision(4) << std::scientific
                          << this->res_.reduction << std::endl
                          << "      new defect:                       "
                          << std::setw(12) << std::setprecision(4) << std::scientific
                          << this->res_.defect << std::endl;
              if (this->verbosity_level_ == 2)
                std::cout << "  Newton iteration " << std::setw(2) << this->res_.iterations
                          << ".  New defect: "
                          << std::setw(12) << std::setprecision(4) << std::scientific
                          << this->res_.defect
                          << ".  Reduction (this): "
                          << std::setw(12) << std::setprecision(4) << std::scientific
                          << this->res_.defect/this->prev_defect_
                          << ".  Reduction (total): "
                          << std::setw(12) << std::setprecision(4) << std::scientific
                          << this->res_.reduction << std::endl;
            } // end while
        } // end try
      catch(...)
        {
          this->res_.elapsed = timer.elapsed();
          throw;
        }
      this->res_.elapsed = timer.elapsed();

      ios_base_all_saver restorer(std::cout); // store old ios flags

      if (this->verbosity_level_ == 1)
        std::cout << "  Newton converged after " << std::setw(2) << this->res_.iterations
                  << " iterations.  Reduction: "
                  << std::setw(12) << std::setprecision(4) << std::scientific
                  << this->res_.reduction
                  << "   (" << std::setprecision(4) << this->res_.elapsed << "s)"
                  << std::endl;

      if(!this->keep_matrix_)
        this->A_.reset();
    } // end apply

    template<class GOS, class TrlV, class TstV>
    class NewtonTerminate : public virtual NewtonBase<GOS,TrlV,TstV>
    {
      typedef GOS GridOperator;
      typedef TrlV TrialVector;

      typedef typename TstV::ElementType RFType;

    public:
      NewtonTerminate(const GridOperator& go, TrialVector& u_)
        : NewtonBase<GOS,TrlV,TstV>(go,u_)
        , maxit_(40)
        , force_iteration_(false)
      {
        this->reduction_ = 1e-8;
        this->abs_limit_ = 1e-12;
      }

      NewtonTerminate(const GridOperator& go)
        : NewtonBase<GOS,TrlV,TstV>(go)
        , maxit_(40)
        , force_iteration_(false)
      {
        this->reduction_ = 1e-8;
        this->abs_limit_ = 1e-12;
      }

      void setReduction(RFType reduction)
      {
        this->reduction_ = reduction;
      }

      void setMaxIterations(unsigned int maxit)
      {
        maxit_ = maxit;
      }

      void setForceIteration(bool force_iteration)
      {
        force_iteration_ = force_iteration;
      }

      void setAbsoluteLimit(RFType abs_limit_)
      {
        this->abs_limit_ = abs_limit_;
      }

      virtual bool terminate()
      {
        if (force_iteration_ && this->res_.iterations == 0)
          return false;
        this->res_.converged = this->res_.defect < this->abs_limit_
          || this->res_.defect < this->res_.first_defect * this->reduction_;
        if (this->res_.iterations >= maxit_ && !this->res_.converged)
          DUNE_THROW(NewtonNotConverged,
                     "NewtonTerminate::terminate(): Maximum iteration count reached");
        return this->res_.converged;
      }

    private:
      unsigned int maxit_;
      bool force_iteration_;
    }; // end class NewtonTerminate

    template<class GOS, class TrlV, class TstV>
    class NewtonPrepareStep : public virtual NewtonBase<GOS,TrlV,TstV>
    {
      typedef GOS GridOperator;
      typedef TrlV TrialVector;

      typedef typename TstV::ElementType RFType;
      typedef typename GOS::Traits::Jacobian Matrix;

    public:
      NewtonPrepareStep(const GridOperator& go, TrialVector& u_)
        : NewtonBase<GOS,TrlV,TstV>(go,u_)
        , min_linear_reduction_(1e-3)
        , fixed_linear_reduction_(0.0)
        , reassemble_threshold_(0.0)
      {}

      NewtonPrepareStep(const GridOperator& go)
        : NewtonBase<GOS,TrlV,TstV>(go)
        , min_linear_reduction_(1e-3)
        , fixed_linear_reduction_(0.0)
        , reassemble_threshold_(0.0)
      {}

      /**\brief set the minimal reduction in the linear solver

         \note with min_linear_reduction > 0, the linear reduction will be
         determined as mininum of the min_linear_reduction and the
         linear_reduction needed to achieve second order
         Newton convergence. */
      void setMinLinearReduction(RFType min_linear_reduction)
      {
        min_linear_reduction_ = min_linear_reduction;
      }

      /**\brief set a fixed reduction in the linear solver (overwrites setMinLinearReduction)

         \note with fixed_linear_reduction > 0, the linear reduction
         rate will always be fixed to min_linear_reduction. */
      void setFixedLinearReduction(bool fixed_linear_reduction)
      {
        fixed_linear_reduction_ = fixed_linear_reduction;
      }

      /**\brief set a threshold, when the linear operator is reassembled

         We allow to keep the linear operator over several newton
         iterations. If the reduction in the newton drops below a
         given threshold the linear operator is reassembled to ensure
         convergence.
       */
      void setReassembleThreshold(RFType reassemble_threshold)
      {
        reassemble_threshold_ = reassemble_threshold;
      }

      virtual void prepare_step(Matrix& A, TstV& )
      {
        this->reassembled_ = false;
        if (this->res_.defect/this->prev_defect_ > reassemble_threshold_)
          {
            if (this->verbosity_level_ >= 3)
              std::cout << "      Reassembling matrix..." << std::endl;
            A = 0.0;
            this->gridoperator_.jacobian(*this->u_, A);
            this->reassembled_ = true;
          }

        if (fixed_linear_reduction_ == true)
          this->linear_reduction_ = min_linear_reduction_;
        else {
          // determine maximum defect, where Newton is converged.
          RFType stop_defect =
            std::max(this->res_.first_defect * this->reduction_,
                     this->abs_limit_);

          /*
            To achieve second order convergence of newton
            we need a linear reduction of at least
            current_defect^2/prev_defect^2.
            For the last newton step a linear reduction of
            1/10*end_defect/current_defect
            is sufficient for convergence.
          */
          if ( stop_defect/(10*this->res_.defect) >
               this->res_.defect*this->res_.defect/(this->prev_defect_*this->prev_defect_) )
            this->linear_reduction_ =
              stop_defect/(10*this->res_.defect);
          else
            this->linear_reduction_ =
              std::min(min_linear_reduction_,this->res_.defect*this->res_.defect/(this->prev_defect_*this->prev_defect_));
        }

        this->prev_defect_ = this->res_.defect;

        ios_base_all_saver restorer(std::cout); // store old ios flags

        if (this->verbosity_level_ >= 3)
          std::cout << "      requested linear reduction:       "
                    << std::setw(12) << std::setprecision(4) << std::scientific
                    << this->linear_reduction_ << std::endl;
      }

    private:
      RFType min_linear_reduction_;
      bool fixed_linear_reduction_;
      RFType reassemble_threshold_;
    }; // end class NewtonPrepareStep

    template<class GOS, class TrlV, class TstV>
    class NewtonLineSearch : public virtual NewtonBase<GOS,TrlV,TstV>
    {
      typedef GOS GridOperator;
      typedef TrlV TrialVector;
      typedef TstV TestVector;

      typedef typename TestVector::ElementType RFType;

    public:
      enum Strategy {
        /** \brief don't do any line search or damping */
        noLineSearch,
        /** \brief perform a linear search for the optimal damping parameter with multiples of damping

         the strategy was described in <a href="http://dx.doi.org/10.1007/BF01406516">[Hackbusch and Reusken, 1989]</a> */
        hackbuschReusken,
        /** \brief same as hackbuschReusken, but doesn't fail if the best update is still not good enough */
        hackbuschReuskenAcceptBest };

      NewtonLineSearch(const GridOperator& go, TrialVector& u_)
        : NewtonBase<GOS,TrlV,TstV>(go,u_)
        , strategy_(hackbuschReusken)
        , maxit_(10)
        , damping_factor_(0.5)
      {}

      NewtonLineSearch(const GridOperator& go)
        : NewtonBase<GOS,TrlV,TstV>(go)
        , strategy_(hackbuschReusken)
        , maxit_(10)
        , damping_factor_(0.5)
      {}

      void setLineSearchStrategy(Strategy strategy)
      {
        strategy_ = strategy;
      }

      void setLineSearchStrategy(std::string strategy)
      {
        strategy_ = strategyFromName(strategy);
      }

      void setLineSearchMaxIterations(unsigned int maxit)
      {
        maxit_ = maxit;
      }

      void setLineSearchDampingFactor(RFType damping_factor)
      {
        damping_factor_ = damping_factor;
      }

      virtual void line_search(TrialVector& z, TestVector& r)
      {
        if (strategy_ == noLineSearch)
          {
            this->u_->axpy(-1.0, z);
            this->defect(r);
            return;
          }

        if (this->verbosity_level_ >= 4)
          std::cout << "      Performing line search..." << std::endl;
        RFType lambda = 1.0;
        RFType best_lambda = 0.0;
        RFType best_defect = this->res_.defect;
        TrialVector prev_u(*this->u_);
        unsigned int i = 0;
        ios_base_all_saver restorer(std::cout); // store old ios flags

        while (1)
          {
            if (this->verbosity_level_ >= 4)
              std::cout << "          trying line search damping factor:   "
                        << std::setw(12) << std::setprecision(4) << std::scientific
                        << lambda
                        << std::endl;

            this->u_->axpy(-lambda, z);
            try {
              this->defect(r);
            }
             catch (NewtonDefectError)
              {
                if (this->verbosity_level_ >= 4)
                  std::cout << "          NaNs detected" << std::endl;
              }       // ignore NaNs and try again with lower lambda

            if (this->res_.defect <= (1.0 - lambda/4) * this->prev_defect_)
              {
                if (this->verbosity_level_ >= 4)
                  std::cout << "          line search converged" << std::endl;
                break;
              }

            if (this->res_.defect < best_defect)
              {
                best_defect = this->res_.defect;
                best_lambda = lambda;
              }

            if (++i >= maxit_)
              {
                if (this->verbosity_level_ >= 4)
                  std::cout << "          max line search iterations exceeded" << std::endl;
                switch (strategy_)
                  {
                  case hackbuschReusken:
                    *this->u_ = prev_u;
                    this->defect(r);
                    DUNE_THROW(NewtonLineSearchError,
                               "NewtonLineSearch::line_search(): line search failed, "
                               "max iteration count reached, "
                               "defect did not improve enough");
                  case hackbuschReuskenAcceptBest:
                    if (best_lambda == 0.0)
                      {
                        *this->u_ = prev_u;
                        this->defect(r);
                        DUNE_THROW(NewtonLineSearchError,
                                   "NewtonLineSearch::line_search(): line search failed, "
                                   "max iteration count reached, "
                                   "defect did not improve in any of the iterations");
                      }
                    if (best_lambda != lambda)
                      {
                        *this->u_ = prev_u;
                        this->u_->axpy(-best_lambda, z);
                        this->defect(r);
                      }
                    break;
                  case noLineSearch:
                    break;
                  }
                break;
              }

            lambda *= damping_factor_;
            *this->u_ = prev_u;
          }
        if (this->verbosity_level_ >= 4)
          std::cout << "          line search damping factor:   "
                    << std::setw(12) << std::setprecision(4) << std::scientific
                    << lambda << std::endl;
      } // end line_search

    protected:
      /** helper function to get the different strategies from their name */
      Strategy strategyFromName(const std::string & s) {
        if (s == "noLineSearch")
          return noLineSearch;
        if (s == "hackbuschReusken")
          return hackbuschReusken;
        if (s == "hackbuschReuskenAcceptBest")
          return hackbuschReuskenAcceptBest;
        DUNE_THROW(Exception, "unknown line search strategy" << s);
      }

    private:
      Strategy strategy_;
      unsigned int maxit_;
      RFType damping_factor_;
    }; // end class NewtonLineSearch

    template<class GOS, class S, class TrlV, class TstV = TrlV>
    class Newton : public NewtonSolver<GOS,S,TrlV,TstV>
                 , public NewtonTerminate<GOS,TrlV,TstV>
                 , public NewtonLineSearch<GOS,TrlV,TstV>
                 , public NewtonPrepareStep<GOS,TrlV,TstV>
    {
      typedef GOS GridOperator;
      typedef S Solver;
      typedef TrlV TrialVector;

    public:
      Newton(const GridOperator& go, TrialVector& u_, Solver& solver_)
        : NewtonBase<GOS,TrlV,TstV>(go,u_)
        , NewtonSolver<GOS,S,TrlV,TstV>(go,u_,solver_)
        , NewtonTerminate<GOS,TrlV,TstV>(go,u_)
        , NewtonLineSearch<GOS,TrlV,TstV>(go,u_)
        , NewtonPrepareStep<GOS,TrlV,TstV>(go,u_)
      {}
      Newton(const GridOperator& go, Solver& solver_)
        : NewtonBase<GOS,TrlV,TstV>(go)
        , NewtonSolver<GOS,S,TrlV,TstV>(go,solver_)
        , NewtonTerminate<GOS,TrlV,TstV>(go)
        , NewtonLineSearch<GOS,TrlV,TstV>(go)
        , NewtonPrepareStep<GOS,TrlV,TstV>(go)
      {}
      //! interpret a parameter tree as a set of options for the newton solver
      /**

         example configuration:

         \code
         [NewtonParameters]

         ReassembleThreshold = 0.1
         LineSearchMaxIterations = 10
         MaxIterations = 7
         AbsoluteLimit = 1e-6
         Reduction = 1e-4
         MinLinearReduction = 1e-3
         LineSearchDamping  = 0.9
         \endcode

         and invocation in the code:
         \code
         newton.setParameters(param.sub("NewtonParameters"));
         \endcode
      */
      void setParameters(const Dune::ParameterTree & param)
      {
        typedef typename TstV::ElementType RFType;
        if (param.hasKey("VerbosityLevel"))
          this->setVerbosityLevel(
            param.get<unsigned int>("VerbosityLevel"));
        if (param.hasKey("Reduction"))
          this->setReduction(
            param.get<RFType>("Reduction"));
        if (param.hasKey("MaxIterations"))
          this->setMaxIterations(
            param.get<unsigned int>("MaxIterations"));
        if (param.hasKey("ForceIteration"))
          this->setForceIteration(
            param.get<bool>("ForceIteration"));
        if (param.hasKey("AbsoluteLimit"))
          this->setAbsoluteLimit(
            param.get<RFType>("AbsoluteLimit"));
        if (param.hasKey("MinLinearReduction"))
          this->setMinLinearReduction(
            param.get<RFType>("MinLinearReduction"));
        if (param.hasKey("FixedLinearReduction"))
          this->setFixedLinearReduction(
            param.get<bool>("FixedLinearReduction"));
        if (param.hasKey("ReassembleThreshold"))
          this->setReassembleThreshold(
            param.get<RFType>("ReassembleThreshold"));
        if (param.hasKey("LineSearchStrategy"))
          this->setLineSearchStrategy(
            param.get<std::string>("LineSearchStrategy"));
        if (param.hasKey("LineSearchMaxIterations"))
          this->setLineSearchMaxIterations(
            param.get<unsigned int>("LineSearchMaxIterations"));
        if (param.hasKey("LineSearchDampingFactor"))
          this->setLineSearchDampingFactor(
            param.get<RFType>("LineSearchDampingFactor"));
        if (param.hasKey("KeepMatrix"))
          this->setKeepMatrix(
            param.get<bool>("KeepMatrix"));
      }
    }; // end class Newton
  } // end namespace PDELab
} // end namespace Dune

#endif // DUNE_PDELAB_NEWTON_NEWTON_HH
