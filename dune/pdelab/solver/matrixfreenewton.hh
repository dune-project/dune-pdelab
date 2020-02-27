// -*- tab-width: 2; indent-tabs-mode: nil -*-
// vi: set et ts=2 sw=2 sts=2:

#ifndef DUNE_PDELAB_SOLVER_MATRIXFREENEWTON_HH
#define DUNE_PDELAB_SOLVER_MATRIXFREENEWTON_HH

//
// Note: This is basically a copy of the old Newton solver. A smarter
// implementation that avoids code duplication needs a change in the solver
// backends.
//

#include <iostream>
#include <iomanip>
#include <cmath>
#include <memory>

#include <math.h>

#include <dune/common/ios_state.hh>
#include <dune/common/timer.hh>
#include <dune/common/parametertree.hh>

#include <dune/pdelab/newton/newton.hh>
#include <dune/pdelab/solver/newtonerrors.hh>
#include <dune/pdelab/solver/newton.hh>

namespace Dune {
  namespace PDELab {

    enum struct MatrixFreeLineSearchStrategy {

      /** \brief don't do any linesearch or damping */
      noLineSearch,

      /** \brief perform a linear search for the optimal damping parameter with multiples of damping

          the strategy was described in <a href="http://dx.doi.org/10.1007/BF01406516">[Hackbusch and Reusken, 1989]</a> */
     hackbuschReusken,

      /** \brief same as hackbuschReusken, but doesn't fail if the best update is still not good enough */
      hackbuschReuskenAcceptBest
      };


    template<typename GO, typename S, typename TrlV, typename TstV = TrlV>
    class MatrixFreeNewton
    {
      typedef GO GridOperator;
      typedef S Solver;
      typedef TrlV TrialVector;
      typedef TstV TestVector;

      typedef typename TestVector::ElementType RFType;

    public :
      //! export result type
      using Result = NewtonResult<RFType>;

      MatrixFreeNewton(const GridOperator& go, TrialVector& u, Solver& solver)
        : gridoperator_(go)
        , u_(&u)
        , solver_(solver)
        , verbosity_level_(1), maxit_(40), force_iteration_(false)
        , min_linear_reduction_(1e-3), fixed_linear_reduction_(false)
        , strategy_(MatrixFreeLineSearchStrategy::hackbuschReusken), linesearch_maxit_(10), damping_factor_(0.5)
        , reduction_(1e-8), abs_limit_(1e-12), result_valid_(false)
      {
        if (gridoperator_.trialGridFunctionSpace().gridView().comm().rank()>0)
          verbosity_level_ = 0;
      }
      MatrixFreeNewton(const GridOperator& go, Solver& solver)
        : gridoperator_(go)
        , u_(static_cast<TrialVector*>(0))
        , solver_(solver)
        , verbosity_level_(1), maxit_(40), force_iteration_(false)
        , min_linear_reduction_(1e-3), fixed_linear_reduction_(false)
        , strategy_(MatrixFreeLineSearchStrategy::hackbuschReusken), linesearch_maxit_(10), damping_factor_(0.5)
        , reduction_(1e-8), abs_limit_(1e-12), result_valid_(false)
      {
        if (gridoperator_.trialGridFunctionSpace().gridView().comm().rank()>0)
          verbosity_level_ = 0;
      }
      //! interpret a parameter tree as a set of options for the newton solver
      /**

         example configuration:

         \code
         [NewtonParameters]

         LineSearchMaxIterations = 10
         MaxIterations = 7
         AbsoluteLimit = 1e-6
         Reduction = 1e-4
         LinearReduction = 1e-3
         LineSearchDamping  = 0.9
         \endcode

         and invocation in the code:
         \code
         newton.setParameters(param.sub("NewtonParameters"));
         \endcode
      */
      void setParameters(Dune::ParameterTree & param)
      {
        typedef typename TstV::ElementType RFType;
        if (param.hasKey("VerbosityLevel"))
          setVerbosityLevel(
            param.get<unsigned int>("VerbosityLevel"));
        if (param.hasKey("Reduction"))
          setReduction(
            param.get<RFType>("Reduction"));
        if (param.hasKey("MaxIterations"))
          setMaxIterations(
            param.get<unsigned int>("MaxIterations"));
        if (param.hasKey("ForceIteration"))
          setForceIteration(
            param.get<bool>("ForceIteration"));
        if (param.hasKey("AbsoluteLimit"))
          setAbsoluteLimit(
            param.get<RFType>("AbsoluteLimit"));
        if (param.hasKey("MinLinearReduction"))
          setMinLinearReduction(
            param.get<RFType>("MinLinearReduction"));
        if (param.hasKey("FixedLinearReduction"))
          setFixedLinearReduction(
            param.get<bool>("FixedLinearReduction"));
        if (param.hasKey("LineSearchStrategy"))
          setLineSearchStrategy(
            param.get<std::string>("LineSearchStrategy"));
        if (param.hasKey("LineSearchMaxIterations"))
          setLineSearchMaxIterations(
            param.get<unsigned int>("LineSearchMaxIterations"));
        if (param.hasKey("LineSearchDampingFactor"))
          setLineSearchDampingFactor(
            param.get<RFType>("LineSearchDampingFactor"));
      }

      //! set verbosity level
      void setVerbosityLevel(unsigned int verbosity_level)
      {
        if (gridoperator_.trialGridFunctionSpace().gridView().comm().rank()>0)
          verbosity_level_ = 0;
        else
          verbosity_level_ = verbosity_level;
      }

      //! set reduction in Newton
      void setReduction(RFType reduction)
      {
        reduction_ = reduction;
      }

      //! set maximum of nonlinear iterations
      void setMaxIterations(unsigned int maxit)
      {
        maxit_ = maxit;
      }

      void setForceIteration(bool force_iteration)
      {
        force_iteration_ = force_iteration;
      }

      //! set absolute limit to achieve for the norm of the residual
      void setAbsoluteLimit(RFType abs_limit)
      {
        abs_limit_ = abs_limit;
      }

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

      //! set line search strategy
      void setLineSearchStrategy(LineSearchStrategy strategy)
      {
        strategy_ = strategy;
      }

      void setLineSearchStrategy(std::string strategy)
      {
        // read strategy from string
        if(strategy == "noLineSearch") {
          strategy_ = MatrixFreeLineSearchStrategy::noLineSearch;
          return;
        }
        if(strategy == "hackbuschReusken") {
          strategy_ = MatrixFreeLineSearchStrategy::hackbuschReusken;
          return;
        }
        if(strategy == "hackbuschReuskenAcceptBest") {
          strategy_ = MatrixFreeLineSearchStrategy::hackbuschReuskenAcceptBest;
          return;
        }
        DUNE_THROW(Exception, "unknown linesearch strategy" << strategy);
      }

      //! set maximum number of line search iterations
      void setLineSearchMaxIterations(unsigned int linesearch_maxit)
      {
        linesearch_maxit_ = linesearch_maxit;
      }

      //! set damping factor for line search parameter
      void setLineSearchDampingFactor(RFType damping_factor)
      {
        damping_factor_ = damping_factor;
      }

      //! stopping criteria for Newton
      bool terminate()
      {
        if(force_iteration_ and res_.iterations == 0)
          return false;

        // standard criterion
        res_.converged = (res_.defect < abs_limit_) or (res_.defect < res_.first_defect * reduction_);

        if(res_.iterations >= maxit_ and (not res_.converged))
          DUNE_THROW(NewtonNotConverged,
                     "NewtonTerminate::terminate(): Maximum iteration count reached");
        return res_.converged;
      }

      //! line search in Newton
      void lineSearch(TrialVector& z, TrialVector& r)
      {
        if (strategy_ == MatrixFreeLineSearchStrategy::noLineSearch)
          {
            u_->axpy(-1.0, z);
            r = 0.0;
            gridoperator_.residual(*u_,r);
            res_.defect = solver_.norm(r);
            if(not std::isfinite(res_.defect))
              DUNE_THROW(NewtonDefectError,
                         "Non linear residual is NaN or Inf!");
            return;
          }
        if (verbosity_level_ >= 4)
          std::cout << "      Performing line search..." << std::endl;
        RFType lambda = 1.0;
        RFType best_lambda = 0.0;
        RFType best_defect = res_.defect;
        TrialVector prev_u(*u_);
        unsigned int i = 0;
        ios_base_all_saver restorer(std::cout); // store old ios flags

        while (1)
          {
            if (verbosity_level_ >= 4)
              std::cout << "          trying line search damping factor:   "
                        << std::setw(12) << std::setprecision(4) << std::scientific
                        << lambda
                        << std::endl;

            u_->axpy(-lambda, z);
            try {
              r = 0.0;
              gridoperator_.residual(*u_,r);
              res_.defect = solver_.norm(r);
              if(not std::isfinite(res_.defect))
                DUNE_THROW(NewtonDefectError,
                           "Non linear residual is NaN or Inf!");
            }
             catch (NewtonDefectError)
              {
                if (verbosity_level_ >= 4)
                  std::cout << "          Nans detected" << std::endl;
              }       // ignore NaNs and try again with lower lambda

            if (res_.defect <= (1.0 - lambda/4) * prev_defect_)
              {
                if (verbosity_level_ >= 4)
                  std::cout << "          line search converged" << std::endl;
                break;
              }

            if (res_.defect < best_defect)
              {
                best_defect = res_.defect;
                best_lambda = lambda;
              }

            if (++i >= maxit_)
              {
                if (verbosity_level_ >= 4)
                  std::cout << "          max line search iterations exceeded" << std::endl;
                switch (strategy_)
                  {
                  case MatrixFreeLineSearchStrategy::hackbuschReusken:
                    *u_ = prev_u;
                    r = 0.0;
                    gridoperator_.residual(*u_,r);
                    res_.defect = solver_.norm(r);
                    if(not std::isfinite(res_.defect))
                      DUNE_THROW(NewtonDefectError,
                                 "Non linear residual is NaN or Inf!");
                    DUNE_THROW(NewtonLineSearchError,
                               "NewtonLineSearch::line_search(): line search failed, "
                               "max iteration count reached, "
                               "defect did not improve enough");
                  case MatrixFreeLineSearchStrategy::hackbuschReuskenAcceptBest:
                    if (best_lambda == 0.0)
                      {
                        *u_ = prev_u;
                        r = 0.0;
                        gridoperator_.residual(*u_,r);
                        res_.defect = solver_.norm(r);
                        if(not std::isfinite(res_.defect))
                          DUNE_THROW(NewtonDefectError,
                                     "Non linear residual is NaN or Inf!");
                        DUNE_THROW(NewtonLineSearchError,
                                   "NewtonLineSearch::line_search(): line search failed, "
                                   "max iteration count reached, "
                                   "defect did not improve in any of the iterations");
                      }
                    if (best_lambda != lambda)
                      {
                        *u_ = prev_u;
                        u_->axpy(-best_lambda, z);
                        r = 0.0;
                        gridoperator_.residual(*u_,r);
                        res_.defect = solver_.norm(r);
                        if(not std::isfinite(res_.defect))
                          DUNE_THROW(NewtonDefectError,
                                     "Non linear residual is NaN or Inf!");
                      }
                    break;
                  case MatrixFreeLineSearchStrategy::noLineSearch:
                    break;
                  }
                break;
              }

            lambda *= damping_factor_;
            *u_ = prev_u;
          }
        if (verbosity_level_ >= 4)
          std::cout << "          line search damping factor:   "
                    << std::setw(12) << std::setprecision(4) << std::scientific
                    << lambda << std::endl;
      } // end lineSearch

      const Result& result() const
      {
        return res_;
      }

      void apply(TrialVector& u)
      {
        u_ = &u;
        apply();
      }

      void apply()
      {
        // reset solver statistics
        res_.iterations = 0;
        res_.converged = false;
        res_.reduction = 1.0;
        res_.conv_rate = 1.0;
        res_.elapsed = 0.0;
        res_.linear_solver_time = 0.0;
        res_.linear_solver_iterations = 0;
        result_valid_ = true;
        Timer timer;

        // allocate only once the residual vector
        if(not r_)
          r_ = std::make_shared<TestVector>(gridoperator_.testGridFunctionSpace());

        //============================================
        // calculate initial residual and its norm
        //============================================
        *r_ = 0.0;
        gridoperator_.residual(*u_,*r_);
        res_.defect = solver_.norm(*r_);
        if(not std::isfinite(res_.defect))
          DUNE_THROW(NewtonDefectError,
                     "Non linear residual is NaN or Inf!");
        res_.first_defect = res_.defect;
        prev_defect_ = res_.defect;

        if(verbosity_level_ >= 2) {
          // store old ios flags
          ios_base_all_saver restorer(std::cout);
          std::cout << "  Initial defect: "
                    << std::setw(12) << std::setprecision(4) << std::scientific
                    << res_.defect << std::endl;
        }

        // allocate the correction vector only one
        if(not z_)
          z_ = std::make_shared<TrialVector>(gridoperator_.trialGridFunctionSpace());

        //============================================
        // the nonlinear iteration loop
        //============================================
        while(not terminate()) {
          if(verbosity_level_ >= 3)
            std::cout << "  Newton iteration " << res_.iterations
                      << " --------------------------------" << std::endl;

          //============================================
          // solution of linear system
          //============================================
          Timer linear_solver_timer;
          if(verbosity_level_ >= 4)
            std::cout << "      Solving linear system..." << std::endl;
          // <<<1>>> set linear reduction
          if(fixed_linear_reduction_ == true)
            linear_reduction_ = min_linear_reduction_;
          else {
            // determine maximum defect, where Newton is converged.
            auto stop_defect =
              std::max(res_.first_defect * reduction_, abs_limit_);

            /*
              To achieve second order convergence of newton
              we need a linear reduction of at least
              current_defect^2/prev_defect^2.
              For the last newton step a linear reduction of
              1/10*end_defect/current_defect
              is sufficient for convergence.
            */
            if ( stop_defect/(10*res_.defect) >
                 res_.defect*res_.defect/(prev_defect_*prev_defect_) )
              linear_reduction_ =
                stop_defect/(10*res_.defect);
            else
              linear_reduction_ =
                std::min(min_linear_reduction_,res_.defect*res_.defect/(prev_defect_*prev_defect_));
          }

          ios_base_all_saver restorer(std::cout); // store old ios flags

          if (verbosity_level_ >= 3)
            std::cout << "      requested linear reduction:       "
                      << std::setw(12) << std::setprecision(4) << std::scientific
                      << linear_reduction_ << std::endl;

          // <<<2>>> set position of jacobian for its application
          solver_.setLinearizationPoint(*u_);

          // <<<3>>> call linear solver
          *z_ = 0.0;
          solver_.apply(*z_, *r_, linear_reduction_);
          res_.linear_solver_time += linear_solver_timer.elapsed();
          res_.linear_solver_iterations += solver_.result().iterations;
          auto linear_solver_time = linear_solver_timer.elapsed();
          res_.linear_solver_time += linear_solver_time;
          res_.linear_solver_iterations += solver_.result().iterations;

          //============================================
          // line search with correction z
          // the undamped version is also covered in here
          //============================================
          lineSearch(*z_,*r_);

          // store statistics and output per nonlinear iteration
          res_.reduction = res_.defect/res_.first_defect;
          res_.iterations++;
          res_.conv_rate = std::pow(res_.reduction, 1.0/res_.iterations);

          if (verbosity_level_ >= 3)
            std::cout << "      linear solver time:               "
                      << std::setw(12) << std::setprecision(4) << std::scientific
                      << linear_solver_time << std::endl
                      << "      defect reduction (this iteration):"
                      << std::setw(12) << std::setprecision(4) << std::scientific
                      << res_.defect/prev_defect_ << std::endl
                      << "      defect reduction (total):         "
                      << std::setw(12) << std::setprecision(4) << std::scientific
                      << res_.reduction << std::endl
                      << "      new defect:                       "
                      << std::setw(12) << std::setprecision(4) << std::scientific
                      << res_.defect << std::endl;
          if (verbosity_level_ == 2)
            std::cout << "  Newton iteration " << std::setw(2) << res_.iterations
                      << ".  New defect: "
                      << std::setw(12) << std::setprecision(4) << std::scientific
                      << res_.defect
                      << ".  Reduction (this): "
                      << std::setw(12) << std::setprecision(4) << std::scientific
                      << res_.defect/prev_defect_
                      << ".  Reduction (total): "
                      << std::setw(12) << std::setprecision(4) << std::scientific
                      << res_.reduction << std::endl;
        } // end while

        res_.elapsed = timer.elapsed();

        ios_base_all_saver restorer(std::cout); // store old ios flags

        if (verbosity_level_ == 1)
          std::cout << "  Newton converged after " << std::setw(2) << res_.iterations
                    << " iterations.  Reduction: "
                    << std::setw(12) << std::setprecision(4) << std::scientific
                    << res_.reduction
                    << "   (" << std::setprecision(4) << res_.elapsed << "s)"
                    << std::endl;
      } // end apply

    private :
      const GridOperator& gridoperator_;
      TrialVector *u_;
      Solver& solver_;
      std::shared_ptr<TrialVector> z_;
      std::shared_ptr<TestVector> r_;
      Result res_;
      unsigned int verbosity_level_;
      unsigned int maxit_;
      bool force_iteration_;
      RFType min_linear_reduction_;
      bool fixed_linear_reduction_;
      MatrixFreeLineSearchStrategy strategy_;
      unsigned int linesearch_maxit_;
      RFType damping_factor_;
      RFType prev_defect_;
      RFType linear_reduction_;
      RFType reduction_;
      RFType abs_limit_;
      bool result_valid_;
    }; // end class MatrixFreeNewton

  } // end namespace PDELab
} // end namespace Dune

#endif
