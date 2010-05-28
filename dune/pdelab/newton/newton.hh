// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=8 sw=4 sts=4:
#ifndef DUNE_PDELAB_NEWTON_HH
#define DUNE_PDELAB_NEWTON_HH

#include <iostream>
#include <iomanip>
#include <cmath>

#include <math.h>

#include <dune/common/exceptions.hh>
#include <dune/common/timer.hh>
#include <dune/common/ios.hh>

namespace Dune
{
    namespace PDELab
    {
        // Exception classes used in NewtonSolver
        class NewtonError : public Exception {};
        class NewtonDefectError : public NewtonError {};
        class NewtonLinearSolverError : public NewtonError {};
        class NewtonLineSearchError : public NewtonError {};
        class NewtonNotConverged : public NewtonError {};

        // Status information of a linear solver
        template<class RFType>
        struct LinearSolverResult
        {
            bool converged;            // Solver converged
            unsigned int iterations;   // number of iterations
            double elapsed;            // total user time in seconds
            RFType reduction;          // defect reduction
        };

        // Status information of Newton's method
        template<class RFType>
        struct NewtonResult : LinearSolverResult<RFType>
        {
            RFType conv_rate;          // average reduction per Newton iteration
            RFType first_defect;       // the first defect
            RFType defect;             // the final defect
        };

        template<class GOS, class TrlV, class TstV>
        class NewtonBase
        {
            typedef GOS GridOperator;
            typedef TrlV TrialVector;
            typedef TstV TestVector;

            typedef typename TestVector::ElementType RFType;
            typedef typename GOS::template MatrixContainer<RFType>::Type Matrix;

            typedef NewtonResult<RFType> Result;

        public:
            void setVerbosityLevel(unsigned int verbosity_level_)
            {
                if (gridoperator.trialGridFunctionSpace().gridview().comm().rank()>0)
                    verbosity_level = 0;
                else
                    verbosity_level = verbosity_level_;
            }

        protected:
            GridOperator& gridoperator;
            TrialVector *u;
            Result res;
            unsigned int verbosity_level;
            RFType prev_defect;
            RFType linear_reduction;
            bool reassembled;

            NewtonBase(GridOperator& go, TrialVector& u_)
                : gridoperator(go)
                , u(&u_)
                , verbosity_level(1)
            {
                if (gridoperator.trialGridFunctionSpace().gridview().comm().rank()>0)
                    verbosity_level = 0;
            }

            NewtonBase(GridOperator& go)
                : gridoperator(go)
                , u(0)
                , verbosity_level(1)
            {
                if (gridoperator.trialGridFunctionSpace().gridview().comm().rank()>0)
                    verbosity_level = 0;
            }

            virtual bool terminate() = 0;
            virtual void prepare_step(Matrix& A) = 0;
            virtual void line_search(TrialVector& z, TestVector& r) = 0;
            virtual void defect(TestVector& r) = 0;
        };

        template<class GOS, class S, class TrlV, class TstV>
        class NewtonSolver : public virtual NewtonBase<GOS,TrlV,TstV>
        {
            typedef S Solver;
            typedef GOS GridOperator;
            typedef TrlV TrialVector;
            typedef TstV TestVector;

            typedef typename TestVector::ElementType RFType;
            typedef typename GOS::template MatrixContainer<RFType>::Type Matrix;

        public:
            typedef NewtonResult<RFType> Result;

            NewtonSolver(GridOperator& go, TrialVector& u_, Solver& solver_)
                : NewtonBase<GOS,TrlV,TstV>(go,u_)
                , solver(solver_)
                , result_valid(false)
            {}

            NewtonSolver(GridOperator& go, Solver& solver_)
                : NewtonBase<GOS,TrlV,TstV>(go)
                , solver(solver_)
                , result_valid(false)
            {}

            void apply();

            void apply(TrialVector& u_);

            const Result& result() const
            {
                if (!result_valid)
                    DUNE_THROW(NewtonError,
                               "NewtonSolver::result() called before NewtonSolver::solve()");
                return this->res;
            }

        protected:
            virtual void defect(TestVector& r)
            {
                r = 0.0;                                        // TODO: vector interface
                this->gridoperator.residual(*this->u, r);
                this->res.defect = this->solver.norm(r);                    // TODO: solver interface
                if (!std::isfinite(this->res.defect))
                    DUNE_THROW(NewtonDefectError,
                               "NewtonSolver::defect(): Non-linear defect is NaN or Inf");
            }


        private:
            void linearSolve(const Matrix& A, TrialVector& z, TestVector& r) const
            {
                if (this->verbosity_level >= 4)
                    std::cout << "      Solving linear system..." << std::endl;
                z = 0.0;                                        // TODO: vector interface
                this->solver.apply(A, z, r, this->linear_reduction);        // TODO: solver interface
                
                IosFlagsRestorer restorer(std::cout); // store old ios flags
                
                if (!this->solver.result().converged)                 // TODO: solver interface
                    DUNE_THROW(NewtonLinearSolverError,
                               "NewtonSolver::linearSolve(): Linear solver did not converge "
                               "in iteration" << this->res.iterations);
                if (this->verbosity_level >= 4)
                    std::cout << "          linear solver iterations:     "
                              << std::setw(12) << solver.result().iterations << std::endl
                              << "          linear defect reduction:      "
                              << std::setw(12) << std::setprecision(4) << std::scientific
                              << solver.result().reduction << std::endl;
            }

            Solver& solver;
            bool result_valid;
        };

        template<class GOS, class S, class TrlV, class TstV>
        void NewtonSolver<GOS,S,TrlV,TstV>::apply(TrlV& u_)
        {
            this->u = &u_;
            apply();
        }

        template<class GOS, class S, class TrlV, class TstV>
        void NewtonSolver<GOS,S,TrlV,TstV>::apply()
        {
            this->res.iterations = 0;
            this->res.converged = false;
            this->res.reduction = 1.0;
            this->res.conv_rate = 1.0;
            this->res.elapsed = 0.0;
            result_valid = true;
            Timer timer;

            try
            {
                TestVector r(this->gridoperator.testGridFunctionSpace());
                this->defect(r);
                this->res.first_defect = this->res.defect;
                this->prev_defect = this->res.defect;

                                

                if (this->verbosity_level >= 2){
                    IosFlagsRestorer restorer(std::cout); // store old ios flags
                    std::cout << "  Initial defect: "
                              << std::setw(12) << std::setprecision(4) << std::scientific
                              << this->res.defect << std::endl;
                }
                
                Matrix A(this->gridoperator);
                TrialVector z(this->gridoperator.trialGridFunctionSpace());

                while (!this->terminate())
                {
                    if (this->verbosity_level >= 3)
                        std::cout << "  Newton iteration " << this->res.iterations
                                  << " --------------------------------" << std::endl;

                    this->prepare_step(A);

                    this->linearSolve(A, z, r);

                    try
                    {
                        this->line_search(z, r);
                    }
                    catch (NewtonLineSearchError)
                    {
                        if (this->reassembled)
                            throw;
                        if (this->verbosity_level >= 3)
                            std::cout << "      line search failed - try again with reassembled matrix" << std::endl;
                        continue;
                    }

                    this->res.reduction = this->res.defect/this->res.first_defect;
                    this->res.iterations++;
                    this->res.conv_rate = std::pow(this->res.reduction, 1.0/this->res.iterations);

                    IosFlagsRestorer restorer(std::cout); // store old ios flags
                                    
                    if (this->verbosity_level >= 3)
                        std::cout << "      defect reduction (this iteration):"
                                  << std::setw(12) << std::setprecision(4) << std::scientific
                                  << this->res.defect/this->prev_defect << std::endl
                                  << "      defect reduction (total):         "
                                  << std::setw(12) << std::setprecision(4) << std::scientific
                                  << this->res.reduction << std::endl
                                  << "      new defect:                       "
                                  << std::setw(12) << std::setprecision(4) << std::scientific
                                  << this->res.defect << std::endl;
                    if (this->verbosity_level == 2)
                        std::cout << "  Newton iteration " << std::setw(2) << this->res.iterations
                                  << ".  Reduction (this): "
                                  << std::setw(12) << std::setprecision(4) << std::scientific
                                  << this->res.defect/this->prev_defect
                                  << ".  Reduction (total): "
                                  << std::setw(12) << std::setprecision(4) << std::scientific
                                  << this->res.reduction << std::endl;
                }
            }
            catch(...)
            {
                this->res.elapsed = timer.elapsed();
                throw;
            }
            this->res.elapsed = timer.elapsed();
            
            IosFlagsRestorer restorer(std::cout); // store old ios flags

            if (this->verbosity_level == 1)
                std::cout << "  Newton converged after " << std::setw(2) << this->res.iterations
                          << " iterations.  Reduction: "
                          << std::setw(12) << std::setprecision(4) << std::scientific
                          << this->res.reduction
                          << "   (" << std::setprecision(4) << this->res.elapsed << "s)"
                          << std::endl;
        }

        template<class GOS, class TrlV, class TstV>
        class NewtonTerminate : public virtual NewtonBase<GOS,TrlV,TstV>
        {
            typedef GOS GridOperator;
            typedef TrlV TrialVector;

            typedef typename TstV::ElementType RFType;

        public:
            NewtonTerminate(GridOperator& go, TrialVector& u_)
                : NewtonBase<GOS,TrlV,TstV>(go,u_)
                , reduction(1e-8)
                , maxit(40)
                , force_iteration(false)
                , abs_limit(1e-12)
            {}

            NewtonTerminate(GridOperator& go)
                : NewtonBase<GOS,TrlV,TstV>(go)
                , reduction(1e-8)
                , maxit(40)
                , force_iteration(false)
                , abs_limit(1e-12)
            {}

            void setReduction(RFType reduction_)
            {
                reduction = reduction_;
            }

            void setMaxIterations(unsigned int maxit_)
            {
                maxit = maxit_;
            }

            void setForceIteration(bool force_iteration_)
            {
                force_iteration = force_iteration_;
            }

            void setAbsoluteLimit(RFType abs_limit_)
            {
                abs_limit = abs_limit_;
            }

            virtual bool terminate()
            {
                if (force_iteration && this->res.iterations == 0)
                    return false;
                this->res.converged = this->res.defect < abs_limit
                    || this->res.defect < this->res.first_defect * reduction;
                if (this->res.iterations >= maxit && !this->res.converged)
                    DUNE_THROW(NewtonNotConverged,
                               "NewtonTerminate::terminate(): Maximum iteration count reached");
                return this->res.converged;
            }

        private:
            RFType reduction;
            unsigned int maxit;
            bool force_iteration;
            RFType abs_limit;
        };

        template<class GOS, class TrlV, class TstV>
        class NewtonPrepareStep : public virtual NewtonBase<GOS,TrlV,TstV>
        {
            typedef GOS GridOperator;
            typedef TrlV TrialVector;

            typedef typename TstV::ElementType RFType;
            typedef typename GOS::template MatrixContainer<RFType>::Type Matrix;

        public:
            NewtonPrepareStep(GridOperator& go, TrialVector& u_)
                : NewtonBase<GOS,TrlV,TstV>(go,u_)
                , min_linear_reduction(1e-3)
                , reassemble_threshold(0.0)
            {}

            NewtonPrepareStep(GridOperator& go)
                : NewtonBase<GOS,TrlV,TstV>(go)
                , min_linear_reduction(1e-3)
                , reassemble_threshold(0.0)
            {}

            void setMinLinearReduction(RFType min_linear_reduction_)
            {
                min_linear_reduction = min_linear_reduction_;
            }

            void setReassembleThreshold(RFType reassemble_threshold_)
            {
                reassemble_threshold = reassemble_threshold_;
            }

            virtual void prepare_step(Matrix& A)
            {
                this->reassembled = false;
                if (this->res.defect/this->prev_defect > reassemble_threshold)
                {
                    if (this->verbosity_level >= 3)
                        std::cout << "      Reassembling matrix..." << std::endl;
                    A = 0.0;                                    // TODO: Matrix interface
                    this->gridoperator.jacobian(*this->u, A);
                    this->reassembled = true;
                }

                this->linear_reduction = std::min(min_linear_reduction,
                                                   this->res.defect*this->res.defect/
                                                   (this->prev_defect*this->prev_defect));

                this->prev_defect = this->res.defect;

                IosFlagsRestorer restorer(std::cout); // store old ios flags

                if (this->verbosity_level >= 4)
                    std::cout << "      requested linear reduction:       "
                              << std::setw(12) << std::setprecision(4) << std::scientific
                              << this->linear_reduction << std::endl;
            }

        private:
            RFType min_linear_reduction;
            RFType reassemble_threshold;
        };

        template<class GOS, class TrlV, class TstV>
        class NewtonLineSearch : public virtual NewtonBase<GOS,TrlV,TstV>
        {
            typedef GOS GridOperator;
            typedef TrlV TrialVector;
            typedef TstV TestVector;

            typedef typename TestVector::ElementType RFType;

        public:
            enum Strategy { noLineSearch,
                            hackbuschReusken,
                            hackbuschReuskenAcceptBest };

            NewtonLineSearch(GridOperator& go, TrialVector& u_)
                : NewtonBase<GOS,TrlV,TstV>(go,u_)
                , strategy(hackbuschReusken)
                , maxit(10)
                , damping_factor(0.5)
            {}

            NewtonLineSearch(GridOperator& go)
                : NewtonBase<GOS,TrlV,TstV>(go)
                , strategy(hackbuschReusken)
                , maxit(10)
                , damping_factor(0.5)
            {}

            void setLineSearchStrategy(Strategy strategy_)
            {
                strategy = strategy_;
            }

            void setLineSearchMaxIterations(unsigned int maxit_)
            {
                maxit = maxit_;
            }

            void setLineSearchDampingFactor(RFType damping_factor_)
            {
                damping_factor = damping_factor_;
            }

            virtual void line_search(TrialVector& z, TestVector& r)
            {
                if (strategy == noLineSearch)
                {
                    this->u->axpy(-1.0, z);                     // TODO: vector interface
                    this->defect(r);
                    return;
                }

                if (this->verbosity_level >= 4)
                    std::cout << "      Performing line search..." << std::endl;
                RFType lambda = 1.0;
                RFType best_lambda = 0.0;
                RFType best_defect = this->res.defect;
                TrialVector prev_u(*this->u);  // TODO: vector interface
                unsigned int i = 0;
                IosFlagsRestorer restorer(std::cout); // store old ios flags

                while (1)
                {
                    if (this->verbosity_level >= 4)
                        std::cout << "          trying line search damping factor:   "
                                  << std::setw(12) << std::setprecision(4) << std::scientific
                                  << lambda 
                                  << std::endl;

                    this->u->axpy(-lambda, z);                  // TODO: vector interface
                    try { 
                        this->defect(r); 
                    }
                    catch (NewtonDefectError) 
                        {
                            if (this->verbosity_level >= 4)
                                std::cout << "          Nans detected" << std::endl;
                        }       // ignore NaNs and try again with lower lambda

                    if (this->res.defect <= (1.0 - lambda/4) * this->prev_defect)
                        {
                            if (this->verbosity_level >= 4)
                                std::cout << "          line search converged" << std::endl;
                            break;
                        }

                    if (this->res.defect < best_defect)
                    {
                        best_defect = this->res.defect;
                        best_lambda = lambda;
                    }

                    if (++i >= maxit)
                    {
                        if (this->verbosity_level >= 4)
                            std::cout << "          max line search iterations exceeded" << std::endl;
                        switch (strategy)
                        {
                        case hackbuschReusken:
                            *this->u = prev_u;
                            this->defect(r);
                            DUNE_THROW(NewtonLineSearchError,
                                       "NewtonLineSearch::line_search(): line search failed, "
                                       "max iteration count reached, "
                                       "defect did not improve enough");
                        case hackbuschReuskenAcceptBest:
                            if (best_lambda == 0.0)
                            {
                                *this->u = prev_u;
                                this->defect(r);
                                DUNE_THROW(NewtonLineSearchError,
                                           "NewtonLineSearch::line_search(): line search failed, "
                                           "max iteration count reached, "
                                           "defect did not improve in any of the iterations");
                            }
                            if (best_lambda != lambda)
                            {
                                *this->u = prev_u;
                                this->u->axpy(-best_lambda, z);
                                this->defect(r);
                            }
                            break;
                        case noLineSearch:
                            break;
                        }
                        break;
                    }

                    lambda *= damping_factor;
                    *this->u = prev_u;                          // TODO: vector interface
                }
                if (this->verbosity_level >= 4)
                    std::cout << "          line search damping factor:   "
                              << std::setw(12) << std::setprecision(4) << std::scientific
                              << lambda << std::endl;
            }

        private:
            Strategy strategy;
            unsigned int maxit;
            RFType damping_factor;
        };

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
            Newton(GridOperator& go, TrialVector& u_, Solver& solver_)
                : NewtonBase<GOS,TrlV,TstV>(go,u_)
                , NewtonSolver<GOS,S,TrlV,TstV>(go,u_,solver_)
                , NewtonTerminate<GOS,TrlV,TstV>(go,u_)
                , NewtonLineSearch<GOS,TrlV,TstV>(go,u_)
                , NewtonPrepareStep<GOS,TrlV,TstV>(go,u_)
            {}
            Newton(GridOperator& go, Solver& solver_)
                : NewtonBase<GOS,TrlV,TstV>(go)
                , NewtonSolver<GOS,S,TrlV,TstV>(go,solver_)
                , NewtonTerminate<GOS,TrlV,TstV>(go)
                , NewtonLineSearch<GOS,TrlV,TstV>(go)
                , NewtonPrepareStep<GOS,TrlV,TstV>(go)
            {}
        };
    }
}

#endif // DUNE_PDELAB_NEWTON_HH
