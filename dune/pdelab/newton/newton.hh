#include <dune/common/exceptions.hh>

namespace Dune
{
    namespace PDELab
    {
        // Exception classes used in NewtonSolver
        class NewtonError : public Exception {};
        class NewtonDefectError : public NewtonError {};
        class NewtonLinearSolverError : public NewtonError {};
        class NewtonLineSearchError : public NewtonError {};

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

        // Traits class for NewtonSolver
        template<class T, class PS, class LS>
        struct NewtonSolverTraits
        {
            typedef T Terminate;
            typedef PS PrepareStep;
            typedef LS LineSearch;
        };
        
        // Default traits for NewtonSolver
        template<class RFType> class NewtonTerminate;
        template<class RFType> class NewtonPrepareStep;
        template<class RFType> class NewtonLineSearch;

        template<class TestVector>
        struct NewtonSolverDefaultTraits
        {
            typedef typename TestVector::ElementType RFType;

            typedef NewtonTerminate<RFType> Terminate;
            typedef NewtonPrepareStep<RFType> PrepareStep;
            typedef NewtonLineSearch<RFType> LineSearch;
        };

        template<class GOS, class S, class TrlV, class TstV = TrlV,
                 class NST = NewtonSolverDefaultTraits<TstV> >
        class NewtonSolver
        {
            typedef GOS GridOperator;
            typedef S Solver;
            typedef TrlV TrialVector;
            typedef TstV TestVector;
            typedef NST Traits;

            typedef typename Traits::Terminate Terminate;
            typedef typename Traits::PrepareStep PrepareStep;
            typedef typename Traits::LineSearch LineSearch;

            friend class Traits::Terminate;
            friend class Traits::PrepareStep;
            friend class Traits::LineSearch;

            typedef typename TestVector::ElementType RFType;

            typedef typename GOS::template MatrixContainer<RFType>::Type Matrix;

        public:
            typedef NewtonResult<RFType> Result;

            NewtonSolver(GridOperator& go, TrialVector& u_, Solver& solver_,
                         const Terminate& terminate_ = Terminate(),
                         const PrepareStep& prepare_step_ = PrepareStep(),
                         const LineSearch& line_search_ = LineSearch())
                : gridoperator(go), u(u_), solver(solver_), terminate(terminate_),
                  prepare_step(prepare_step_), line_search(line_search_),
                  result_valid(false), verbosity_level(1)
            {}

            void apply();

            const Result& result() const
            {
                if (!result_valid)
                    DUNE_THROW(NewtonError,
                               "NewtonSolver::result() called before NewtonSolver::solve()");
                return res;
            }

        private:
            void defect(TestVector& r)
            {
                r = 0.0;                                        // TODO: vector interface
                gridoperator.residual(u, r);
                res.defect = solver.norm(r);                    // TODO: solver interface
                if (!std::isfinite(res.defect))
                    DUNE_THROW(NewtonDefectError,
                               "NewtonSolver::defect(): Non-linear defect is NaN or Inf");
            }

            void linearSolve(const Matrix& A, TrialVector& z, TestVector& r) const
            {
                if (verbosity_level >= 1)
                    std::cout << "      Solving linear system..." << std::endl;
                solver.apply(A, z, r, linear_reduction);        // TODO: solver interface
                if (!solver.result().converged)                 // TODO: solver interface
                    DUNE_THROW(NewtonLinearSolverError,
                               "NewtonSolver::linearSolve(): Linear solver did not converge "
                               "in iteration" << res.iterations);
            }

            GridOperator& gridoperator;
            TrialVector& u;
            Solver& solver;
            Terminate terminate;
            PrepareStep prepare_step;
            LineSearch line_search;
            Result res;
            bool result_valid;
            unsigned int verbosity_level;

            RFType prev_defect;
            RFType linear_reduction;
        };

        template<class GOS, class S, class TrlV, class TstV, class NST>
        void NewtonSolver<GOS,S,TrlV,TstV,NST>::apply()
        {
            res.iterations = 0;
            res.converged = false;
            res.reduction = 1.0;
            res.conv_rate = 1.0;
            res.elapsed = 0.0;
            result_valid = true;
            Timer timer;

            try
            {
                TestVector r(gridoperator.testGridFunctionSpace());
                defect(r);
                res.first_defect = res.defect;
                prev_defect = res.defect;

                if (verbosity_level >= 1)
                    std::cout << "  Initial defect: "
                              << std::setw(12) << std::setprecision(4) << std::scientific
                              << res.defect << std::endl;

                Matrix A(gridoperator);
                TrialVector z(gridoperator.trialGridFunctionSpace());

                while (!terminate(*this))
                {
                    if (verbosity_level >= 1)
                        std::cout << "  Newton iteration " << res.iterations
                                  << " --------------------------" << std::endl;

                    prepare_step(*this, A);

                    linearSolve(A, z, r);

                    prev_defect = res.defect;

                    line_search(*this, z, r);

                    res.reduction = res.defect/res.first_defect;
                    res.iterations++;
                    res.conv_rate = std::pow(res.reduction, 1.0/res.iterations);

                    if (verbosity_level >= 1)
                        std::cout << "      defect reduction (this iteration):"
                                  << std::setw(12) << std::setprecision(4) << std::scientific
                                  << res.defect/prev_defect << std::endl
                                  << "      defect reduction (total):         "
                                  << std::setw(12) << std::setprecision(4) << std::scientific
                                  << res.reduction << std::endl;
                }
            }
            catch(...)
            {
                res.elapsed = timer.elapsed();
                throw;
            }
            res.elapsed = timer.elapsed();
        }

        template<class RFType>
        class NewtonTerminate
        {
        public:
            NewtonTerminate(RFType reduction_ = 1e-8, unsigned int maxit_ = 40,
                            bool force_iteration_ = false, RFType abs_limit_ = 1e-12)
                : reduction(reduction_), maxit(maxit_),
                  force_iteration(force_iteration_), abs_limit(abs_limit_)
            {}

            template<class Newton>
            bool operator()(Newton& newton)
            {
                if (force_iteration && newton.res.iterations == 0)
                    return false;
                newton.res.converged = newton.res.defect < abs_limit
                    || newton.res.defect < newton.res.first_defect * reduction;
                return newton.res.iterations >= maxit || newton.res.converged;
            }

        private:
            RFType reduction;
            unsigned int maxit;
            bool force_iteration;
            RFType abs_limit;
        };

        template<class RFType>
        class NewtonPrepareStep
        {
        public:
            NewtonPrepareStep(RFType min_linear_reduction_ = 1e-3,
                              RFType reassemble_threshold_ = 0.8)
                : min_linear_reduction(min_linear_reduction_),
                  reassemble_threshold(reassemble_threshold_)
            {}

            template<class Newton>
            void operator()(Newton& newton, typename Newton::Matrix& A)
            {
                if (newton.res.defect/newton.prev_defect > reassemble_threshold)
                {
                    if (newton.verbosity_level >= 1)
                        std::cout << "      Reassembling matrix..." << std::endl;
                    A = 0.0;                                    // TODO: Matrix interface
                    newton.gridoperator.jacobian(newton.u, A);
                }

                newton.linear_reduction = std::min(min_linear_reduction,
                                                   newton.res.defect*newton.res.defect/
                                                   (newton.prev_defect*newton.prev_defect));
                if (newton.verbosity_level >= 1)
                    std::cout << "      linear reduction:                 "
                              << std::setw(12) << std::setprecision(4) << std::scientific
                              << newton.linear_reduction << std::endl;
            }

        private:
            RFType min_linear_reduction;
            RFType reassemble_threshold;
        };

        template<class RFType>
        class NewtonLineSearch
        {
        public:
            enum Strategy { noLineSearch,
                            hackbuschReusken,
                            hackbuschReuskenAcceptBest };

            NewtonLineSearch(Strategy strategy_ = hackbuschReusken, unsigned int maxit_ = 10,
                             RFType damping_factor_ = 0.5)
                : strategy(strategy_), maxit(maxit_), damping_factor(damping_factor_)
            {}

            template<class Newton>
            void operator()(Newton& newton, typename Newton::TrialVector& z,
                            typename Newton::TestVector& r)
            {
                if (strategy == noLineSearch)
                {
                    newton.u.axpy(-1.0, z);                     // TODO: vector interface
                    newton.defect(r);
                    return;
                }

                RFType lambda = 1.0;
                RFType best_lambda = 0.0;
                RFType best_defect = newton.res.defect;
                typename Newton::TrialVector prev_u(newton.u);  // TODO: vector interface
                unsigned int i = 0;
                while (1)
                {
                    newton.u.axpy(-lambda, z);                  // TODO: vector interface
                    try { newton.defect(r); }
                    catch (NewtonDefectError) {}

                    if (newton.res.defect <= (1.0 - lambda/4) * newton.prev_defect)
                        break;

                    if (newton.res.defect < best_defect)
                    {
                        best_defect = newton.res.defect;
                        best_lambda = lambda;
                    }

                    if (++i >= maxit)
                    {
                        switch (strategy)
                        {
                        case hackbuschReusken:
                            newton.u = prev_u;
                            newton.defect(r);
                            DUNE_THROW(NewtonLineSearchError,
                                       "NewtonLineSearch::operator(): line search failed, "
                                       "max iteration count reached, "
                                       "defect did not improve enough");
                        case hackbuschReuskenAcceptBest:
                            if (best_lambda == 0.0)
                                DUNE_THROW(NewtonLineSearchError,
                                           "NewtonLineSearch::operator(): line search failed, "
                                           "max iteration count reached, "
                                           "defect did not improve in any of the iterations");
                            if (best_lambda != lambda)
                            {
                                newton.u = prev_u;
                                newton.u.axpy(-best_lambda, z);
                                newton.defect(r);
                            }
                            break;
                        case noLineSearch:
                            break;
                        }
                        break;
                    }

                    lambda *= damping_factor;
                    newton.u = prev_u;                          // TODO: vector interface
                }
            }

        private:
            Strategy strategy;
            unsigned int maxit;
            RFType damping_factor;
        };
    }
}
