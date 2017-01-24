// -*- tab-width: 2; indent-tabs-mode: nil -*-
// vi: set et ts=2 sw=2 sts=2:

#ifndef DUNE_PDELAB_INSTATIONARY_IMPLICITONESTEP_HH
#define DUNE_PDELAB_INSTATIONARY_IMPLICITONESTEP_HH

#include <iostream>
#include <iomanip>

#include <dune/common/ios_state.hh>
#include <dune/pdelab/instationary/onestepparameter.hh>

namespace Dune {
  namespace PDELab {

    /**
     *  @addtogroup OneStepMethod
     *  @{
     */

    // Status information of Newton's method
    struct OneStepMethodPartialResult
    {
      unsigned int timesteps;
      double assembler_time;
      double linear_solver_time;
      int linear_solver_iterations;
      int nonlinear_solver_iterations;

      OneStepMethodPartialResult() :
        timesteps(0),
        assembler_time(0.0),
        linear_solver_time(0.0),
        linear_solver_iterations(0),
        nonlinear_solver_iterations(0)
      {}
    };

    struct OneStepMethodResult
    {
      OneStepMethodPartialResult total;
      OneStepMethodPartialResult successful;
      OneStepMethodResult() : total(), successful()
      {}
    };

    //! Do one step of a time-stepping scheme
    /**
     * \tparam T          type to represent time values
     * \tparam IGOS       assembler for instationary problems
     * \tparam PDESOLVER  solver problem in each step (typically Newton)
     * \tparam TrlV       vector type to represent coefficients of solutions
     * \tparam TstV       vector type to represent residuals
     */
    template<class T, class IGOS, class PDESOLVER, class TrlV, class TstV = TrlV>
    class OneStepMethod
    {
      typedef typename PDESOLVER::Result PDESolverResult;

    public:
      typedef OneStepMethodResult Result;

      //! construct a new one step scheme
      /**
       * \param method_    Parameter object. This chooses the actual method
       *                   used.
       * \param igos_      Assembler object (instationary grid operator space).
       * \param pdesolver_ solver object (typically Newton).
       *
       * The contructed method object stores references to the object it is
       * constructed with, so these objects should be valid for as long as the
       * constructed object is used (or until setMethod() is called, see
       * there).
       */
      OneStepMethod(const TimeSteppingParameterInterface<T>& method_,
                    IGOS& igos_, PDESOLVER& pdesolver_)
        : method(&method_), igos(igos_), pdesolver(pdesolver_), verbosityLevel(1), step(1), res()
      {
        if (igos.trialGridFunctionSpace().gridView().comm().rank()>0)
          verbosityLevel = 0;
      }

      //! change verbosity level; 0 means completely quiet
      void setVerbosityLevel (int level)
      {
        if (igos.trialGridFunctionSpace().gridView().comm().rank()>0)
          verbosityLevel = 0;
        else
          verbosityLevel = level;
      }

      //! change number of current step
      void setStepNumber(int newstep) { step = newstep; }

      //! Access to the (non) linear solver
      const PDESOLVER & getPDESolver() const
      {
        return pdesolver;
      }

      //! Access to the (non) linear solver
      PDESOLVER & getPDESolver()
      {
        return pdesolver;
      }

      const Result& result() const
      {
        return res;
      }

      //! Set a new result
      /**
       *  \param result_ OneStepMethodResult object
       *
       *  Set the step number to the next timestep according to the result.
       */
      void setResult (const OneStepMethodResult& result_)
      {
        res = result_;
        setStepNumber(res.successful.timesteps+1);
      }

      //! redefine the method to be used; can be done before every step
      /**
       * \param method_ Parameter object.
       *
       * The OneStepMethod object stores a reference to the method_ object.
       * The old method object is no longer referenced after this member
       * function returns.
       */
      void setMethod (const TimeSteppingParameterInterface<T>& method_)
      {
        method = &method_;
      }

      /*! \brief do one step;
       * \param[in]  time start of time step
       * \param[in]  dt suggested time step size
       * \param[in]  xold value at begin of time step
       * \param[in,out] xnew value at end of time step; contains initial guess for first substep on entry
       * \return selected time step size
       */
      T apply (T time, T dt, TrlV& xold, TrlV& xnew)
      {
        // save formatting attributes
        ios_base_all_saver format_attribute_saver(std::cout);

        // do statistics
        OneStepMethodPartialResult step_result;

        std::vector<TrlV*> x(1); // vector of pointers to all steps
        x[0] = &xold;            // initially we have only one

        if (verbosityLevel>=1){
          std::ios_base::fmtflags oldflags = std::cout.flags();
          std::cout << "TIME STEP [" << method->name() << "] "
                    << std::setw(6) << step
                    << " time (from): "
                    << std::setw(12) << std::setprecision(4) << std::scientific
                    << time
                    << " dt: "
                    << std::setw(12) << std::setprecision(4) << std::scientific
                    << dt
                    << " time (to): "
                    << std::setw(12) << std::setprecision(4) << std::scientific
                    << time+dt
                    << std::endl;
          std::cout.flags(oldflags);
        }

        // prepare assembler
        igos.preStep(*method,time,dt);

        // loop over all stages
        for (unsigned r=1; r<=method->s(); ++r)
          {
            if (verbosityLevel>=2){
              std::ios_base::fmtflags oldflags = std::cout.flags();
              std::cout << "STAGE "
                        << r
                        << " time (to): "
                        << std::setw(12) << std::setprecision(4) << std::scientific
                        << time+method->d(r)*dt
                        << "." << std::endl;
              std::cout.flags(oldflags);
            }

            // prepare stage
            igos.preStage(r,x);

            // get vector for current stage
            if (r==method->s())
              {
                // last stage
                x.push_back(&xnew);
                if (r>1) xnew = *(x[r-1]); // if r=1 then xnew has already initial guess
              }
            else
              {
                // intermediate step
                x.push_back(new TrlV(igos.trialGridFunctionSpace()));
                if (r>1)
                  *(x[r]) = *(x[r-1]); // use result of last stage as initial guess
                else
                  *(x[r]) = xnew;
              }

            // solve stage
            try {
              pdesolver.apply(*x[r]);
            }
            catch (...)
              {
                // time step failed -> accumulate to total only
                PDESolverResult pderes = pdesolver.result();
                step_result.assembler_time += pderes.assembler_time;
                step_result.linear_solver_time += pderes.linear_solver_time;
                step_result.linear_solver_iterations += pderes.linear_solver_iterations;
                step_result.nonlinear_solver_iterations += pderes.iterations;
                res.total.assembler_time += step_result.assembler_time;
                res.total.linear_solver_time += step_result.linear_solver_time;
                res.total.linear_solver_iterations += step_result.linear_solver_iterations;
                res.total.nonlinear_solver_iterations += step_result.nonlinear_solver_iterations;
                res.total.timesteps += 1;
                throw;
              }
            PDESolverResult pderes = pdesolver.result();
            step_result.assembler_time += pderes.assembler_time;
            step_result.linear_solver_time += pderes.linear_solver_time;
            step_result.linear_solver_iterations += pderes.linear_solver_iterations;
            step_result.nonlinear_solver_iterations += pderes.iterations;

            // stage cleanup
            igos.postStage();
          }

        // delete intermediate steps
        for (unsigned i=1; i<method->s(); ++i) delete x[i];

        // step cleanup
        igos.postStep();

        // update statistics
        res.total.assembler_time += step_result.assembler_time;
        res.total.linear_solver_time += step_result.linear_solver_time;
        res.total.linear_solver_iterations += step_result.linear_solver_iterations;
        res.total.nonlinear_solver_iterations += step_result.nonlinear_solver_iterations;
        res.total.timesteps += 1;
        res.successful.assembler_time += step_result.assembler_time;
        res.successful.linear_solver_time += step_result.linear_solver_time;
        res.successful.linear_solver_iterations += step_result.linear_solver_iterations;
        res.successful.nonlinear_solver_iterations += step_result.nonlinear_solver_iterations;
        res.successful.timesteps += 1;
        if (verbosityLevel>=1){
          std::ios_base::fmtflags oldflags = std::cout.flags();
          std::cout << "::: timesteps      " << std::setw(6) << res.successful.timesteps
                    << " (" << res.total.timesteps << ")" << std::endl;
          std::cout << "::: nl iterations  " << std::setw(6) << res.successful.nonlinear_solver_iterations
                    << " (" << res.total.nonlinear_solver_iterations << ")" << std::endl;
          std::cout << "::: lin iterations " << std::setw(6) << res.successful.linear_solver_iterations
                    << " (" << res.total.linear_solver_iterations << ")" << std::endl;
          std::cout << "::: assemble time  " << std::setw(12) << std::setprecision(4) << std::scientific
                    << res.successful.assembler_time << " (" << res.total.assembler_time << ")" << std::endl;
          std::cout << "::: lin solve time " << std::setw(12) << std::setprecision(4) << std::scientific
                    << res.successful.linear_solver_time << " (" << res.total.linear_solver_time << ")" << std::endl;
          std::cout.flags(oldflags);
        }

        step++;
        return dt;
      }

      /*! \brief do one step;
       * This is a version which interpolates constraints at the start of each stage
       *
       * \param[in]  time start of time step
       * \param[in]  dt suggested time step size
       * \param[in]  xold value at begin of time step
       * \param[in]  f function to interpolate boundary conditions from
       * \param[in,out] xnew value at end of time step; contains initial guess for first substep on entry
       * \return selected time step size
       */
      template<typename F>
      T apply (T time, T dt, TrlV& xold, F& f, TrlV& xnew)
      {
        // do statistics
        OneStepMethodPartialResult step_result;

        // save formatting attributes
        ios_base_all_saver format_attribute_saver(std::cout);

        std::vector<TrlV*> x(1); // vector of pointers to all steps
        x[0] = &xold;            // initially we have only one

        if (verbosityLevel>=1){
          std::ios_base::fmtflags oldflags = std::cout.flags();
          std::cout << "TIME STEP [" << method->name() << "] "
                    << std::setw(6) << step
                    << " time (from): "
                    << std::setw(12) << std::setprecision(4) << std::scientific
                    << time
                    << " dt: "
                    << std::setw(12) << std::setprecision(4) << std::scientific
                    << dt
                    << " time (to): "
                    << std::setw(12) << std::setprecision(4) << std::scientific
                    << time+dt
                    << std::endl;
          std::cout.flags(oldflags);
        }

        // prepare assembler
        igos.preStep(*method,time,dt);

        // loop over all stages
        for (unsigned r=1; r<=method->s(); ++r)
          {
            if (verbosityLevel>=2){
              std::ios_base::fmtflags oldflags = std::cout.flags();
              std::cout << "STAGE "
                        << r
                        << " time (to): "
                        << std::setw(12) << std::setprecision(4) << std::scientific
                        << time+method->d(r)*dt
                        << "." << std::endl;
              std::cout.flags(oldflags);
            }

            // prepare stage
            igos.preStage(r,x);

            // get vector for current stage
            if (r==method->s())
              {
                // last stage
                x.push_back(&xnew);
              }
            else
              {
                // intermediate step
                x.push_back(new TrlV(igos.trialGridFunctionSpace()));
              }

            // set boundary conditions and initial value
            igos.interpolate(r,*x[r-1],f,*x[r]);

            // solve stage
            try {
              pdesolver.apply(*x[r]);
            }
            catch (...)
              {
                // time step failed -> accumulate to total only
                PDESolverResult pderes = pdesolver.result();
                step_result.assembler_time += pderes.assembler_time;
                step_result.linear_solver_time += pderes.linear_solver_time;
                step_result.linear_solver_iterations += pderes.linear_solver_iterations;
                step_result.nonlinear_solver_iterations += pderes.iterations;
                res.total.assembler_time += step_result.assembler_time;
                res.total.linear_solver_time += step_result.linear_solver_time;
                res.total.linear_solver_iterations += step_result.linear_solver_iterations;
                res.total.nonlinear_solver_iterations += step_result.nonlinear_solver_iterations;
                res.total.timesteps += 1;
                throw;
              }
            PDESolverResult pderes = pdesolver.result();
            step_result.assembler_time += pderes.assembler_time;
            step_result.linear_solver_time += pderes.linear_solver_time;
            step_result.linear_solver_iterations += pderes.linear_solver_iterations;
            step_result.nonlinear_solver_iterations += pderes.iterations;

            // stage cleanup
            igos.postStage();
          }

        // delete intermediate steps
        for (unsigned i=1; i<method->s(); ++i) delete x[i];

        // step cleanup
        igos.postStep();

        // update statistics
        res.total.assembler_time += step_result.assembler_time;
        res.total.linear_solver_time += step_result.linear_solver_time;
        res.total.linear_solver_iterations += step_result.linear_solver_iterations;
        res.total.nonlinear_solver_iterations += step_result.nonlinear_solver_iterations;
        res.total.timesteps += 1;
        res.successful.assembler_time += step_result.assembler_time;
        res.successful.linear_solver_time += step_result.linear_solver_time;
        res.successful.linear_solver_iterations += step_result.linear_solver_iterations;
        res.successful.nonlinear_solver_iterations += step_result.nonlinear_solver_iterations;
        res.successful.timesteps += 1;
        if (verbosityLevel>=1){
          std::ios_base::fmtflags oldflags = std::cout.flags();
          std::cout << "::: timesteps      " << std::setw(6) << res.successful.timesteps
                    << " (" << res.total.timesteps << ")" << std::endl;
          std::cout << "::: nl iterations  " << std::setw(6) << res.successful.nonlinear_solver_iterations
                    << " (" << res.total.nonlinear_solver_iterations << ")" << std::endl;
          std::cout << "::: lin iterations " << std::setw(6) << res.successful.linear_solver_iterations
                    << " (" << res.total.linear_solver_iterations << ")" << std::endl;
          std::cout << "::: assemble time  " << std::setw(12) << std::setprecision(4) << std::scientific
                    << res.successful.assembler_time << " (" << res.total.assembler_time << ")" << std::endl;
          std::cout << "::: lin solve time " << std::setw(12) << std::setprecision(4) << std::scientific
                    << res.successful.linear_solver_time << " (" << res.total.linear_solver_time << ")" << std::endl;
          std::cout.flags(oldflags);
        }

        step++;
        return dt;
      }

    private:
      const TimeSteppingParameterInterface<T> *method;
      IGOS& igos;
      PDESOLVER& pdesolver;
      int verbosityLevel;
      int step;
      Result res;
    };

    /** @} */
  } // end namespace PDELab
} // end namespace Dune
#endif // DUNE_PDELAB_INSTATIONARY_IMPLICITONESTEP_HH
