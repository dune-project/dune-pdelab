// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_MULTISTEP_METHOD_HH
#define DUNE_PDELAB_MULTISTEP_METHOD_HH

#include <iomanip>
#include <iostream>
#include <memory>

#include <dune/common/ios_state.hh>
#include <dune/common/timer.hh>

#include <dune/pdelab/multistep/parameter.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup MultiStepMethods
    //! \{

    //! Do one step of a multi-step time-stepping scheme
    /**
     * \tparam T          Type to represent time values
     * \tparam MGOS       Assembler for multi-step instationary problems
     * \tparam PDESOLVER  Solver problem in each step (typically Newton)
     * \tparam TrialV     Vector type to represent coefficients of solutions
     * \tparam TestV      Vector type to represent residuals
     */
    template<class T, class MGOS, class PDESolver, class TrialV,
             class TestV = TrialV>
    class MultiStepMethod {
      typedef MultiStepParameterInterface<T, MGOS::order> Parameters;

      const Parameters* parameters;
      MGOS& mgos;
      PDESolver& pdeSolver;
      unsigned verbosity;
      unsigned step;

    public:
      //! construct a new multi-step scheme
      /**
       * \param parameters_ Parameter object.
       * \param mgos_       Assembler object (MultiStepGridOperatorSpace).
       * \param pdesolver_  Solver object (typically Newton).
       *
       * The contructed method object stores references to the object it is
       * constructed with, so these objects should be valid for as long as the
       * constructed object is used.
       */
      MultiStepMethod(const Parameters& parameters_,
                      MGOS& mgos_, PDESolver& pdeSolver_)
        : parameters(&parameters_), mgos(mgos_), pdeSolver(pdeSolver_),
          verbosity(1), step(1)
      {
        if(mgos.trialGridFunctionSpace().gridView().comm().rank() > 0)
          verbosity = 0;
      }

      //! change verbosity level; 0 means completely quiet
      void setVerbosityLevel (int level)
      {
        if(mgos.trialGridFunctionSpace().gridView().comm().rank() == 0)
          verbosity = level;
      }

      //! redefine the method to be used; can be done before every step
      void setMethod(const Parameters &parameters_) {
        parameters = &parameters_;
      }

      //! do one step;
      /*
       * \param time      Start of time step
       * \param dt        Suggested time step size
       * \param oldValues Vector of pointers to the old values.  Must support
       *                  the expression *oldvalues[i], which should yield a
       *                  reference to the old value at time time-i*dt.
       * \param xnew      Where to store the new value.
       *
       * \return Actual time step size
       */
      template<class OldValues>
      T apply(T time, T dt, const OldValues& oldValues, TrialV& xnew)
      {
        Timer allTimer;
        Timer subTimer;

        // save formatting attributes
        ios_base_all_saver format_attribute_saver(std::cout);

        if(verbosity >= 1)
          std::cout << "TIME STEP [" << parameters->name() << "] "
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

        // prepare assembler
        if(verbosity >= 2) {
          std::cout << "== prepare assembler" << std::endl;
          subTimer.reset();
        }
        mgos.preStep(time,dt, oldValues);
        if(verbosity >= 2)
          std::cout << "== prepare assembler (" << subTimer.elapsed() << "s)"
                    << std::endl;

        // solve
        if(verbosity >= 2) {
          std::cout << "== apply solver" << std::endl;
          subTimer.reset();
        }
        pdeSolver.apply(xnew);
        if(verbosity >= 2)
          std::cout << "== apply solver (" << subTimer.elapsed() << "s)"
                    << std::endl;

        // postprocessing in the assembler
        if(verbosity >= 2) {
          std::cout << "== cleanup assembler" << std::endl;
          subTimer.reset();
        }
        mgos.postStep();
        if(verbosity >= 2)
          std::cout << "== cleanup assembler (" << subTimer.elapsed() << "s)"
                    << std::endl;

        ++step;

        if(verbosity >= 2)
          std::cout << "== time step done (" << allTimer.elapsed() << "s)"
                    << std::endl;

        return dt;
      }

      //! do one step;
      /* This is a version which interpolates constraints at the start of each
       * stage
       *
       * \param time      Start of time step
       * \param dt        Suggested time step size
       * \param oldValues Vector of pointers to the old values.  Must support
       *                  the expression *oldValues[i], which should yield a
       *                  reference to the old value at time time-i*dt.
       * \param f         Function to interpolate boundary conditions from.
       *                  Should support the method setTime().
       * \param xnew      Where to store the new value.
       *
       * \return Actual time step size
       */
      template<typename OldValues, typename F>
      T apply(T time, T dt, const OldValues& oldValues, F& f, TrialV& xnew)
      {
        Timer allTimer;
        Timer subTimer;

        // save formatting attributes
        ios_base_all_saver format_attribute_saver(std::cout);

        if(verbosity >= 1)
          std::cout << "TIME STEP [" << parameters->name() << "] "
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

        // prepare assembler
        if(verbosity >= 2) {
          std::cout << "== prepare assembler" << std::endl;
          subTimer.reset();
        }
        mgos.preStep(time, dt, oldValues);
        if(verbosity >= 2)
          std::cout << "== prepare assembler (" << subTimer.elapsed() << "s)"
                    << std::endl;

        // set boundary conditions and initial value
        if(verbosity >= 2) {
          std::cout << "== setup result vector" << std::endl;
          subTimer.reset();
        }
        f.setTime(time+dt);
        mgos.interpolate(*oldValues[0],f,xnew);
        if(verbosity >= 2)
          std::cout << "== setup result vector (" << subTimer.elapsed() << "s)"
                    << std::endl;

        // solve stage
        if(verbosity >= 2) {
          std::cout << "== apply solver" << std::endl;
          subTimer.reset();
        }
        pdeSolver.apply(xnew);
        if(verbosity >= 2)
          std::cout << "== apply solver (" << subTimer.elapsed() << "s)"
                    << std::endl;

        // postprocessing in the assembler
        if(verbosity >= 2) {
          std::cout << "== cleanup assembler" << std::endl;
          subTimer.reset();
        }
        mgos.postStep();
        if(verbosity >= 2)
          std::cout << "== cleanup assembler (" << subTimer.elapsed() << "s)"
                    << std::endl;

        step++;

        if(verbosity >= 2)
          std::cout << "== time step done (" << allTimer.elapsed() << "s)"
                    << std::endl;

        return dt;
      }

      //! do one step (with caching)
      /**
       * \param time Start of time step
       * \param dt   Time step size
       *
       * \return A shared_ptr to the new value
       *
       * The old values are expected in the cache of the GridOperatorSpace.
       * The computed value is store in the cache as well.
       */
      std::shared_ptr<const TrialV> apply(T time, T dt)
      {
        Timer allTimer;
        Timer subTimer;

        // save formatting attributes
        ios_base_all_saver format_attribute_saver(std::cout);

        if(verbosity >= 1)
          std::cout << "TIME STEP [" << parameters->name() << "] "
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

        // prepare assembler
        if(verbosity >= 2) {
          std::cout << "== prepare assembler" << std::endl;
          subTimer.reset();
        }
        mgos.preStep(step, time, dt);
        if(verbosity >= 2)
          std::cout << "== prepare assembler (" << subTimer.elapsed() << "s)"
                    << std::endl;

        // create vector, using last time step as start
        if(verbosity >= 2) {
          std::cout << "== setup result vector" << std::endl;
          subTimer.reset();
        }
        std::shared_ptr<TrialV>
          xnew(new TrialV(*mgos.getCache()->getUnknowns(step-1)));
        if(verbosity >= 2)
          std::cout << "== setup result vector (" << subTimer.elapsed() << "s)"
                    << std::endl;

        // solve
        if(verbosity >= 2) {
          std::cout << "== apply solver" << std::endl;
          subTimer.reset();
        }
        pdeSolver.apply(*xnew);
        if(verbosity >= 2)
          std::cout << "== apply solver (" << subTimer.elapsed() << "s)"
                    << std::endl;

        // postprocessing in the assembler
        if(verbosity >= 2) {
          std::cout << "== cleanup assembler" << std::endl;
          subTimer.reset();
        }
        mgos.postStep();
        if(verbosity >= 2)
          std::cout << "== cleanup assembler (" << subTimer.elapsed() << "s)"
                    << std::endl;

        // store result for next step
        mgos.getCache()->setUnknowns(step, xnew);

        ++step;

        if(verbosity >= 2)
          std::cout << "== time step done (" << allTimer.elapsed() << "s)"
                    << std::endl;

        return xnew;
      }

      //! do one step (with caching)
      /**
       * This is a version which interpolates constraints at the start of each
       * stage
       *
       * \param time Start of time step
       * \param dt   Time step size
       * \param f    Function to interpolate boundary conditions from.
       *             Should support the method setTime().
       *
       * \return A shared_ptr to the new value
       *
       * The old values are expected in the cache of the GridOperatorSpace.
       * The computed value is store in the cache as well.
       */
      template<typename F>
      std::shared_ptr<const TrialV> apply(T time, T dt, F& f)
      {
        Timer allTimer;
        Timer subTimer;

        // save formatting attributes
        ios_base_all_saver format_attribute_saver(std::cout);

        if(verbosity >= 1)
          std::cout << "TIME STEP [" << parameters->name() << "] "
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

        // prepare assembler
        if(verbosity >= 2) {
          std::cout << "== prepare assembler" << std::endl;
          subTimer.reset();
        }
        mgos.preStep(step, time, dt);
        if(verbosity >= 2)
          std::cout << "== prepare assembler (" << subTimer.elapsed() << "s)"
                    << std::endl;

        // setup vector
        if(verbosity >= 2) {
          std::cout << "== setup result vector" << std::endl;
          subTimer.reset();
        }
        std::shared_ptr<TrialV> xnew(new TrialV(mgos.trialGridFunctionSpace()));
        // set boundary conditions and initial value
        f.setTime(time+dt);
        mgos.interpolate(*mgos.getCache()->getUnknowns(step-1),f,*xnew);
        if(verbosity >= 2)
          std::cout << "== setup result vector (" << subTimer.elapsed() << "s)"
                    << std::endl;

        // solve stage
        if(verbosity >= 2) {
          std::cout << "== apply solver" << std::endl;
          subTimer.reset();
        }
        pdeSolver.apply(*xnew);
        if(verbosity >= 2)
          std::cout << "== apply solver (" << subTimer.elapsed() << "s)"
                    << std::endl;

        // postprocessing in the assembler
        if(verbosity >= 2) {
          std::cout << "== cleanup assembler" << std::endl;
          subTimer.reset();
        }
        mgos.postStep();
        if(verbosity >= 2)
          std::cout << "== cleanup assembler (" << subTimer.elapsed() << "s)"
                    << std::endl;

        // store result for next step
        mgos.getCache()->setUnknowns(step, xnew);

        ++step;

        if(verbosity >= 2)
          std::cout << "== time step done (" << allTimer.elapsed() << "s)"
                    << std::endl;

        return xnew;
      }
    };

    //! \} group MultiStepMethods
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_MULTISTEP_METHOD_HH
