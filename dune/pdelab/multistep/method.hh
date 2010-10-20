// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_MULTISTEP_METHOD_HH
#define DUNE_PDELAB_MULTISTEP_METHOD_HH

#include <iomanip>
#include <iostream>

#include <dune/common/ios_state.hh>
#include <dune/common/shared_ptr.hh>

#include <dune/pdelab/multistep/parameter.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup MultiStepMethods Multi-Step Methods
    //! \ingroup PDELab
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
        if(mgos.trialGridFunctionSpace().gridview().comm().rank() > 0)
          verbosity = 0;
      }

      //! change verbosity level; 0 means completely quiet
      void setVerbosityLevel (int level)
      {
        if(mgos.trialGridFunctionSpace().gridview().comm().rank() == 0)
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
        mgos.preStep(time,dt, oldValues);

        // solve
        pdeSolver.apply(xnew);

        // postprocessing in the assembler
        mgos.postStep();

        ++step;

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
        mgos.preStep(time, dt, oldValues);

        // set boundary conditions and initial value
        f.setTime(time+dt);
        mgos.interpolate(*oldValues[0],f,xnew);

        // solve stage
        pdeSolver.apply(xnew);

        // postprocessing in the assembler
        mgos.postStep();

        step++;
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
      shared_ptr<const TrialV> apply(T time, T dt)
      {
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

        // create vector, using last time step as start
        shared_ptr<TrialV>
          xnew(new TrialV(*mgos.getCache()->getUnknowns(step-1)));

        // prepare assembler
        mgos.preStep(step, time, dt);

        // solve
        pdeSolver.apply(*xnew);

        // postprocessing in the assembler
        mgos.postStep();

        // store result for next step
        mgos.getCache()->setUnknowns(step, xnew);

        ++step;

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
      shared_ptr<const TrialV> apply(T time, T dt, F& f)
      {
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

        // create vector
        shared_ptr<TrialV> xnew(new TrialV(mgos.trialGridFunctionSpace()));

        // prepare assembler
        mgos.preStep(step, time, dt);

        // set boundary conditions and initial value
        f.setTime(time+dt);
        mgos.interpolate(*mgos.getCache()->getUnknowns(step-1),f,*xnew);

        // solve stage
        pdeSolver.apply(*xnew);

        // postprocessing in the assembler
        mgos.postStep();

        // store result for next step
        mgos.getCache()->setUnknowns(step, xnew);

        ++step;

        return xnew;
      }
    };

    //! \} group MultiStepMethods
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_MULTISTEP_METHOD_HH
