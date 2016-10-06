// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=8 sw=4 sts=4:
#ifndef DUNE_PDELAB_LOCALOPERATOR_IDEFAULT_HH
#define DUNE_PDELAB_LOCALOPERATOR_IDEFAULT_HH

namespace Dune
{
    namespace PDELab
    {
        //! \addtogroup LocalOperator
        //! \ingroup PDELab
        //! \{

        //! Default class for additional methods in instationary local operators
        /**
         * The algorithm used by the InstationaryGridOperator and the
         * OneStepMethod to advance one step in time is as follows:
         * <ol>
         * <li>Call the method preStep(start_time_of_step, step_size, nstages)
         *     on each localoperator.
         * <li>For each stage1 in [1, nstages]
         *     <ol>
         *     <li>Call preStage(time_of_stage1, number_of_stage1) on each
         *         localoperator.
         *     <li>Assemble constant part of residual: For each stage2 in [0,
         *         number_of_stage1-1]
         *         <ol>
         *         <li>Call setTime(time_of_stage2) for each local operator.
         *         <li>Iterate over grid and evaluate the local operator.
         *         </ol>
         *     <li>Call setTime(time_of_stage1) for each local operator.
         *     <li>Apply solver: Iterate over grid (possibly multiple times
         *         and evaluate local operator.
         *     <li>Call postStage() on each local operator.
         *     </ol>
         * <li>Call postStep() on the temporal local operator.
         * </ol>
         *
         * The algorithm used by the InstationaryGridOperator and the
         * ExplicitOneStepMethod to advance one step in time is as follows:
         * <ol>
         * <li>Call the method preStep(start_time_of_step, step_size, nstages)
         *     on each localoperator.
         * <li>For each stage1 in [1, nstages]
         *     <ol>
         *     <li>Call preStage(time_of_stage1, number_of_stage1) on each
         *         localoperator.
         *     <li>Assemble constant part of residual: For each stage2 in [0,
         *         number_of_stage1-1]
         *         <ol>
         *         <li>Call setTime(time_of_stage2) for each local operator.
         *         <li>Iterate over grid and evaluate the local operator.
         *         </ol>
         *     <li>Call setTime(time_of_stage1) for each local operator.
         *     <li>Assemble matrices to be solved: Iterate over the grid and
         *         evaluate the local operator.
         *     <li>If stage1 == 1:
         *         <ol>
         *         <li>Possibly call suggestTimeStep() and adjust the time
         *             step accordingly.
         *         </ol>
         *     <li>Call postStage() on each local operator.
         *     </ol>
         * <li>Call postStep() on the temporal local operator.
         * </ol>
         *
         * The algorithm used by the MultiStepGridOperator and the
         * MultiStepMethod to advance one step in time is as follows:
         * <ol>
         * <li>Call the method preStep(start_time_of_step, step_size, 1)
         *     on each localoperator.
         * <li>Call the method preStage(end_time_of_step, 1) on each local
         *     operator.
         * <li>Assemble constant part of residual: For each step in [1,
         *     number_of_steps]
         *     <ol>
         *     <li>Call setTime(end_time_of_step-step_size*number_of_step) for
         *         each local operator.
         *     <li>Iterate over grid and evaluate the local operators.
         *     </ol>
         * <li>Call setTime(end_time_of_step) for each local operator.
         * <li>Apply solver: Iterate over grid (possibly multiple times and
         *     evaluate local operator.
         * <li>Call postStage() on each local operator.
         * <li>Call postStep() on each local operator.
         * </ol>
         */
        template<class R = double>
        class InstationaryLocalOperatorDefaultMethods
        {
        public:
            typedef R RealType;

            //! set time for subsequent evaluation
            /**
             * This method set the time for subsequent calls to the alpha_*(),
             * lambda_*(), jacobian_*() and jacobian_apply_*() methods.
             *
             * \note For ExplicitOneStepMethod the time given here in the
             *       first stage may be incorrect, since the time step size is
             *       only finally determined after the first stage has been
             *       assembled.
             */
            void setTime (R t_)
            {
                t = t_;
            }

            //! get current time
            /**
             * \return The time previously set by setTime().
             */
            R getTime () const
            {
                return t;
            }

            //! to be called once before each time step
            /**
             * \param time   Time at beginning of the step.
             * \param dt     Size of time step.
             * \param stages Number of stages to do in the step.  For the
             *               MultiStepMethod this is always 1.
             *
             * \note For ExplicitOneStepMethod the dt given here may be
             *       incorrect, since the time step size is only finally
             *       determined after the first stage has been assembled.
             *
             * \note For the MultiStepMethod the number of stages is given as
             *       1.  Since there are no times of evaluation in the middle
             *       of the step, a multi-step method is similar to a one step
             *       method with one stage.
             */
            void preStep (RealType time, RealType dt, int stages)
            {
            }

            //! to be called once at the end of each time step
            /**
             * \note With the OneStepMethod and the ExplicitOneStepMetod, for
             *       reasons unknown this is only called for temporal but not
             *       for spatial local operators.  With the MultiStepMethod
             *       this *is* called for all local operators.
             */
            void postStep ()
            {
            }

            //! to be called once before each stage
            /**
             * \param time Time of the stage
             * \param r    Number of the stage, r ∈ [1, nstages] inclusive,
             *             where nstages is the number of stage in the step
             *             given in the previous call to preStep()
             *
             * \note For ExplicitOneStepMethod the time given here for stage 1
             *       may be incorrect, since the time step size is only
             *       finally determined after the first stage has been
             *       assembled.
             *
             * \note For the MultiStepMethod, this is called once after
             *       preStep() with r=1.
             */
            void preStage (RealType time, int r)
            {
                stage = r;
            }

            //! get current stage
            /**
             * \return The current stage number previously set by preStage().
             */
            int getStage () const
            {
                return stage;
            }

            //! to be called once at the end of each stage
            void postStage ()
            {
            }

            //! to be called after stage 1
            /**
             * \note Only used by the ExplicitOneStepMethod.
             *
             * This may be called on the spatial local operator in the case of
             * an explicit one step scheme.  It is called after stage 1 has
             * been assembled (so the time given to preStep() may not apply
             * anymore in this case).  All the alpha_*() and lambda_*()
             * methods should have been called, so they are a good place to
             * generate the information returned here.
             */
            RealType suggestTimestep (RealType dt) const
            {
                return dt;
            }

        private:
            int stage;
            R t;
        };
        //! \} group LocalOperator
    }
}

#endif // DUNE_PDELAB_LOCALOPERATOR_IDEFAULT_HH
