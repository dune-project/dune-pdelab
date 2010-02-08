// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=8 sw=4 sts=4:
#ifndef DUNE_PDELAB_IDEFAULT_HH
#define DUNE_PDELAB_IDEFAULT_HH

namespace Dune
{
    namespace PDELab
    {
        //! \addtogroup LocalOperator
        //! \ingroup PDELab
        //! \{

        //! Default class for additional methods in instationary local operators
        template<class R>
        class InstationaryLocalOperatorDefaultMethods
        {
        public:
            //! set stage for subsequent evaluation
            void set_stage (int r)
            {
                stage = r;
            }

            //! get current stage
            int get_stage () const
            {
                return stage;
            }

            //! set time for subsequent evaluation
            void set_time (R t_)
            {
                t = t_;
            }

            //! get current time
            R get_time () const
            {
                return t;
            }

            //! to be called once before each time step
            void pre_step ()
            {
            }

            //! to be called once at the end of each time step
            void post_step ()
            {
            }

            //! to be called once before each stage
            void pre_stage ()
            {
            }

            //! to be called once at the end of each stage
            void post_stage ()
            {
            }

        private:
            int stage;
            R t;
        };
        //! \} group LocalOperator
    }
}

#endif
