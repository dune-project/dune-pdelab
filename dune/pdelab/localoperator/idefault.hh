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
        template<class R = double>
        class InstationaryLocalOperatorDefaultMethods
        {
        public:
            typedef R RealType;

            //! set time for subsequent evaluation
            void setTime (R t_)
            {
                t = t_;
            }

            //! get current time
            R getTime () const
            {
                return t;
            }

            //! to be called once before each time step
            void preStep (RealType time, int stages)
            {
            }

            //! to be called once at the end of each time step
            void postStep ()
            {
            }

            //! to be called once before each stage
            void preStage (RealType time, int r)
            {
                stage = r;
            }

            //! get current stage
            int getStage () const
            {
                return stage;
            }

            //! to be called once at the end of each stage
            void postStage ()
            {
            }

            //! to be called once before each stage
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

#endif
