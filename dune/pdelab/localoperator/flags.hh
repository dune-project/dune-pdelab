// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=8 sw=4 sts=4:
#ifndef DUNE_PDELAB_LOCALOPERATOR_FLAGS_HH
#define DUNE_PDELAB_LOCALOPERATOR_FLAGS_HH

namespace Dune
{
    namespace PDELab
    {
        //! \addtogroup LocalOperator
        //! \ingroup PDELab
        //! \{

        //! Default flags for all local operators
        class LocalOperatorDefaultFlags
        {
        public:
            // pattern assembly flags
            //! \todo
            enum { /*! \hideinitializer */ doPatternVolume = false };
            //! \todo
            enum { /*! \hideinitializer */ doPatternSkeleton = false };

            // residual assembly flags
            //! \brief Whether not to call the local operators alpha_volume()
            //!        and jacobian_volume()
            enum { /*! \hideinitializer */ doAlphaVolume = false };
            //! \todo
            enum { /*! \hideinitializer */ doAlphaVolumePostSkeleton = false };
            //! \todo
            enum { /*! \hideinitializer */ doAlphaSkeleton = false };
            //! \todo
            enum { /*! \hideinitializer */ doAlphaBoundary = false };
            //! \todo
            enum { /*! \hideinitializer */ doLambdaVolume = false };
            //! \todo
            enum { /*! \hideinitializer */ doLambdaVolumePostSkeleton = false };
            //! \todo
            enum { /*! \hideinitializer */ doLambdaSkeleton = false };
            //! \todo
            enum { /*! \hideinitializer */ doLambdaBoundary = false };
            //! \todo
            enum { /*! \hideinitializer */ doSkeletonTwoSided = false };
        };
        //! \} group LocalOperator
    }
}

#endif
