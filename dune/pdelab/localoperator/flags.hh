#ifndef DUNE_PDELAB_LOCALOPERATOR_FLAGS_HH
#define DUNE_PDELAB_LOCALOPERATOR_FLAGS_HH

namespace Dune
{
    namespace PDELab
    {
        class LocalOperatorDefaultFlags
        {
        public:
            // pattern assembly flags
            enum { doPatternVolume = false };
            enum { doPatternSkeleton = false };

            // residual assembly flags
            enum { doAlphaVolume = false };
            enum { doAlphaVolumePostSkeleton = false };
            enum { doAlphaSkeleton = false };
            enum { doAlphaBoundary = false };
            enum { doLambdaVolume = false };
            enum { doLambdaVolumePostSkeleton = false };
            enum { doLambdaSkeleton = false };
            enum { doLambdaBoundary = false };
            enum { doSkeletonTwoSided = false };
        };
    }
}

#endif
