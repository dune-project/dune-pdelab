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
        /**
         * \nosubgrouping
         */
        class LocalOperatorDefaultFlags
        {
        public:
            //! \name Flags for the sparsity pattern
            //! \{

            //! \brief Whether to assemble the pattern on the elements,
            //!        i.e. whether or not pattern_volume() should be called.
            enum { /*! \hideinitializer */ doPatternVolume = false };
            //! \brief Whether to assemble the pattern on the elements after
            //!        the skeleton has been handled, i.e. whether or not
            //!        pattern_volume_post_skeleton() should be called.
            enum { /*! \hideinitializer */ doPatternVolumePostSkeleton = false };
            //! \brief Whether to assemble the pattern on the interior
            //!        intersections, i.e. whether or not pattern_skeleton()
            //!        should be called.
            enum { /*! \hideinitializer */ doPatternSkeleton = false };
            //! \brief Whether to assemble the pattern on the boundary
            //!        intersections, i.e. whether or not pattern_boundary()
            //!        should be called.
            enum { /*! \hideinitializer */ doPatternBoundary = false };

            //! \} Flags for the sparsity pattern

            //! \name Flags for the non-constant part of the residual and the jacobian
            //! \{

            //! \brief Whether to call the local operator's alpha_volume(),
            //!        jacobian_apply_volume() and jacobian_volume().
            enum { /*! \hideinitializer */ doAlphaVolume = false };
            //! \brief Whether to call the local operator's
            //!        alpha_volume_post_skeleton(),
            //!        jacobian_apply_volume_post_skeleton() and
            //!        jacobian_volume_post_skeleton().
            enum { /*! \hideinitializer */ doAlphaVolumePostSkeleton = false };
            //! \brief Whether to call the local operator's alpha_skeleton(),
            //!        jacobian_apply_skeleton() and jacobian_skeleton().
            enum { /*! \hideinitializer */ doAlphaSkeleton = false };
            //! \brief Whether to call the local operator's alpha_boundary(),
            //!        jacobian_apply_boundary() and jacobian_boundary().
            enum { /*! \hideinitializer */ doAlphaBoundary = false };

            //! \} Flags for the non-constant part of the residual and the jacobian

            //! \name Flags for the constant part of the residual
            //! \{

            //! \brief Whether to call the local operator's lambda_volume().
            enum { /*! \hideinitializer */ doLambdaVolume = false };
            //! \brief Whether to call the local operator's
            //!        lambda_volume_post_skeleton().
            enum { /*! \hideinitializer */ doLambdaVolumePostSkeleton = false };
            //! \brief Whether to call the local operator's lambda_skeleton().
            enum { /*! \hideinitializer */ doLambdaSkeleton = false };
            //! \brief Whether to call the local operator's lambda_boundary().
            enum { /*! \hideinitializer */ doLambdaBoundary = false };

            //! \} Flags for the constant part of the residual

            //! \name Special flags
            //! \{

            //! \brief Whether to visit the skeleton methods from both sides
            enum { /*! \hideinitializer */ doSkeletonTwoSided = false };

            //! \brief Wheter the local operator describes a linear problem
            enum { /*! \hideinitializer */ isLinear = true };

            //! \} Special flags
        };


        //! Namespace with decorator classes that influence assembler behavior.
        namespace lop {

            //! Decorator base class for local operators that have a diagonal jacobian matrix.
            /**
             * By inheriting from this decorator class, local operators assert that
             * their jacobian is completely diagonal.
             * This information can be used by the assembler to e.g. switch to a diagonal
             * local matrix that gets passed to the LocalOperator, which can save a lot of
             * memory bandwidth for large local function spaces.
             */
            struct DiagonalJacobian
            {};

        }

        //! \} group LocalOperator
    }
}

#endif // DUNE_PDELAB_LOCALOPERATOR_FLAGS_HH
