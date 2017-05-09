// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_GRIDOPERATOR_COMMON_LOCALASSEMBLERENGINEBASE_HH
#define DUNE_PDELAB_GRIDOPERATOR_COMMON_LOCALASSEMBLERENGINEBASE_HH

namespace Dune {
  namespace PDELab {

      /** \addtogroup GridOperator
       *  \{
       */

      //! Base class for LocalAssemblerEngine implementations to avoid boilerplate code.
      /**
       * This class is a very simple LocalAssemblerEngine that causes the Assembler to
       * iterate over the whole grid, but will not do anything. In particular, all the
       * require...() methods return false, all assembly methods are empty and assembleCell()
       * simply returns false to make sure each cell gets processed by the Assembler.
       */
      class LocalAssemblerEngineBase
      {

      public:

        struct Traits
        {};

        //! @name Query methods - return false by default
        //! @{

        bool requireSkeleton() const
        {
          return false;
        }

        bool requireSkeletonTwoSided() const
        {
          return false;
        }

        bool requireUVVolume() const
        {
          return false;
        }

        bool requireVVolume() const
        {
          return false;
        }

        bool requireUVSkeleton() const
        {
          return false;
        }

        bool requireVSkeleton() const
        {
          return false;
        }

        bool requireUVBoundary() const
        {
          return false;
        }

        bool requireVBoundary() const
        {
          return false;
        }

        bool requireUVProcessor() const
        {
          return false;
        }

        bool requireVProcessor() const
        {
          return false;
        }

        bool requireUVEnrichedCoupling() const
        {
          return false;
        }

        bool requireVEnrichedCoupling() const
        {
          return false;
        }

        bool requireUVVolumePostSkeleton() const
        {
          return false;
        }

        bool requireVVolumePostSkeleton() const
        {
          return false;
        }

        auto partition() const
        {
          return Partitions::interiorBorder;
        }

        //! @}

        //! @name Callbacks for LocalFunctionSpace binding and unbinding events
        //! @{

        template<typename EG,
                 typename LFSU, typename LFSV>
        void onBindLFSUV(const EG& eg,
                         const LFSU& lfsu, const LFSV& lfsv)
        {
        }

        template<typename EG,
                 typename LFSV>
        void onBindLFSV(const EG& eg,
                        const LFSV& lfsv)
        {
        }

        template<typename EG,
                 typename LFSU, typename LFSV>
        void onUnbindLFSUV(const EG& eg,
                           const LFSU& lfsu, const LFSV& lfsv)
        {
        }

        template<typename EG,
                 typename LFSV_S>
        void onUnbindLFSV(const EG& eg,
                          const LFSV_S& lfsv_s)
        {
        }


        template<typename IG,
                 typename LFSU, typename LFSV>
        void onBindLFSUVInside(const IG& ig,
                               const LFSU& lfsu, const LFSV& lfsv)
        {
        }

        template<typename IG,
                 typename LFSV>
        void onBindLFSVInside(const IG& ig,
                              const LFSV& lfsv)
        {
        }

        template<typename IG,
                 typename LFSU, typename LFSV>
        void onUnbindLFSUVInside(const IG& ig,
                                 const LFSU& lfsu, const LFSV& lfsv)
        {
        }

        template<typename IG,
                 typename LFSV_S>
        void onUnbindLFSVInside(const IG& ig,
                                const LFSV_S& lfsv_s)
        {
        }


        template<typename IG,
                 typename LFSU_S, typename LFSV_S,
                 typename LFSU_N, typename LFSV_N>
        void onBindLFSUVOutside(const IG& ig,
                                const LFSU_S& lfsu_s, const LFSV_S& lfsv_s,
                                const LFSU_N& lfsu_n, const LFSV_N& lfsv_n)
        {
        }

        template<typename IG,
                 typename LFSV_S,
                 typename LFSV_N>
        void onBindLFSVOutside(const IG& ig,
                               const LFSV_S& lfsv_s,
                               const LFSV_N& lfsv_n)
        {
        }

        template<typename IG,
                 typename LFSU_S, typename LFSV_S,
                 typename LFSU_N, typename LFSV_N>
        void onUnbindLFSUVOutside(const IG& ig,
                                  const LFSU_S& lfsu_s, const LFSV_S& lfsv_s,
                                  const LFSU_N& lfsu_n, const LFSV_N& lfsv_n)
        {
        }

        template<typename IG,
                 typename LFSV_S,
                 typename LFSV_N>
        void onUnbindLFSVOutside(const IG& ig,
                                 const LFSV_S& lfsv_s,
                                 const LFSV_N& lfsv_n)
        {
        }


        template<typename IG,
                 typename LFSU_S, typename LFSV_S,
                 typename LFSU_N, typename LFSV_N,
                 typename LFSU_C, typename LFSV_C>
        void onBindLFSUVCoupling(const IG& ig,
                                 const LFSU_S& lfsu_s, const LFSV_S& lfsv_s,
                                 const LFSU_N& lfsu_n, const LFSV_N& lfsv_n,
                                 const LFSU_C& lfsu_c, const LFSV_C& lfsv_c)
        {
        }

        template<typename IG,
                 typename LFSV_S,
                 typename LFSV_N,
                 typename LFSV_C>
        void onBindLFSVCoupling(const IG& ig,
                                const LFSV_S& lfsv_s,
                                const LFSV_N& lfsv_n,
                                const LFSV_C& lfsv_c)
        {
        }

        template<typename IG,
                 typename LFSU_S, typename LFSV_S,
                 typename LFSU_N, typename LFSV_N,
                 typename LFSU_C, typename LFSV_C>
        void onUnbindLFSUVCoupling(const IG& ig,
                                   const LFSU_S& lfsu_s, const LFSV_S& lfsv_s,
                                   const LFSU_N& lfsu_n, const LFSV_N& lfsv_n,
                                   const LFSU_C& lfsu_c, const LFSV_C& lfsv_c)
        {
        }

        template<typename IG,
                 typename LFSV_S,
                 typename LFSV_N,
                 typename LFSV_C>
        void onUnbindLFSVCoupling(const IG& ig,
                                  const LFSV_S& lfsv_s,
                                  const LFSV_N& lfsv_n,
                                  const LFSV_C& lfsv_c)
        {
        }

        template<typename LFSU>
        void loadCoefficientsLFSUInside(const LFSU& lfsu_s)
        {
        }

        template<typename LFSU_N>
        void loadCoefficientsLFSUOutside(const LFSU_N& lfsu_n)
        {
        }

        template<typename LFSU_C>
        void loadCoefficientsLFSUCoupling(const LFSU_C& lfsu_c)
        {
        }

        //! @}

        //! @name Assembly methods
        //! @{

        //! Method for per-cell assembly setup and possibly aborting assembly of current cell
        //! - returns false by default to continue cell assembly.
        template<typename EG>
        bool assembleCell(const EG & eg)
        {
          return false;
        }

        template<typename EG, typename LFSU, typename LFSV>
        void assembleUVVolume(const EG& eg, const LFSU& lfsu, const LFSV& lfsv)
        {
        }

        template<typename EG, typename LFSV>
        void assembleVVolume(const EG& eg, const LFSV& lfsv)
        {
        }


        template<typename IG, typename LFSU_S, typename LFSV_S, typename LFSU_N, typename LFSV_N>
        void assembleUVSkeleton(const IG& ig,
                                const LFSU_S& lfsu_s, const LFSV_S& lfsv_s,
                                const LFSU_N& lfsu_n, const LFSV_N& lfsv_n)
        {
        }

        template<typename IG, typename LFSV_S, typename LFSV_N>
        void assembleVSkeleton(const IG& ig,
                               const LFSV_S& lfsv_s,
                               const LFSV_N& lfsv_n)
        {
        }


        template<typename IG, typename LFSU, typename LFSV>
        void assembleUVBoundary(const IG& ig, const LFSU& lfsu, const LFSV& lfsv)
        {
        }

        template<typename IG, typename LFSV>
        void assembleVBoundary(const IG& ig, const LFSV& lfsv)
        {
        }


        template<typename IG, typename LFSU, typename LFSV>
        void assembleUVProcessor(const IG& ig, const LFSU& lfsu, const LFSV& lfsv)
        {
        }

        template<typename IG, typename LFSV>
        void assembleVProcessor(const IG& ig, const LFSV& lfsv)
        {
        }


        template<typename IG,
                 typename LFSU_S, typename LFSV_S,
                 typename LFSU_N, typename LFSV_N,
                 typename LFSU_C, typename LFSV_C>
        void assembleUVEnrichedCoupling(const IG& ig,
                                        const LFSU_S& lfsu_s, const LFSV_S& lfsv_s,
                                        const LFSU_N& lfsu_n, const LFSV_N& lfsv_n,
                                        const LFSU_C& lfsu_c, const LFSV_C& lfsv_c)
        {
        }

        template<typename IG,
                 typename LFSV_S,
                 typename LFSV_N,
                 typename LFSV_C>
        void assembleVEnrichedCoupling(const IG& ig,
                                       const LFSV_S& lfsv_s,
                                       const LFSV_N& lfsv_n,
                                       const LFSV_C& lfsv_c)
        {
        }


        template<typename EG, typename LFSU, typename LFSV>
        void assembleUVVolumePostSkeleton(const EG& eg, const LFSU& lfsu, const LFSV& lfsv)
        {
        }

        template<typename EG, typename LFSV>
        void assembleVVolumePostSkeleton(const EG& eg, const LFSV& lfsv)
        {
        }

        //! @}

        //! @name Global assembly preparation and finalization methods
        //! @{

        void preAssembly()
        {
        }

        template<typename GFSU, typename GFSV>
        void postAssembly(const GFSU& gfsu, const GFSV& gfsv)
        {
        }

        //! @}

      };


      //! \} group GridOperator

  } // namespace PDELab
} //namespace Dune

#endif // DUNE_PDELAB_GRIDOPERATOR_COMMON_LOCALASSEMBLERENGINEBASE_HH
