// -*- tab-width: 2; indent-tabs-mode: nil -*-
// vi: set et ts=2 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDOPERATOR_FASTDG_ASSEMBLER_HH
#define DUNE_PDELAB_GRIDOPERATOR_FASTDG_ASSEMBLER_HH

#include <dune/common/typetraits.hh>
#include <dune/pdelab/gridoperator/common/assemblerutilities.hh>
#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>
#include <dune/pdelab/common/geometrywrapper.hh>
#include <dune/pdelab/common/intersectiontype.hh>

namespace Dune{
  namespace PDELab{

    /**
       \brief The fast DG assembler for standard DUNE grid.

       * \tparam GFSU GridFunctionSpace for ansatz functions
       * \tparam GFSV GridFunctionSpace for test functions
       * \tparam nonoverlapping_mode Indicates whether assembling is done for overlap cells
       */

    template<typename GFSU, typename GFSV, typename CU, typename CV, bool nonoverlapping_mode=false>
    class FastDGAssembler {
    public:

      //! Types related to current grid view
      //! @{
      using EntitySet = typename GFSU::Traits::EntitySet;
      using Element = typename EntitySet::Element;
      using Intersection = typename EntitySet::Intersection;
      //! @}

      //! Grid function spaces
      //! @{
      typedef GFSU TrialGridFunctionSpace;
      typedef GFSV TestGridFunctionSpace;
      //! @}

      //! Size type as used in grid function space
      typedef typename GFSU::Traits::SizeType SizeType;

      //! Static check on whether this is a Galerkin method
      static const bool isGalerkinMethod = std::is_same<GFSU,GFSV>::value;

      FastDGAssembler (const GFSU& gfsu_, const GFSV& gfsv_, const CU& cu_, const CV& cv_)
        : gfsu(gfsu_)
        , gfsv(gfsv_)
        , cu(cu_)
        , cv(cv_)
        , lfsu(gfsu_)
        , lfsv(gfsv_)
        , lfsun(gfsu_)
        , lfsvn(gfsv_)
      { }

      FastDGAssembler (const GFSU& gfsu_, const GFSV& gfsv_)
        : gfsu(gfsu_)
        , gfsv(gfsv_)
        , cu()
        , cv()
        , lfsu(gfsu_)
        , lfsv(gfsv_)
        , lfsun(gfsu_)
        , lfsvn(gfsv_)
      { }

      //! Get the trial grid function space
      const GFSU& trialGridFunctionSpace() const
      {
        return gfsu;
      }

      //! Get the test grid function space
      const GFSV& testGridFunctionSpace() const
      {
        return gfsv;
      }

      // Assembler (const GFSU& gfsu_, const GFSV& gfsv_)
      //   : gfsu(gfsu_), gfsv(gfsv_), lfsu(gfsu_), lfsv(gfsv_),
      //     lfsun(gfsu_), lfsvn(gfsv_),
      //     sub_triangulation(ST(gfsu_.gridview(),Dune::PDELab::NoSubTriangulationImp()))
      // { }

      //! do the assembly
      /**
       * \param engineFactory Factory object used to get the engine.
       * \param la            Local assembler to get the engine from.
       */
      template<class EngineFactory, class LocalAssembler>
      void assemble(const EngineFactory &engineFactory,
                    LocalAssembler &la) const
      {
        assemble(engineFactory(la));
      }

      template<class LocalAssemblerEngine>
      void assemble(LocalAssemblerEngine & assembler_engine) const
      {
        const bool fast = true;
        typedef LFSIndexCache<LFSU,CU,fast> LFSUCache;

        typedef LFSIndexCache<LFSV,CV,fast> LFSVCache;

        const bool needs_constraints_caching = assembler_engine.needsConstraintsCaching(cu,cv);

        LFSUCache lfsu_cache(lfsu,cu,needs_constraints_caching);
        LFSVCache lfsv_cache(lfsv,cv,needs_constraints_caching);
        LFSUCache lfsun_cache(lfsun,cu,needs_constraints_caching);
        LFSVCache lfsvn_cache(lfsvn,cv,needs_constraints_caching);

        // Notify assembler engine about oncoming assembly
        assembler_engine.preAssembly();

        // Extract integration requirements from the local assembler
        const bool require_uv_skeleton = assembler_engine.requireUVSkeleton();
        const bool require_v_skeleton = assembler_engine.requireVSkeleton();
        const bool require_uv_boundary = assembler_engine.requireUVBoundary();
        const bool require_v_boundary = assembler_engine.requireVBoundary();
        const bool require_uv_processor = assembler_engine.requireUVBoundary();
        const bool require_v_processor = assembler_engine.requireVBoundary();
        const bool require_uv_post_skeleton = assembler_engine.requireUVVolumePostSkeleton();
        const bool require_v_post_skeleton = assembler_engine.requireVVolumePostSkeleton();
        const bool require_skeleton_two_sided = assembler_engine.requireSkeletonTwoSided();

        auto entity_set = gfsu.entitySet();
        auto& index_set = entity_set.indexSet();

        // Traverse grid view
        for (const auto& element : elements(entity_set))
          {
            // Compute unique id
            auto ids = index_set.uniqueIndex(element);

            ElementGeometry<Element> eg(element);

            if(assembler_engine.assembleCell(eg))
              continue;

            // Bind local test function space to element
            lfsv.bind(element,std::integral_constant<bool,fast>{});
            lfsv_cache.update();

            // Notify assembler engine about bind
            assembler_engine.onBindLFSV(eg,lfsv_cache);

            // Volume integration
            assembler_engine.assembleVVolume(eg,lfsv_cache);

            // Bind local trial function space to element
            lfsu.bind(element,std::integral_constant<bool,fast>{});
            lfsu_cache.update();

            // Notify assembler engine about bind
            assembler_engine.onBindLFSUV(eg,lfsu_cache,lfsv_cache);

            // Load coefficients of local functions
            assembler_engine.loadCoefficientsLFSUInside(lfsu_cache);

            // Volume integration
            assembler_engine.assembleUVVolume(eg,lfsu_cache,lfsv_cache);

            // Skip if no intersection iterator is needed
            if (require_uv_skeleton || require_v_skeleton ||
                require_uv_boundary || require_v_boundary ||
                require_uv_processor || require_v_processor)
              {
                // Traverse intersections
                unsigned int intersection_index = 0;
                for(const auto& intersection : intersections(entity_set,element))
                  {

                    IntersectionGeometry<Intersection> ig(intersection,intersection_index);

                    auto intersection_data = classifyIntersection(entity_set,intersection);
                    auto intersection_type = std::get<0>(intersection_data);
                    auto& outside_element = std::get<1>(intersection_data);

                    switch (intersection_type)
                      {
                      case IntersectionType::skeleton:
                        // the specific ordering of the if-statements in the old code caused periodic
                        // boundary intersection to be handled the same as skeleton intersections
                      case IntersectionType::periodic:
                        if (require_uv_skeleton || require_v_skeleton)
                          {
                            // compute unique id for neighbor

                            auto idn = index_set.uniqueIndex(outside_element);

                            // Visit face if id is bigger
                            bool visit_face = ids > idn || require_skeleton_two_sided;

                            // unique vist of intersection
                            if (visit_face)
                              {
                                // Bind local test space to neighbor element
                                lfsvn.bind(outside_element,std::integral_constant<bool,fast>{});
                                lfsvn_cache.update();

                                // Notify assembler engine about binds
                                assembler_engine.onBindLFSVOutside(ig,lfsv_cache,lfsvn_cache);

                                // Skeleton integration
                                assembler_engine.assembleVSkeleton(ig,lfsv_cache,lfsvn_cache);

                                if(require_uv_skeleton){

                                  // Bind local trial space to neighbor element
                                  lfsun.bind(outside_element,std::integral_constant<bool,fast>{});
                                  lfsun_cache.update();

                                  // Notify assembler engine about binds
                                  assembler_engine.onBindLFSUVOutside(ig,
                                                                      lfsu_cache,lfsv_cache,
                                                                      lfsun_cache,lfsvn_cache);

                                  // Load coefficients of local functions
                                  assembler_engine.loadCoefficientsLFSUOutside(lfsun_cache);

                                  // Skeleton integration
                                  assembler_engine.assembleUVSkeleton(ig,lfsu_cache,lfsv_cache,lfsun_cache,lfsvn_cache);

                                  // Notify assembler engine about unbinds
                                  assembler_engine.onUnbindLFSUVOutside(ig,
                                                                        lfsu_cache,lfsv_cache,
                                                                        lfsun_cache,lfsvn_cache);
                                }

                                // Notify assembler engine about unbinds
                                assembler_engine.onUnbindLFSVOutside(ig,lfsv_cache,lfsvn_cache);
                              }
                          }
                        break;

                      case IntersectionType::boundary:
                        if(require_uv_boundary || require_v_boundary )
                          {

                            // Boundary integration
                            assembler_engine.assembleVBoundary(ig,lfsv_cache);

                            if(require_uv_boundary){
                              // Boundary integration
                              assembler_engine.assembleUVBoundary(ig,lfsu_cache,lfsv_cache);
                            }
                          }
                        break;

                      case IntersectionType::processor:
                        if(require_uv_processor || require_v_processor )
                          {

                            // Processor integration
                            assembler_engine.assembleVProcessor(ig,lfsv_cache);

                            if(require_uv_processor){
                              // Processor integration
                              assembler_engine.assembleUVProcessor(ig,lfsu_cache,lfsv_cache);
                            }
                          }
                        break;
                      } // switch

                    ++intersection_index;
                  } // iit
              } // do skeleton

            if(require_uv_post_skeleton || require_v_post_skeleton){
              // Volume integration
              assembler_engine.assembleVVolumePostSkeleton(eg,lfsv_cache);

              if(require_uv_post_skeleton){
                // Volume integration
                assembler_engine.assembleUVVolumePostSkeleton(eg,lfsu_cache,lfsv_cache);
              }
            }

            // Notify assembler engine about unbinds
            assembler_engine.onUnbindLFSUV(eg,lfsu_cache,lfsv_cache);

            // Notify assembler engine about unbinds
            assembler_engine.onUnbindLFSV(eg,lfsv_cache);

          } // it

        // Notify assembler engine that assembly is finished
        assembler_engine.postAssembly(gfsu,gfsv);

      }

    private:

      /* global function spaces */
      const GFSU& gfsu;
      const GFSV& gfsv;

      typename std::conditional<
        std::is_same<CU,EmptyTransformation>::value,
        const CU,
        const CU&
        >::type cu;
      typename std::conditional<
        std::is_same<CV,EmptyTransformation>::value,
        const CV,
        const CV&
        >::type cv;

      /* local function spaces */
      typedef LocalFunctionSpace<GFSU, TrialSpaceTag> LFSU;
      typedef LocalFunctionSpace<GFSV, TestSpaceTag> LFSV;
      // local function spaces in local cell
      mutable LFSU lfsu;
      mutable LFSV lfsv;
      // local function spaces in neighbor
      mutable LFSU lfsun;
      mutable LFSV lfsvn;

    }; // end class FastDGAssembler
  } // end namespace PDELab
} // end namespace Dune
#endif
