#ifndef DUNE_PDELAB_DEFAULT_ASSEMBLER_HH
#define DUNE_PDELAB_DEFAULT_ASSEMBLER_HH

#include<dune/pdelab/gridoperator/common/assemblerutilities.hh>

namespace Dune{
  namespace PDELab{

    /**
       \brief The assembler for standard DUNE grid

       * \tparam GFSU GridFunctionSpace for ansatz functions
       * \tparam GFSV GridFunctionSpace for test functions
       * \tparam nonoverlapping_mode Indicates whether assembling is done for overlap cells
       */

    template<typename GFSU, typename GFSV, bool nonoverlapping_mode=false>
    class DefaultAssembler {
    public:

      //! Types related to current grid view
      //! @{
      typedef typename GFSU::Traits::GridViewType GV;
      typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
      typedef typename GV::Traits::template Codim<0>::Entity Element;
      typedef typename GV::IntersectionIterator IntersectionIterator;
      typedef typename IntersectionIterator::Intersection Intersection;
      //! @}

      //! Grid function spaces
      //! @{
      typedef GFSU TrialGridFunctionSpace;
      typedef GFSV TestGridFunctionSpace;
      //! @}

      //! Size type as used in grid function space
      typedef typename GFSU::Traits::SizeType SizeType;

      //! Static check on whether this is a Galerkin method
      static const bool isGalerkinMethod = Dune::is_same<GFSU,GFSV>::value;

      DefaultAssembler (const GFSU& gfsu_, const GFSV& gfsv_) 
        : gfsu(gfsu_), gfsv(gfsv_), lfsu(gfsu_), lfsv(gfsv_), 
          lfsun(gfsu_), lfsvn(gfsv_)
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

      template<class LocalAssemblerEngine>
      void assemble(LocalAssemblerEngine & assembler_engine) const
      {
        typedef typename GV::Traits::template Codim<0>::Entity Element;

        // Notify assembler engine about oncoming assembly
        assembler_engine.preAssembly();

        // Map each cell to unique id
        Dune::PDELab::MultiGeomUniqueIDMapper<GV> cell_mapper(gfsu.gridview());

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

        // Traverse grid view
        for (ElementIterator it = gfsu.gridview().template begin<0>();
             it!=gfsu.gridview().template end<0>(); ++it)
          {
            // Compute unique id
            const typename GV::IndexSet::IndexType ids = cell_mapper.map(*it);

            ElementGeometry<Element> eg(*it);

            if(assembler_engine.assembleCell(eg))
              continue;

            // Bind local test function space to element
            lfsv.bind( *it );

            // Notify assembler engine about bind
            assembler_engine.onBindLFSV(eg,lfsv);

            // Volume integration
            assembler_engine.assembleVVolume(eg,lfsv);

            // Bind local trial function space to element
            lfsu.bind( *it );

            // Notify assembler engine about bind
            assembler_engine.onBindLFSUV(eg,lfsu,lfsv);

            // Load coefficients of local functions
            assembler_engine.loadCoefficientsLFSUInside(lfsu);

            // Volume integration
            assembler_engine.assembleUVVolume(eg,lfsu,lfsv);

            // Skip if no intersection iterator is needed
            if (require_uv_skeleton || require_v_skeleton ||
                require_uv_boundary || require_v_boundary ||
                require_uv_processor || require_v_processor)
              {
                // Traverse intersections
                unsigned int intersection_index = 0;
                IntersectionIterator endit = gfsu.gridview().iend(*it);
                IntersectionIterator iit = gfsu.gridview().ibegin(*it);
                for(; iit!=endit; ++iit, ++intersection_index)
                  {

                    IntersectionGeometry<Intersection> ig(*iit,intersection_index);

                    switch (IntersectionType::get(*iit))
                      {
                      case IntersectionType::skeleton:
                        // the specific ordering of the if-statements in the old code caused periodic
                        // boundary intersection to be handled the same as skeleton intersections
                      case IntersectionType::periodic:
                        if (require_uv_skeleton || require_v_skeleton)
                          {
                            // compute unique id for neighbor

                            const typename GV::IndexSet::IndexType idn = cell_mapper.map(*(iit->outside()));

                            // Visit face if id is bigger
                            bool visit_face = ids > idn || require_skeleton_two_sided;
                            // or interior is a ghost
                            visit_face |= (nonoverlapping_mode &&
                                           (iit->inside())->partitionType()!=Dune::InteriorEntity);

                            // unique vist of intersection
                            if (visit_face)
                              {
                                // Bind local test space to neighbor element
                                lfsvn.bind(*(iit->outside()));

                                // Notify assembler engine about binds
                                assembler_engine.onBindLFSVOutside(ig,lfsvn);

                                // Skeleton integration
                                assembler_engine.assembleVSkeleton(ig,lfsv,lfsvn);

                                if(require_uv_skeleton){

                                  // Bind local trial space to neighbor element
                                  lfsun.bind(*(iit->outside()));

                                  // Notify assembler engine about binds
                                  assembler_engine.onBindLFSUVOutside(ig,lfsun,lfsvn);

                                  // Load coefficients of local functions
                                  assembler_engine.loadCoefficientsLFSUOutside(lfsun);

                                  // Skeleton integration
                                  assembler_engine.assembleUVSkeleton(ig,lfsu,lfsv,lfsun,lfsvn);

                                  // Notify assembler engine about unbinds
                                  assembler_engine.onUnbindLFSUVOutside(ig,lfsun,lfsvn);
                                }

                                // Notify assembler engine about unbinds
                                assembler_engine.onUnbindLFSVOutside(ig,lfsvn);
                              }
                          }
                        break;

                      case IntersectionType::boundary:
                        if(require_uv_boundary || require_v_boundary )
                          {

                            // Boundary integration
                            assembler_engine.assembleVBoundary(ig,lfsv);

                            if(require_uv_boundary){
                              // Boundary integration
                              assembler_engine.assembleUVBoundary(ig,lfsu,lfsv);
                            }
                          }
                        break;

                      case IntersectionType::processor:
                        if(require_uv_processor || require_v_processor )
                          {

                            // Processor integration
                            assembler_engine.assembleVProcessor(ig,lfsv);

                            if(require_uv_processor){
                              // Processor integration
                              assembler_engine.assembleUVProcessor(ig,lfsu,lfsv);
                            }
                          }
                        break;
                      } // switch

                  } // iit
              } // do skeleton
            
            if(require_uv_post_skeleton || require_v_post_skeleton){
              // Volume integration
              assembler_engine.assembleVVolumePostSkeleton(eg,lfsv);

              if(require_uv_post_skeleton){
                // Volume integration
                assembler_engine.assembleUVVolumePostSkeleton(eg,lfsu,lfsv);
              }
            }

            // Notify assembler engine about unbinds
            assembler_engine.onUnbindLFSUV(eg,lfsu,lfsv);

            // Notify assembler engine about unbinds
            assembler_engine.onUnbindLFSV(eg,lfsv);

          } // it

        // Notify assembler engine that assembly is finished
        assembler_engine.postAssembly();

      }

    private:

      /* global function spaces */
      const GFSU& gfsu;
      const GFSV& gfsv;

      /* local function spaces */
      typedef Dune::PDELab::LocalFunctionSpace<GFSU, Dune::PDELab::TrialSpaceTag> LFSU;
      typedef Dune::PDELab::LocalFunctionSpace<GFSV, Dune::PDELab::TestSpaceTag> LFSV;
      // local function spaces in local cell
      mutable LFSU lfsu;
      mutable LFSV lfsv;
      // local function spaces in neighbor
      mutable LFSU lfsun;
      mutable LFSV lfsvn;
    };

  };
};
#endif
