#ifndef DUNE_PDELAB_TBB_ASSEMBLER_HH
#define DUNE_PDELAB_TBB_ASSEMBLER_HH

#include <cstddef>
#include <iostream>
#include <ostream>
#include <vector>
#include <type_traits>

#include <tbb/blocked_range.h>
#include <tbb/parallel_reduce.h>
#include <tbb/tbb_stddef.h>

#include <dune/common/classname.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/common/typetraits.hh>

#include <dune/grid/utility/intersectionrange.hh>

#include <dune/pdelab/gridoperator/common/assemblerutilities.hh>
#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>
#include <dune/pdelab/common/elementmapper.hh>
#include <dune/pdelab/common/geometrywrapper.hh>

namespace Dune{
  namespace PDELab{

    template<class T, bool = std::is_constructible<T, T&, tbb::split>::value>
    struct TBBAssemblerSplit
    {
      static constexpr bool value = false;
      static shared_ptr<T> split(const shared_ptr<T> &other)
      {
        DUNE_THROW(NotImplemented, className<T>() << " does not support "
                   "splitting");
      }
      static void join(T &self, T &other)
      {
        DUNE_THROW(NotImplemented, className<T>() << " does not support "
                   "splitting");
      }
    };
    template<class T>
    struct TBBAssemblerSplit<T, true>
    {
      static constexpr bool value = true;
      static shared_ptr<T> split(const shared_ptr<T> &other)
      {
        return make_shared<T>(*other, tbb::split());
      }
      static void join(T &self, T &other)
      {
        self.join(other);
      }
    };

    /**
       \brief The assembler for standard DUNE grid

       * \tparam Partitioning Partitioning to use for parallel traversal.
       * \tparam GFSU GridFunctionSpace for ansatz functions
       * \tparam GFSV GridFunctionSpace for test functions
       * \tparam nonoverlapping_mode Indicates whether assembling is done for overlap cells
       */

    template<typename Partitioning, typename GFSU, typename GFSV, typename CU,
             typename CV, bool nonoverlapping_mode=false>
    class TBBAssembler {
    protected:
      template<class EngineFactory, class LocalAssembler>
      class AssembleBody;

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

      TBBAssembler (const GFSU& gfsu_, const GFSV& gfsv_, const CU& cu_, const CV& cv_)
        : gfsu(gfsu_)
        , gfsv(gfsv_)
        , cu(cu_)
        , cv(cv_)
        , verbosity(0)
      { }

      TBBAssembler (const GFSU& gfsu_, const GFSV& gfsv_)
        : gfsu(gfsu_)
        , gfsv(gfsv_)
        , cu()
        , cv()
        , verbosity(0)
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

      void setPartitioning(const shared_ptr<const Partitioning> &partitioning_)
      {
        partitioning = partitioning_;
      }

      template<class Coloring>
      void setColoring(const shared_ptr<Coloring> &coloring)
      { }

      void setVerbosity(unsigned verbosity_)
      {
        verbosity = verbosity_;
      }

      //! do the assembly
      /**
       * \param engineFactory Factory object used to get the engine.
       * \param la            Local assembler to get the engine from.
       */
      template<class EngineFactory, class LocalAssembler>
      void assemble(const EngineFactory &engineFactory,
                    LocalAssembler &la) const
      {
        typedef typename GV::Traits::template Codim<0>::Entity Element;

        auto &assembler_engine = engineFactory(la);

        // Notify assembler engine about oncoming assembly
        assembler_engine.preAssembly();

        // Map each cell to unique id
        ElementMapper<GV> cell_mapper(gfsu.gridView());

        AssembleBody<EngineFactory, LocalAssembler>
          body(*this, engineFactory, la, assembler_engine, cell_mapper);

        tbb::blocked_range<std::size_t> range(0, partitioning->partitions());
        if(TBBAssemblerSplit<LocalAssembler>::value &&
           assembler_engine.threadSafe())
        {
          if(verbosity > 0 && gfsu.gridView().comm().rank() == 0)
            std::cout << "TBBAssembler: parallel iteration" << std::endl;
          tbb::parallel_reduce(range, body);
        }
        else
        {
          if(verbosity > 0 && gfsu.gridView().comm().rank() == 0)
            std::cout << "TBBAssembler: sequential iteration" << std::endl;
          body(range);
        }

        // Notify assembler engine that assembly is finished
        assembler_engine.postAssembly(gfsu,gfsv);

      }

    protected:
      shared_ptr<const Partitioning> partitioning;

      /* global function spaces */
      const GFSU& gfsu;
      const GFSV& gfsv;

      typename conditional<
        is_same<CU,EmptyTransformation>::value,
        const CU,
        const CU&
        >::type cu;
      typename conditional<
        is_same<CV,EmptyTransformation>::value,
        const CV,
        const CV&
        >::type cv;

      unsigned verbosity;
    };

    template<typename Partitioning, typename GFSU, typename GFSV, typename CU,
             typename CV, bool nonoverlapping_mode>
    template<class EngineFactory, class LocalAssembler>
    class TBBAssembler<Partitioning, GFSU, GFSV, CU, CV,
                       nonoverlapping_mode>::AssembleBody
    {
    protected:
      typedef TBBAssembler Assembler;

    private:
      /* local function spaces */
      typedef LocalFunctionSpace<GFSU, TrialSpaceTag> LFSU;
      typedef LocalFunctionSpace<GFSV, TestSpaceTag> LFSV;

      typedef LFSIndexCache<LFSU,CU> LFSUCache;

      typedef LFSIndexCache<LFSV,CV> LFSVCache;

      const Assembler &assembler;

      EngineFactory engineFactory;
      std::shared_ptr<LocalAssembler> localAssembler;
    protected:
      typedef decltype(engineFactory(*localAssembler)) LocalAssemblerEngine;
    private:
      LocalAssemblerEngine &engine;

      // Map each cell to unique id
      const ElementMapper<GV> &cell_mapper;

      // local function spaces in local cell
      LFSU lfsu;
      LFSV lfsv;
      // local function spaces in neighbor
      LFSU lfsun;
      LFSV lfsvn;

      LFSUCache lfsu_cache;
      LFSVCache lfsv_cache;
      LFSUCache lfsun_cache;
      LFSVCache lfsvn_cache;

      // Extract integration requirements from the local assembler
      const bool require_uv_skeleton;
      const bool require_v_skeleton;
      const bool require_uv_boundary;
      const bool require_v_boundary;
      const bool require_uv_processor;
      const bool require_v_processor;
      const bool require_uv_post_skeleton;
      const bool require_v_post_skeleton;
      const bool require_skeleton_two_sided;

    public:
      //! construct assembly body
      /**
       * This constructor is used to construct the initial body.  \c
       * preAssemble() is not automatically called on the engine passed in
       * here; it should be called externally before the call to operator().
       * Also, \c postAssemble() should be called on the engine after the call
       * to operator().
       */
      AssembleBody(const Assembler &assembler_,
                   const EngineFactory &engineFactory_,
                   LocalAssembler &localAssembler_,
                   LocalAssemblerEngine &engine_,
                   const ElementMapper<GV> &cell_mapper_) :
        assembler(assembler_),
        engineFactory(engineFactory_),
        localAssembler(stackobject_to_shared_ptr(localAssembler_)),
        engine(engine_),
        cell_mapper(cell_mapper_),
        lfsu(assembler.gfsu),
        lfsv(assembler.gfsv),
        lfsun(assembler.gfsu),
        lfsvn(assembler.gfsv),
        lfsu_cache(lfsu, assembler.cu),
        lfsv_cache(lfsv, assembler.cv),
        lfsun_cache(lfsun, assembler.cu),
        lfsvn_cache(lfsvn, assembler.cv),
        // Extract integration requirements from the local assembler
        require_uv_skeleton(engine.requireUVSkeleton()),
        require_v_skeleton(engine.requireVSkeleton()),
        require_uv_boundary(engine.requireUVBoundary()),
        require_v_boundary(engine.requireVBoundary()),
        require_uv_processor(engine.requireUVBoundary()),
        require_v_processor(engine.requireVBoundary()),
        require_uv_post_skeleton(engine.requireUVVolumePostSkeleton()),
        require_v_post_skeleton(engine.requireVVolumePostSkeleton()),
        require_skeleton_two_sided(engine.requireSkeletonTwoSided())
      { }

      //! splitting constructor
      /**
       * split another body from a template.
       *
       * The local assemblers are split too; if they don't have a splitting
       * constructor an exception is thrown.  An assembler engine is obtained
       * using the engine factory from the newly constructed local assembler.
       * \c other.engine.split(this->engine) is called to notify the new
       * engine that it is for a worker thread.
       */
      AssembleBody(AssembleBody &other, tbb::split) :
        assembler(other.assembler),
        engineFactory(other.engineFactory),
        localAssembler(TBBAssemblerSplit<LocalAssembler>
                       ::split(other.localAssembler)),
        engine(engineFactory(*localAssembler)),
        cell_mapper(other.cell_mapper),
        lfsu(assembler.gfsu),
        lfsv(assembler.gfsv),
        lfsun(assembler.gfsu),
        lfsvn(assembler.gfsv),
        lfsu_cache(lfsu, assembler.cu),
        lfsv_cache(lfsv, assembler.cv),
        lfsun_cache(lfsun, assembler.cu),
        lfsvn_cache(lfsvn, assembler.cv),
        // Extract integration requirements
        require_uv_skeleton(other.require_uv_skeleton),
        require_v_skeleton(other.require_v_skeleton),
        require_uv_boundary(other.require_uv_boundary),
        require_v_boundary(other.require_v_boundary),
        require_uv_processor(other.require_uv_processor),
        require_v_processor(other.require_v_processor),
        require_uv_post_skeleton(other.require_uv_post_skeleton),
        require_v_post_skeleton(other.require_v_post_skeleton),
        require_skeleton_two_sided(other.require_skeleton_two_sided)
      {
        other.engine.split(engine);
      }

      void operator()(const tbb::blocked_range<std::size_t> &range)
      {
        // for all partitions in tbb-range
        for(std::size_t p = range.begin(); p != range.end(); ++p)
          (*this)(p);
      }

    protected:
      void operator()(std::size_t p)
      {
        // ICC has problems when the partition is a temporary
        auto partition = assembler.partitioning->partition(p);
        // Traverse partition
        for(auto it = partition.begin(); it != partition.end(); ++it)
        {
          const auto &e = *it;
          // Compute unique id
          const typename GV::IndexSet::IndexType ids = cell_mapper.map(e);

          ElementGeometry<Element> eg(e);

          if(engine.assembleCell(eg))
            continue;

          // Bind local test function space to element
          lfsv.bind( e );
          lfsv_cache.update();

          // Notify assembler engine about bind
          engine.onBindLFSV(eg,lfsv_cache);

          // Volume integration
          engine.assembleVVolume(eg,lfsv_cache);

          // Bind local trial function space to element
          lfsu.bind( e );
          lfsu_cache.update();

          // Notify assembler engine about bind
          engine.onBindLFSUV(eg,lfsu_cache,lfsv_cache);

          // Load coefficients of local functions
          engine.loadCoefficientsLFSUInside(lfsu_cache);

          // Volume integration
          engine.assembleUVVolume(eg,lfsu_cache,lfsv_cache);

          // Skip if no intersection iterator is needed
          if (require_uv_skeleton || require_v_skeleton ||
              require_uv_boundary || require_v_boundary ||
              require_uv_processor || require_v_processor)
          {
            // Traverse intersections
            unsigned int intersection_index = 0;

            const auto &irange = intersections(assembler.gfsu.gridView(), e);
            for(auto iit = irange.begin(); iit != irange.end(); ++iit)
            {
              const auto &is = *iit;
              IntersectionGeometry<Intersection> ig(is,intersection_index);

              switch (IntersectionType::get(is))
              {
              case IntersectionType::skeleton:
                // the specific ordering of the if-statements in the old code
                // caused periodic boundary intersection to be handled the
                // same as skeleton intersections
              case IntersectionType::periodic:
                if (require_uv_skeleton || require_v_skeleton)
                {
                  // compute unique id for neighbor
                  const typename GV::IndexSet::IndexType idn =
                    cell_mapper.map(*is.outside());

                  // Visit face if id is bigger
                  bool visit_face = ids > idn || require_skeleton_two_sided;

                  // unique vist of intersection
                  if (visit_face)
                  {
                    // Bind local test space to neighbor element
                    lfsvn.bind(*is.outside());
                    lfsvn_cache.update();

                    // Notify assembler engine about binds
                    engine.onBindLFSVOutside(ig,lfsv_cache,lfsvn_cache);

                    // Skeleton integration
                    engine.assembleVSkeleton(ig,lfsv_cache,lfsvn_cache);

                    if(require_uv_skeleton)
                    {
                      // Bind local trial space to neighbor element
                      lfsun.bind(*is.outside());
                      lfsun_cache.update();

                      // Notify assembler engine about binds
                      engine.onBindLFSUVOutside(ig, lfsu_cache, lfsv_cache,
                                                lfsun_cache, lfsvn_cache);

                      // Load coefficients of local functions
                      engine.loadCoefficientsLFSUOutside(lfsun_cache);

                      // Skeleton integration
                      engine.assembleUVSkeleton(ig, lfsu_cache, lfsv_cache,
                                                lfsun_cache, lfsvn_cache);

                      // Notify assembler engine about unbinds
                      engine.onUnbindLFSUVOutside(ig, lfsu_cache, lfsv_cache,
                                                  lfsun_cache, lfsvn_cache);
                    }

                    // Notify assembler engine about unbinds
                    engine.onUnbindLFSVOutside(ig,lfsv_cache,lfsvn_cache);
                  }
                }
                break;

              case IntersectionType::boundary:
                if(require_uv_boundary || require_v_boundary )
                {
                  // Boundary integration
                  engine.assembleVBoundary(ig,lfsv_cache);

                  if(require_uv_boundary)
                  {
                    // Boundary integration
                    engine.assembleUVBoundary(ig,lfsu_cache,lfsv_cache);
                  }
                }
                break;

              case IntersectionType::processor:
                if(require_uv_processor || require_v_processor )
                {
                  // Processor integration
                  engine.assembleVProcessor(ig,lfsv_cache);

                  if(require_uv_processor)
                  {
                    // Processor integration
                    engine.assembleUVProcessor(ig,lfsu_cache,lfsv_cache);
                  }
                }
                break;
              } // switch

              ++intersection_index;
            } // intersections
          } // do skeleton

          if(require_uv_post_skeleton || require_v_post_skeleton)
          {
            // Volume integration
            engine.assembleVVolumePostSkeleton(eg,lfsv_cache);

            if(require_uv_post_skeleton)
            {
              // Volume integration
              engine.assembleUVVolumePostSkeleton(eg, lfsu_cache, lfsv_cache);
            }
          }

          // Notify assembler engine about unbinds
          engine.onUnbindLFSUV(eg,lfsu_cache,lfsv_cache);

          // Notify assembler engine about unbinds
          engine.onUnbindLFSV(eg,lfsv_cache);

        } // partition
      }

    public:
      //! join another body
      /**
       * This calls join() first on the engines and then on the local
       * assemblers.
       */
      void join(AssembleBody &other)
      {
        other.engine.join(engine);
        TBBAssemblerSplit<LocalAssembler>::join(*localAssembler,
                                                *other.localAssembler);
      }
    };


    //! Assembler for TBB with coloring
    /**
     * \tparam Coloring            Coloring information for the partitioning.
     * \tparam Partitioning        Partitioning to use for parallel traversal.
     * \tparam GFSU                GridFunctionSpace for ansatz functions
     * \tparam GFSV                GridFunctionSpace for test functions
     * \tparam CU                  Constraints of ansatz functions
     * \tparam CV                  Constraints of test functions
     * \tparam nonoverlapping_mode Indicates whether assembling is done for
     *                             overlap cells
     */
    template<typename Coloring, typename Partitioning, typename GFSU,
             typename GFSV, typename CU, typename CV,
             bool nonoverlapping_mode=false>
    class ColoredTBBAssembler :
      public TBBAssembler<Partitioning, GFSU, GFSV, CU, CV,
                          nonoverlapping_mode>
    {
      typedef TBBAssembler<Partitioning, GFSU, GFSV, CU, CV,
                           nonoverlapping_mode> Base;
      template<class EngineFactory, class LocalAssembler>
      class AssembleBody;

    public:
      ColoredTBBAssembler(const GFSU& gfsu_, const GFSV& gfsv_, const CU& cu_,
                          const CV& cv_) :
        Base(gfsu_, gfsv_, cu_, cv_)
      { }

      ColoredTBBAssembler(const GFSU& gfsu_, const GFSV& gfsv_) :
        Base(gfsu_, gfsv_)
      { }

      void setColoring(const shared_ptr<const Coloring> &coloring)
      {
        coloring_ = coloring;
      }

      template<class EngineFactory, class LocalAssembler>
      void assemble(const EngineFactory &engineFactory,
                    LocalAssembler &la) const
      {
        typedef typename Base::Element Element;

        auto &assembler_engine = engineFactory(la);

        // Notify assembler engine about oncoming assembly
        assembler_engine.preAssembly();

        // Map each cell to unique id
        ElementMapper<typename Base::GV> cell_mapper(this->gfsu.gridView());

        // enumerate representatives
        std::vector<std::vector<std::size_t> >
          representatives(coloring_->colors());
        {
          std::vector<std::size_t> colorsize(coloring_->colors(), 0);
          for(typename Partitioning::Size p = 0;
              p < this->partitioning->partitions();
              ++p)
          {
            ++colorsize[coloring_->color(p)];
          }
          for(std::size_t c = 0; c < colorsize.size(); ++c)
            representatives[c].reserve(colorsize[c]);
        }
        for(typename Partitioning::Size p = 0;
            p < this->partitioning->partitions();
            ++p)
        {
          representatives[coloring_->color(p)].push_back(p);
        }

        bool parallel = TBBAssemblerSplit<LocalAssembler>::value
          && assembler_engine.threadSafe();
        if(this->verbosity > 0 && this->gfsu.gridView().comm().rank() == 0)
        {
          if(parallel)
            std::cout << "ColoredTBBAssembler: parallel iteration"
                      << std::endl;
          else
            std::cout << "ColoredTBBAssembler: sequential iteration"
                      << std::endl;
        }

        // loop over colors
        for(const auto &cr : representatives)
        {
          AssembleBody<EngineFactory, LocalAssembler>
            body(*this, engineFactory, la, assembler_engine, cell_mapper, cr);
          tbb::blocked_range<std::size_t> range(0, cr.size());

          if(parallel)
            tbb::parallel_reduce(range, body);
          else
            body(range);
        }

        // Notify assembler engine that assembly is finished
        assembler_engine.postAssembly(this->gfsu,this->gfsv);
      }

    private:
      shared_ptr<const Coloring> coloring_;
    };

    template<typename Coloring, typename Partitioning, typename GFSU,
             typename GFSV, typename CU, typename CV,
             bool nonoverlapping_mode>
    template<class EngineFactory, class LocalAssembler>
    class ColoredTBBAssembler<Coloring, Partitioning, GFSU, GFSV, CU, CV,
                              nonoverlapping_mode>::AssembleBody :
      public ColoredTBBAssembler::Base::
    template AssembleBody<EngineFactory, LocalAssembler>
    {
      typedef typename ColoredTBBAssembler::Base::
        template AssembleBody<EngineFactory, LocalAssembler> Base;

      typedef typename Base::LocalAssemblerEngine LocalAssemblerEngine;

    public:
      typedef typename ColoredTBBAssembler::GV GV;
      typedef typename Base::Assembler Assembler;
      AssembleBody(const Assembler &assembler_,
                   const EngineFactory &engineFactory_,
                   LocalAssembler &localAssembler_,
                   LocalAssemblerEngine &engine_,
                   const ElementMapper<GV> &cell_mapper_,
                   const std::vector<std::size_t> &representatives) :
        Base(assembler_, engineFactory_, localAssembler_, engine_,
             cell_mapper_),
        representatives_(representatives)
      { }

      AssembleBody(AssembleBody &other, tbb::split split) :
        Base(other, split),
        representatives_(other.representatives_)
      { }

      void operator()(const tbb::blocked_range<std::size_t> &range)
      {
        // for all partitions in tbb-range
        for(std::size_t p = range.begin(); p != range.end(); ++p)
          Base::operator()(representatives_[p]);
      }

    private:
      const std::vector<std::size_t> &representatives_;
    };

  }
}
#endif //  DUNE_PDELAB_TBB_ASSEMBLER_HH
