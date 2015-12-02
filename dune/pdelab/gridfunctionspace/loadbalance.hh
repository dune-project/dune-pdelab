// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_LOADBALANCE_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_LOADBALANCE_HH

#include<dune/pdelab/common/polymorphicbufferwrapper.hh>
#include <dune/pdelab/gridfunctionspace/entityindexcache.hh>

namespace Dune {
  namespace PDELab {

    /*! \brief Data Handle for dofs communication while load blancing.
     *
     * \tparam GFS Type of ansatz space we need to restore after load balancing.
     * \tparam MAP Type of map that stores ids and dofs.
     */
    template<typename GFS, typename MAP>
    class LoadBalanceDataHandle
      : public Dune::CommDataHandleIF< LoadBalanceDataHandle<GFS,MAP>, double>
    {
    public:

      using R = typename MAP::mapped_type::value_type;
      using LFS = Dune::PDELab::LocalFunctionSpace<GFS>;
      using LFSCache = Dune::PDELab::LFSIndexCache<LFS>;
      using EntitySet = typename GFS::Traits::EntitySet;
      using IDSet = typename EntitySet::Traits::GridView::Grid::LocalIdSet;

      LoadBalanceDataHandle(const GFS& gfs, MAP& map)
        : gfs_(gfs)
        , lfs_(gfs)
        , lfsCache_(lfs_)
        , map_(map)
        , idSet_(gfs.entitySet().gridView().grid().localIdSet())
      {}

      //! Returns true if data for this codim should be communicated.
      bool contains (int dim, int codim) const
      {
        return true;
      }

      //! Returns true if size per entity of given dim and codim is a constant.
      bool fixedsize (int dim, int codim) const
      {
        return gfs_.finiteElementMap().fixedSize();
      }

      //! Size of data we communicate for a given entity.
      template<class EntityType>
      size_t size (EntityType& e) const
      {
        size_t commSize(0.0);

        // Only communicate if there are dofs for this codimension.
        if (gfs_.finiteElementMap().hasDOFs(e.codimension)){
          // Only communicate if e is in entitySet.
          if (gfs_.entitySet().contains(e)){
            // If size is for geometry type is fixed the fem knows the number
            // of dofs.  Otherwise we need to compute it.
            if (gfs_.finiteElementMap().fixedSize()){
              commSize = gfs_.finiteElementMap().size(e.geometry().type())*sizeof(R);
            }
            else{
              // Get size from map.
              auto find = map_.find(idSet_.id(e));
              commSize = find->second.size()*sizeof(double);
            }
          }
        }
        return commSize;
      }

      //! Send dofs if entitiy is in entitySet.
      template<class MessageBuffer, class EntityType>
      void gather (MessageBuffer& buff, const EntityType& e) const
      {
        // Communicate different things than char.
        using Buf = Dune::PDELab::PolymorphicBufferWrapper<MessageBuffer>;
        Buf bufWrapper(buff);

        // Communication is called for every level and every entity of
        // every codim of all entities where load balance will result in
        // changes.   We only send dofs for current gridView/entitySet.
        if (gfs_.finiteElementMap().hasDOFs(e.codimension)){
          if (gfs_.entitySet().contains(e)){
            // Find entity id in map.
            auto find = map_.find(idSet_.id(e));
            assert(find!=map_.end());

            // Send dofs.
            for (size_t i=0; i<find->second.size(); ++i){
              bufWrapper.write(find->second[i]);
            }
          }
        }
      }

      /*! \brief Receive dofs and store them in map.
       *
       * \param buff The message buffer.
       * \param e An entity.
       * \param n Size of data we receive.
       */
      template<class MessageBuffer, class EntityType>
      void scatter (MessageBuffer& buff, const EntityType& e, size_t n)
      {
        // Communicate different things than char.
        using Buf = Dune::PDELab::PolymorphicBufferWrapper<MessageBuffer>;
        Buf bufWrapper(buff);

        // Number of R's we receive.
        size_t numberOfEntries = n/sizeof(R);
        assert(sizeof(R)*numberOfEntries==n);

        // Create vector of dofs.
        std::vector<R> dofs(numberOfEntries);
        for (size_t i=0; i<numberOfEntries; ++i){
          bufWrapper.read(dofs[i]);
        }

        // Store id and dofs in map.
        const auto& id_set = gfs_.entitySet().grid().localIdSet();
        map_.insert({{id_set.id(e),dofs}});
      }

    private:
      const GFS& gfs_;
      LFS lfs_;
      LFSCache lfsCache_;
      MAP& map_;
      const IDSet& idSet_;
    };


    //! Fill map with id and correspding dofs for every entity of codim.
    template <typename GFS, typename V, typename MAP, int codim>
    void loadBalanceMapFiller (const GFS& gfs, V& v, MAP& map)
    {
      using IndexCache = Dune::PDELab::EntityIndexCache<GFS>;
      using LocalView =  typename V::template LocalView<IndexCache>;
      IndexCache indexCache(gfs);
      LocalView localView(v);
      const auto& id_set = gfs.entitySet().grid().localIdSet();

      // Iterate over all interiorBorder entities of codim.
      for (const auto& e : entities(gfs.entitySet(),Dune::Codim<codim>(),Dune::Partitions::interiorBorder)){
        // Bind cache to entity.
        indexCache.update(e);
        localView.bind(indexCache);

        // Store dofs of entity in std::vector.
        std::vector<typename LocalView::ElementType> dofs;
        for (std::size_t i=0; i<localView.size(); ++i){
          dofs.push_back(localView[i]);
        }

        // Insert id and vector in map.
        map.insert ( {{id_set.id(e),dofs}});

        // Unbind cache.
        localView.unbind();
      }
    }

    //! For every codim: Fill map with id of every entity and corresponding dofs.
    template <int codim>
    struct FillLoadBalanceDOFMap
    {
      template <typename GFS, typename V, typename MAP>
      static void fillMap (const GFS& gfs, V& v, MAP& map)
      {
        if (gfs.finiteElementMap().hasDOFs(codim)){
          loadBalanceMapFiller<GFS,V,MAP,codim>(gfs,v,map);
        }
        FillLoadBalanceDOFMap<codim-1>::fillMap(gfs,v,map);
      }
    };
    template <>
    struct FillLoadBalanceDOFMap<0>
    {
      template <typename GFS, typename V, typename MAP>
      static void fillMap (const GFS& gfs, V& v, MAP& map)
      {
        if (gfs.finiteElementMap().hasDOFs(0)){
          loadBalanceMapFiller<GFS,V,MAP,0>(gfs,v,map);
        }
      }
    };

    //! Store dofs for every entity of codim in map.
    template <typename GFS, typename V, typename MAP, int codim>
    void loadBalanceMapReader (const GFS& gfs, V& v, MAP& map)
    {
      using IndexCache = Dune::PDELab::EntityIndexCache<GFS>;
      using LocalView =  typename V::template LocalView<IndexCache>;
      IndexCache indexCache(gfs);
      LocalView localView(v);
      const auto& id_set = gfs.entitySet().grid().localIdSet();

      // Iterate over all interiorBorder entities of codim.
      for (const auto& e : entities(gfs.entitySet(),Dune::Codim<codim>(),Dune::Partitions::interiorBorder)){
        // Bind cache to entity.
        indexCache.update(e);
        localView.bind(indexCache);

        // Find key in map and get vector of dofs.
        auto find = map.find(id_set.id(e));
        auto& dofs(find->second);

        // Assert that we found element and that sizes of dof vector and local view match.
        assert(find!=map.end());
        assert(dofs.size()==localView.size());

        // Store dofs in local view.
        for (std::size_t i=0; i<dofs.size(); ++i){
          localView[i]=dofs[i];
        }

        // Write changes to underlying container.
        localView.commit();

        // Unbind cache.
        localView.unbind();
      }
    }

    //! For every codim: Update dofs for every entity from map.
    template <int codim>
    struct ReadLoadBalanceDOFMap
    {
      template <typename GFS, typename V, typename MAP>
      static void readMap (const GFS& gfs, V& v, MAP& map)
      {
        if (gfs.finiteElementMap().hasDOFs(codim)){
          loadBalanceMapReader<GFS,V,MAP,codim>(gfs,v,map);
        }
        ReadLoadBalanceDOFMap<codim-1>::readMap(gfs,v,map);
      }
    };
    template <>
    struct ReadLoadBalanceDOFMap<0>
    {
      template <typename GFS, typename V, typename MAP>
      static void readMap (const GFS& gfs, V& v, MAP& map)
      {
        if (gfs.finiteElementMap().hasDOFs(0)){
          loadBalanceMapReader<GFS,V,MAP,0>(gfs,v,map);
        }
      }
    };

    /*! \brief Load balance grid and restore gridfunctionspace and gridfunction.
     *
     * \note The grid must allow load balancing after the initial balance.
     *
     * \tparam Grid Type of the grid.
     * \tparam GFS Type of ansatz space we need to restore after load balancing.
     * \tparam X Container class for DOF vector.
     */
    template <typename Grid, typename GFS, typename X>
    void loadBalanceGrid(Grid& grid, GFS& gfs, X& x)
    {
      using IDSet = typename Grid::LocalIdSet;
      using ID = typename IDSet::IdType;
      using R = typename X::field_type;

      // Store vector of dofs for every entitiy of every codim.
      using MAP = std::unordered_map<ID,std::vector<R>>;
      MAP map;
      FillLoadBalanceDOFMap<GFS::Traits::GridView::dimension>::fillMap(gfs,x,map);

      // Create data handle and use load balancing with communication.
      LoadBalanceDataHandle<GFS,MAP> dh(gfs,map);
      grid.loadBalance(dh);

      // After load balancing the gfs has to be updated. We need to force
      // an update of the entity set, therefore: true.
      gfs.update(true);

      // Restore solution from dof map.
      x=X(gfs,0.0);
      ReadLoadBalanceDOFMap<GFS::Traits::GridView::dimension>::readMap(gfs,x,map);
    }

  } // namespace PDELab
} // namespace Dune

#endif
