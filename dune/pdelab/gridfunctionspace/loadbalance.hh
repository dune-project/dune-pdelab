// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_LOADBALANCE_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_LOADBALANCE_HH

#include <dune/geometry/dimension.hh>
#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/common/partitionset.hh>

#include<dune/pdelab/common/polymorphicbufferwrapper.hh>
#include <dune/pdelab/gridfunctionspace/entityindexcache.hh>

namespace Dune {
  namespace PDELab {

#ifndef DOXYGEN
    namespace impl{

      /*! \brief Data Handle for dofs communication while load blancing.
       *
       * \tparam T... GFSandMap structs passed as tuple elements
       */
      template<typename... T>
      class LoadBalanceDataHandle
        : public Dune::CommDataHandleIF<LoadBalanceDataHandle<T...>,
                                        typename std::decay<typename std::tuple_element<0,std::tuple<T...>>::type>::type
                                        ::Map::mapped_type::value_type>
      {
      public:

        using R0 = typename std::decay<typename std::tuple_element<0,std::tuple<T...>>::type>::type
          ::Map::mapped_type::value_type;

        LoadBalanceDataHandle(std::tuple<T&...>&& mapTuple)
          : _mapTuple(mapTuple)
        {
        }

        //! Returns true if data for this codim should be communicated.
        bool contains (int dim, int codim) const
        {
          return true;
        }

        //! Returns true if size per entity of given dim and codim is a constant.
        bool fixedsize (int dim, int codim) const
        {
          // We return false here, since gather and scatter is called
          // for all entities of all levels of the grid but we only
          // communicate data for a given gridView. This means there
          // will be gather and scatter calls for cells that don't exist
          // in our current grid view and we won't communicate anything
          // for these cells.
          return false;
        }

        // End of tmp recursion
        template <std::size_t I, typename EntityType>
        inline typename std::enable_if<I==sizeof...(T),void>::type
        sizeTMP(std::size_t& commSize, EntityType& e) const
        {
        }

        // Go through tuple and accumulate communication size
        template <std::size_t I=0, typename EntityType>
        inline typename std::enable_if<I<sizeof...(T),void>::type
        sizeTMP(std::size_t& commSize, EntityType& e) const
        {
          // Compare if data type for this entry equals data typo of first entry
          using R = typename std::decay<typename std::tuple_element<I,std::tuple<T...>>::type>::type
            ::Map::mapped_type::value_type;
          static_assert(std::is_same<R,R0>::value,"Different field type of vectors not supported by DH");

          // Current function space
          auto& gfs = std::get<I>(_mapTuple)._gfs;

          // Only communicate if there are dofs for this codimension.
          if (gfs.finiteElementMap().hasDOFs(e.codimension)){

            // We first have to communicate the size of data we send for this particular vector
            commSize += sizeof(R);

            // Only communicate data if entity lies in our entitySet
            if (gfs.entitySet().contains(e)){
              // Get some types
              using GFS = typename std::decay<decltype(gfs)>::type;
              using EntitySet = typename GFS::Traits::EntitySet;
              using IDSet = typename EntitySet::Traits::GridView::Grid::LocalIdSet;

              // Find element id in map
              const IDSet& idSet(gfs.entitySet().gridView().grid().localIdSet());
              auto& map = std::get<I>(_mapTuple)._map;
              auto find = map.find(idSet.id(e));
              assert (find!=map.end());

              // Add size of degrees of freedom vetor
              commSize += find->second.size()*sizeof(R);
            }
          }

          // Get size for next vector
          sizeTMP<I+1>(commSize,e);
        }

        //! Size of data we communicate for a given entity.
        template<class EntityType>
        std::size_t size (EntityType& e) const
        {
          // Use tmp to sum up sizes for all vectors
          std::size_t commSize(0.0);
          sizeTMP<0>(commSize,e);
          return commSize;
        }


        // End of tmp recursion
        template <std::size_t I, typename Buf, typename Entity>
        inline typename std::enable_if<I==sizeof...(T),void>::type
        gatherTMP(Buf& buf, const Entity& e) const
        {
        }

        // Go through tuple and call gather for all vectors
        template <std::size_t I=0, typename Buf, typename Entity>
        inline typename std::enable_if<I<sizeof...(T),void>::type
        gatherTMP(Buf& buf, const Entity& e) const
        {
          // Current gridFunctionSpace and dof map for the vector
          auto& gfs = std::get<I>(_mapTuple)._gfs;
          auto& map = std::get<I>(_mapTuple)._map;

          // Communication is called for every level and every entity of
          // every codim of all entities where load balance will result in
          // changes.   We only send dofs for current gridView/entitySet.
          if (gfs.finiteElementMap().hasDOFs(e.codimension)){
            if (gfs.entitySet().contains(e)){
              // Get important types
              using GFS = typename std::decay<decltype(gfs)>::type;
              using EntitySet = typename GFS::Traits::EntitySet;
              using IDSet = typename EntitySet::Traits::GridView::Grid::LocalIdSet;

              // Find element id in map
              const IDSet& idSet(gfs.entitySet().gridView().grid().localIdSet());

              // Find entity id in map.
              auto find = map.find(idSet.id(e));
              assert (find!=map.end());

              // Send size we need to communicate for this vector
              buf.write (static_cast<R0>(find->second.size()));

              // Send dofs
              for (size_t i=0; i<find->second.size(); ++i){
                buf.write(find->second[i]);
              }
            }
            else {
              // Only communicate that we don't communicate any DOFs
              R0 tmp(0);
              buf.write (tmp);
            }
          } // hasDOFs

          // Call gather for next vector
          gatherTMP<I+1> (buf,e);
        }

        //! Send dofs for all vectors from the tuple using tmp
        template<class MessageBuffer, class EntityType>
        void gather (MessageBuffer& buff, const EntityType& e) const
        {
          // Communicate different things than char.
          using Buf = Dune::PDELab::PolymorphicBufferWrapper<MessageBuffer>;
          Buf bufWrapper(buff);

          // Call gather for all vectors using tmp
          gatherTMP<0> (bufWrapper,e);
        }

        // End of tmp reucrsion
        template <std::size_t I, typename Buf, typename Entity>
        inline typename std::enable_if<I==sizeof...(T),void>::type
        scatterTMP(Buf& buf, const Entity& e) const
        {
        }

        // Go through tuple receive DOFs
        template <std::size_t I=0, typename Buf, typename Entity>
        inline typename std::enable_if<I<sizeof...(T),void>::type
        scatterTMP(Buf& buf, const Entity& e) const
        {
          auto& gfs = std::get<I>(_mapTuple)._gfs;
          auto& map = std::get<I>(_mapTuple)._map;
          if (gfs.finiteElementMap().hasDOFs(e.codimension)){

            // Receive number of DOFs for this vector
            R0 tmp;
            buf.read(tmp);
            std::size_t numberOfEntries(0);
            numberOfEntries = (size_t) tmp;

            // Create vector of DOFs and receive DOFs
            std::vector<R0> dofs(numberOfEntries);
            for (size_t i=0; i<numberOfEntries; ++i){
              buf.read(dofs[i]);
            }

            // Store id and dofs in map
            const auto& id_set = gfs.entitySet().grid().localIdSet();
            map.insert({{id_set.id(e),dofs}});
          }

          // Call scatter for next vetcor using tmp
          scatterTMP<I+1> (buf,e);
        }


        /*! \brief Receive dofs and store them in map.
         *
         * \param buff The message buffer.
         * \param e An entity.
         * \param n Size of data we receive for all vectors together.
         */
        template<class MessageBuffer, class EntityType>
        void scatter (MessageBuffer& buff, const EntityType& e, size_t n)
        {
          // Communicate different things than char.
          using Buf = Dune::PDELab::PolymorphicBufferWrapper<MessageBuffer>;
          Buf bufWrapper(buff);

          // Call scatter for all vectors using tmp
          scatterTMP<0> (bufWrapper, e);
        }

      private:
        // Tuple storing GFSandMap structs for every vector that should get adapted
        std::tuple<T&...>& _mapTuple;
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


      // Store reference to function space and map
      template <typename G, typename M>
      struct GFSAndMap
      {
        // Export types
        using GFS = G;
        using Map = M;

        GFSAndMap (GFS& gfs, Map& m) : _gfs(gfs), _map(m)
        {
        }

        GFS& _gfs;
        Map& _map;
      };

      // Create a GFSAndMap struct
      template <typename GFS, typename M>
      GFSAndMap<GFS,M> packGFSAndMap(GFS& gfs, M& m)
      {
        GFSAndMap<GFS,M> pack(gfs,m);
        return pack;
      }


      // Forward declarations needed for the tmp recursion
      template <typename Grid, typename... T>
      void iteratePacks(Grid& grid, std::tuple<T&...>&& mapTuple);
      template <typename Grid, typename... T, typename X, typename... XS>
      void iteratePacks(Grid& grid, std::tuple<T&...>&& mapTuple, X& x, XS&... xs);

      // This function is called after the last vector of the tuple.  Here
      // the next pack is called.  On the way back we update the current
      // function space.
      template<std::size_t I = 0, typename Grid, typename... T, typename X, typename... XS>
      inline typename std::enable_if<I == std::tuple_size<typename X::Tuple>::value, void>::type
      iterateTuple(Grid& grid, std::tuple<T&...>&& mapTuple, X& x, XS&... xs)
      {
        // Iterate next pack
        iteratePacks(grid,std::move(mapTuple),xs...);

        // On our way back we need to update the current function space
        x._gfs.update(true);
      }

      /* In this function we store the data of the current vector (indicated
       * by template parameter I) of the current pack. After recursively
       * iterating through the other packs and vectors we replay the data.
       *
       * @tparam I      std:size_t used for tmp
       * @tparam Grid   Grid type
       * @tparam T...   Types of tuple elements (for storing data transfer maps)
       * @tparam X      Current  pack
       * @tparam ...XS  Remaining packs
       */
      template<std::size_t I = 0, typename Grid, typename... T, typename X, typename... XS>
      inline typename std::enable_if<I < std::tuple_size<typename X::Tuple>::value, void>::type
      iterateTuple(Grid& grid, std::tuple<T&...>&& mapTuple, X& x, XS&... xs)
      {
        // Get some basic types
        using GFS = typename X::GFS;
        using Tuple = typename X::Tuple;
        using V  = typename std::decay<typename std::tuple_element<I,Tuple>::type>::type;
        using IDSet = typename Grid::LocalIdSet;
        using ID = typename IDSet::IdType;
        using R = typename V::field_type;

        // Store vector of dofs for every entitiy of every codim.
        using MAP = std::unordered_map<ID,std::vector<R>>;
        MAP map;
        FillLoadBalanceDOFMap<GFS::Traits::GridView::dimension>::fillMap(x._gfs,std::get<I>(x._tuple),map);

        // Pack gfs and tuple
        auto mapPack = packGFSAndMap(x._gfs,map);
        auto newMapTuple = std::tuple_cat(mapTuple,std::tie(mapPack));

        // Recursively iterate through remaining vectors (and packs). Load
        // balancing will be done at the end of recursion.
        iterateTuple<I+1>(grid,std::move(newMapTuple),x,xs...);

        // Restore solution from dof map.
        std::get<I>(x._tuple) = V(x._gfs,0.0);
        ReadLoadBalanceDOFMap<GFS::Traits::GridView::dimension>::readMap(x._gfs,std::get<I>(x._tuple),map);
      }

      template <typename... T>
      LoadBalanceDataHandle<T...> createLoadBalanceDataHandle (std::tuple<T&...>&& mapTuple)
      {
        LoadBalanceDataHandle<T...> dh(std::move(mapTuple));
        return dh;
      }

      // This gets called after the last pack.  After this function call we
      // have visited every vector of every pack and we will go back through
      // the recursive function calls.
      template <typename Grid, typename... T>
      void iteratePacks(Grid& grid, std::tuple<T&...>&& mapTuple)
      {
        // Create data handle and use load balancing with communication.
        auto dh = createLoadBalanceDataHandle(std::move(mapTuple));

        std::cout << "Calling load balance with data communication" << std::endl;
        grid.loadBalance(dh);
      }

      /* Use template meta programming to iterate over packs at compile time
       *
       * In order to adapt our grid and all vectors of all packs we need to
       * do the following:
       * - Iterate over all vectors of all packs.
       * - Store the data from the vectors where things could change.
       * - Load Balance our grid.
       * - Update function spaces and restore data.
       *
       * The key point is that we need the object that stores the data to
       * replay it.  Because of that we can not just iterate over the packs
       * and within each pack iterate over the vectors but we have to make
       * one big recursion.  Therefore we iterate over the vectors of the
       * current pack.
       */
      template <typename Grid, typename... T, typename X, typename... XS>
      void iteratePacks(Grid& grid, std::tuple<T&...>&& mapTuple, X& x, XS&... xs)
      {
        iterateTuple(grid,std::move(mapTuple),x,xs...);
      }

    } // namespace impl
#endif // DOXYGEN


    /*! \brief Load balance grid and restore gridfunctionspace and gridfunctions given as packs.
     *
     * \note This only works for scalar gridfunctionspaces.
     * \note The grid must allow load balancing after the initial balance.
     * \note A pack can be created using the transferSolution function.
     *
     * @tparam Grid   Type of the Grid
     * @tparam X      Packed GFS with vectors that should be adapted
     */
    template <typename Grid, typename... X>
    void loadBalanceGrid(Grid& grid, X&... x)
    {
      // Create tuple where all data transfer maps get stored
      std::tuple<> mapTuple;

      // Iterate over packs
      impl::iteratePacks(grid,std::move(mapTuple),x...);
    }

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_LOADBALANCE_HH
