// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_COMMON_BORDERINDEXIDCACHE_HH
#define DUNE_PDELAB_COMMON_BORDERINDEXIDCACHE_HH

#include <vector>
#include <utility>
#include <unordered_map>

#include <dune/common/typetraits.hh>
#include <dune/geometry/dimension.hh>
#include <dune/geometry/typeindex.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/partitionset.hh>

namespace Dune {
  namespace PDELab {


    //! \addtogroup Backend
    //! \ingroup PDELab
    //! \{


    template<typename GFS>
    struct BorderIndexIdCache
    {

      typedef GFS GridFunctionSpace;
      using EntitySet = typename GridFunctionSpace::Traits::EntitySet;
      typedef typename GFS::Traits::GridView GridView;
      typedef typename GridView::Grid Grid;

      typedef std::size_t size_type;
      using index_type = typename EntitySet::Traits::Index;
      typedef typename GFS::Traits::GridView::Grid::GlobalIdSet::IdType id_type;


      struct EntityIndex
        : public std::pair<std::size_t,std::size_t>
      {

        typedef std::size_t size_type;

        EntityIndex()
        {}

        EntityIndex(size_type gt_index, size_type entity_index)
          : std::pair<size_type,size_type>(gt_index,entity_index)
        {}

        size_type geometryTypeIndex() const
        {
          return this->first;
        }

        size_type entityIndex() const
        {
          return this->second;
        }

      };


      typedef std::vector<
        std::vector<
          bool
          >
        > BorderEntitySet;

      typedef std::vector<
        std::unordered_map<
          index_type,
          id_type
          >
        > IndexToIdMap;

      typedef std::unordered_map<
        id_type,
        EntityIndex
        > IdToIndexMap;

      BorderIndexIdCache(const GFS& gfs)
        : _gfs(gfs)
        , _entity_set(gfs.entitySet())
      {
        update();
      }

      void update()
      {
        _border_entities.resize(GlobalGeometryTypeIndex::size(Grid::dimension));
        _index_to_id.resize(GlobalGeometryTypeIndex::size(Grid::dimension));

        auto& index_set = _entity_set.indexSet();

        // Skip codim 0 - cells can't ever be border entities
        for (int codim = 1; codim <= Grid::dimension; ++codim)
          {
            if (!_gfs.ordering().contains(codim))
              continue;

            for (auto gt : index_set.types(codim))
              {
                _border_entities[GlobalGeometryTypeIndex::index(gt)].resize(index_set.size(gt));
                _index_to_id[GlobalGeometryTypeIndex::index(gt)];
              }
          }
        create_for_codim<Grid::dimension>();
      }

      bool isBorderEntity(std::size_t gt_index, std::size_t entity_index) const
      {
        return _border_entities[gt_index][entity_index];
      }

      id_type id(std::size_t gt_index,index_type entity_index) const
      {
        typename IndexToIdMap::value_type::const_iterator it = _index_to_id[gt_index].find(entity_index);
        if (it == _index_to_id[gt_index].end())
          {
            DUNE_THROW(Dune::Exception,"invalid argument (entity not in map)");
          }
        return it->second;
      }

      EntityIndex index(id_type entity_id) const
      {
        typename IdToIndexMap::const_iterator it = _id_to_index.find(entity_id);
        if (it == _id_to_index.end())
          {
            DUNE_THROW(Dune::Exception,"invalid argument (entity not in map)");
          }
        return it->second;
      }

      std::pair<bool,EntityIndex> findIndex(id_type entity_id) const
      {
        typename IdToIndexMap::const_iterator it = _id_to_index.find(entity_id);
        if (it == _id_to_index.end())
          return std::make_pair(false,EntityIndex());
        else
          return std::make_pair(true,it->second);
      }

    private:

      const GFS& _gfs;
      EntitySet _entity_set;
      BorderEntitySet _border_entities;
      IndexToIdMap _index_to_id;
      IdToIndexMap _id_to_index;

      template<int codim>
      typename std::enable_if<
        (codim > 0) && Capabilities::hasEntity<Grid,codim>::v
        >::type
      create_for_codim()
      {
        auto& index_set = _entity_set.indexSet();
        auto& id_set = _entity_set.gridView().grid().globalIdSet();

        if (_gfs.ordering().contains(codim))
          {
            for (const auto& e : entities(_entity_set,Codim<codim>{},Partitions::interiorBorder))
              {
                index_type index = index_set.index(e);
                size_type gt_index = GlobalGeometryTypeIndex::index(e.type());

                bool border_entity = _border_entities[gt_index][index] = (e.partitionType() == BorderEntity);
                if (!border_entity)
                  continue;

                id_type id = id_set.id(e);

                _index_to_id[gt_index][index] = id;
                _id_to_index[id] = EntityIndex(gt_index,index);
              }
          }
        create_for_codim<codim-1>();
      }

      template<int codim>
      typename std::enable_if<
        (codim > 0) && !Capabilities::hasEntity<Grid,codim>::v
        >::type
      create_for_codim()
      {
        if (_gfs.ordering().contains(codim))
          DUNE_THROW(Dune::Exception,"Required codim " << codim << " not supported by grid!");
        create_for_codim<codim-1>();
      }

      template<int codim>
      typename std::enable_if<
        (codim == 0)
        >::type
      create_for_codim()
      {}

    };

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_COMMON_BORDERINDEXIDCACHE_HH
