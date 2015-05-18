// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_COMMON_BORDERINDEXIDCACHE_HH
#define DUNE_PDELAB_COMMON_BORDERINDEXIDCACHE_HH

#include <vector>
#include <utility>
#include <unordered_map>

#include <dune/common/typetraits.hh>
#include <dune/geometry/typeindex.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/capabilities.hh>

namespace Dune {
  namespace PDELab {


    //! \addtogroup Backend
    //! \ingroup PDELab
    //! \{


    template<typename GFS>
    struct BorderIndexIdCache
    {

      typedef GFS GridFunctionSpace;
      typedef typename GFS::Traits::GridView GridView;
      typedef typename GridView::Grid Grid;

      typedef std::size_t size_type;
      typedef typename GFS::Traits::GridView::IndexSet::IndexType index_type;
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
        , _grid_view(gfs.gridView())
      {
        update();
      }

      void update()
      {
        _border_entities.resize(GlobalGeometryTypeIndex::size(Grid::dimension));
        _index_to_id.resize(GlobalGeometryTypeIndex::size(Grid::dimension));

        const typename GridView::IndexSet& index_set = _grid_view.indexSet();

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
      GridView _grid_view;
      BorderEntitySet _border_entities;
      IndexToIdMap _index_to_id;
      IdToIndexMap _id_to_index;

      template<int codim>
      typename enable_if<
        (codim > 0) && Capabilities::hasEntity<Grid,codim>::v
        >::type
      create_for_codim()
      {
        const typename GridView::IndexSet& index_set = _grid_view.indexSet();
        const typename Grid::GlobalIdSet& id_set = _grid_view.grid().globalIdSet();

        if (_gfs.ordering().contains(codim))
          {
            typedef typename GridView::template Codim<codim>::template Partition<InteriorBorder_Partition>::Iterator EntityIterator;
            for (EntityIterator it = _grid_view.template begin<codim,InteriorBorder_Partition>(),
                   end_it = _grid_view.template end<codim,InteriorBorder_Partition>();
                 it != end_it;
                 ++it)
              {
                index_type index = index_set.index(*it);
                size_type gt_index = GlobalGeometryTypeIndex::index(it->type());

                bool border_entity = _border_entities[gt_index][index] = (it->partitionType() == BorderEntity);
                if (!border_entity)
                  continue;

                id_type id = id_set.id(*it);

                _index_to_id[gt_index][index] = id;
                _id_to_index[id] = EntityIndex(gt_index,index);
              }
          }
        create_for_codim<codim-1>();
      }

      template<int codim>
      typename enable_if<
        (codim > 0) && !Capabilities::hasEntity<Grid,codim>::v
        >::type
      create_for_codim()
      {
        if (_gfs.ordering().contains(codim))
          DUNE_THROW(Dune::Exception,"Required codim " << codim << " not supported by grid!");
        create_for_codim<codim-1>();
      }

      template<int codim>
      typename enable_if<
        (codim == 0)
        >::type
      create_for_codim()
      {}

    };

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_COMMON_BORDERINDEXIDCACHE_HH
