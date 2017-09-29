#ifndef DUNE_PDELAB_COMMON_PARTITIONVIEWENTITYSET_HH
#define DUNE_PDELAB_COMMON_PARTITIONVIEWENTITYSET_HH

#include <cassert>
#include <vector>
#include <bitset>
#include <memory>
#include <algorithm>
#include <numeric>
#include <type_traits>

#include <dune/common/version.hh>
#include <dune/common/iteratorrange.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/typeindex.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/partitionset.hh>
#include <dune/grid/common/rangegenerators.hh>

namespace Dune {

  namespace PDELab {

    // import range generators to make sure they work with PartitionViewEntitySet
    using Dune::entities;
    using Dune::elements;
    using Dune::facets;
    using Dune::edges;
    using Dune::vertices;
    using Dune::descendantElements;
    using Dune::intersections;

    template<typename GV, typename P>
    class PartitionViewEntitySet;

    template<typename GV, typename P>
    class PartitionViewEntitySetIndexSet;

    template<typename GV, typename P>
    struct PartitionViewEntitySetTraits
    {

      using Partitions = typename std::decay<P>::type;

      using Grid = typename GV::Traits::Grid;
      using GridView = GV;
      using EntitySet = Dune::PDELab::PartitionViewEntitySet<GV,P>;
      using IndexSet = PartitionViewEntitySetIndexSet<GV,Partitions>;
      using BaseIndexSet = typename GV::Traits::IndexSet;

      using Element = typename GV::template Codim<0>::Entity;

      using Intersection = typename GV::Traits::Intersection;

      using IntersectionIterator = typename GV::Traits::IntersectionIterator;

      using CollectiveCommunication = typename GV::Traits::CollectiveCommunication;

      using size_type = std::size_t;
      using dim_type = int;

      using Index = typename BaseIndexSet::IndexType;

      using Types = IteratorRange<std::vector<GeometryType>::const_iterator>;

      using CodimMask = std::bitset<GV::dimension + 1>;

      using CoordinateField = typename Grid::ctype;

      constexpr static Index invalidIndex()
      {
        return ~static_cast<Index>(0ull);
      }

      static const bool conforming = GV::Traits::conforming;

      static const dim_type dimension = GV::dimension;

      static const dim_type dimensionworld = GV::dimensionworld;

      template<dim_type codim>
      struct Codim
      {

        using Iterator = typename GV::template Codim<codim>::template Partition<Partitions::partitionIterator()>::Iterator;

        using Entity = typename GV::template Codim<codim>::Entity;

        using Geometry = typename GV::template Codim<codim>::Geometry;

        using LocalGeometry = typename GV::template Codim<codim>::LocalGeometry;

        template<PartitionIteratorType pitype>
        struct Partition
        {

          using Iterator = typename GV::template Codim<codim>::template Partition<pitype>::Iterator;

        };

      };

    };


    template<typename GV, typename P>
    class PartitionViewEntitySet
    {

    public:

      using Traits = PartitionViewEntitySetTraits<GV,P>;

      using Partitions = typename Traits::Partitions;
      using Grid      = typename Traits::Grid;
      using GridView = typename Traits::GridView;
      using IndexSet  = typename Traits::IndexSet;
      using BaseIndexSet = typename Traits::BaseIndexSet;
      using Element = typename Traits::Element;
      using Intersection = typename Traits::Intersection;
      using IntersectionIterator = typename Traits::IntersectionIterator;
      using CollectiveCommunication = typename Traits::CollectiveCommunication;
      using CodimMask = typename Traits::CodimMask;
      using CoordinateField = typename Traits::CoordinateField;
      using size_type = typename Traits::size_type;
      using dim_type = typename Traits::dim_type;

      using ctype = CoordinateField;

      static const bool conforming = Traits::conforming;
      static const dim_type dimension = Traits::dimension;
      static const dim_type dimensionworld = Traits::dimensionworld;

      template<dim_type codim>
      using Codim = typename Traits::template Codim<codim>;

      constexpr static Partitions partitions()
      {
        return {};
      }

      constexpr static CodimMask allCodims()
      {
        return {~0ull};
      }

      const Grid& grid() const
      {
        return gridView().grid();
      }

      //! Returns the IndexSet of this EntitySet.
      const IndexSet& indexSet() const
      {
        return *_index_set;
      }

      //! Returns the IndexSet of the underlying GridView.
      const BaseIndexSet& baseIndexSet() const
      {
        return indexSet().baseIndexSet();
      }

      template<dim_type codim>
      typename Codim<codim>::Iterator
      begin() const
      {
        return gridView().template begin<codim,Partitions::partitionIterator()>();
      }

      template<dim_type codim>
      typename Codim<codim>::Iterator
      end() const
      {
        return gridView().template end<codim,Partitions::partitionIterator()>();
      }

      template<dim_type codim, PartitionIteratorType pitype>
      typename GV::template Codim<codim>::template Partition<pitype>::Iterator
      begin() const
      {
        return gridView().template begin<codim,pitype>();
      }

      template<dim_type codim, PartitionIteratorType pitype>
      typename GV::template Codim<codim>::template Partition<pitype>::Iterator
      end() const
      {
        return gridView().template end<codim,pitype>();
      }

      size_type size(dim_type codim) const
      {
        return indexSet().size(codim);
      }

      size_type size(const GeometryType& gt) const
      {
        return indexSet().size(gt);
      }

      template<typename Entity>
      bool contains(const Entity& e) const
      {
        return indexSet().contains(e);
      }

      bool contains(dim_type codim) const
      {
        return indexSet().contains(codim);
      }

      bool contains(const GeometryType& gt) const
      {
        return indexSet().contains(gt);
      }

      IntersectionIterator ibegin(const typename Codim<0>::Entity& entity) const
      {
        return gridView().ibegin(entity);
      }

      IntersectionIterator iend(const typename Codim<0>::Entity& entity) const
      {
        return gridView().iend(entity);
      }

      const CollectiveCommunication& comm() const
      {
        return gridView().comm();
      }

      //! Returns the overlap size of this EntitySet, which depends on its PartitionSet.
      size_type overlapSize(dim_type codim) const
      {
        return Partitions::contains(Dune::Partitions::overlap) ? gridView().overlapSize(codim) : 0;
      }

      //! Returns the ghost size of this EntitySet, which depends on its PartitionSet.
      size_type ghostSize(dim_type codim) const
      {
        return Partitions::contains(Dune::Partitions::ghost) ? gridView().ghostSize(codim) : 0;
      }

      template<typename DataHandle>
      void communicate(DataHandle& data, InterfaceType iftype, CommunicationDirection dir) const
      {
        gridView().communicate(data,iftype,dir);
      }

      //! Returns the underlying GridView.
      const GridView& gridView() const
      {
        return indexSet().gridView();
      }

      PartitionViewEntitySet(const GridView& gv, CodimMask supported_codims)
        : _index_set(std::make_shared<IndexSet>(gv,supported_codims,true))
      {}

      explicit PartitionViewEntitySet(const GridView& gv, bool initialize = true)
        : _index_set(std::make_shared<IndexSet>(gv,CodimMask(initialize ? ~0ull : 0ull),initialize))
      {}

      //! Reset this EntitySet, which removes all entities from it.
      void reset()
      {
        _index_set->reset();
      }

      //! Add all entities of the given codim to this EntitySet.
      void addCodim(dim_type codim)
      {
        _index_set->addCodim(codim);
      }

      //! Remove all entities of the given codim from this EntitySet.
      void removeCodim(dim_type codim)
      {
        _index_set->removeCodim(codim);
      }

      //! Returns true if you need to call update on this EntitySet before using it.
      bool needsUpdate() const
      {
        return _index_set->needsUpdate();
      }

      //! Update the internal state of this EntitySet.
      /**
       *
       * \param force   If true, forces an update even if the EntitySet parameters have not
       *                changed. This is e.g. required if the underlying grid has changed due
       *                to adaptivity.
       *
       * \return  Returns true if the state of the EntitySet was changed by this method.
       */
      bool update(bool force = false)
      {
        return _index_set->update(force);
      }

    private:

      std::shared_ptr<IndexSet> _index_set;

    };

    template<typename GV, typename P>
    class PartitionViewEntitySetIndexSetBase
    {

      template<typename,typename>
      friend class PartitionViewEntitySet;

    public:

      using Traits = PartitionViewEntitySetTraits<GV,P>;

      using Partitions = typename Traits::Partitions;
      using Grid      = typename Traits::Grid;
      using GridView = typename Traits::GridView;
      using BaseIndexSet = typename Traits::BaseIndexSet;
      using size_type = typename Traits::size_type;
      using dim_type = typename Traits::dim_type;
      using Index = typename Traits::Index;
      using Types = typename Traits::Types;
      using CodimMask = typename Traits::CodimMask;

      using IndexType = Index;

      constexpr static Index invalidIndex()
      {
        return Traits::invalidIndex();
      }

      template<dim_type codim>
      using Codim = typename Traits::template Codim<codim>;

      PartitionViewEntitySetIndexSetBase(const PartitionViewEntitySetIndexSetBase&) = delete;
      PartitionViewEntitySetIndexSetBase& operator=(const PartitionViewEntitySetIndexSetBase&) = delete;


    protected:

      bool update(bool force)
      {
        if (!(_needs_update || force))
          return false;
        std::fill(_gt_offsets.begin(),_gt_offsets.end(),0);
        std::fill(_mapped_gt_offsets.begin(),_mapped_gt_offsets.end(),0);
        _active_geometry_types.reset();
        _geometry_types.resize(0);
        for (dim_type codim = 0; codim <= GV::dimension; ++codim)
          {
            if (!_wanted_codims.test(codim))
              continue;
            for (const auto& gt : baseIndexSet().types(codim))
              {
                auto gt_index = GlobalGeometryTypeIndex::index(gt);
                _gt_offsets[gt_index + 1] = baseIndexSet().size(gt);
                _geometry_types.push_back(gt);
                _active_geometry_types.set(gt_index);
              }
          }
        for (dim_type codim = 0; codim <= GV::dimension; ++codim)
          {
            auto range = std::equal_range(
              _geometry_types.begin(),
              _geometry_types.end(),
              GeometryType(GeometryType::none,GV::dimension - codim),
              [](const GeometryType& x, const GeometryType& y)
              {
                // reverse order because we store in ascending order with regard to the codim, not the dim
                return y.dim() < x.dim();
              });
            _per_codim_geometry_types[codim] = {range.first,range.second};
          }

        std::partial_sum(_gt_offsets.begin(),_gt_offsets.end(),_gt_offsets.begin());
        _active_codims = _wanted_codims;
        _needs_update = false;
        return true;
      }

    public:

      size_type size(GeometryType gt) const
      {
        assert(!needsUpdate());
        auto gt_index = GlobalGeometryTypeIndex::index(gt);
        return _mapped_gt_offsets[gt_index + 1] - _mapped_gt_offsets[gt_index];
      }

      size_type size(dim_type codim) const
      {
        assert(!needsUpdate());
        auto dim = GV::dimension;
#if DUNE_VERSION_NEWER_REV(DUNE_GEOMETRY,2,4,1)
        return _mapped_gt_offsets[GlobalGeometryTypeIndex::offset(dim-codim+1)] -
          _mapped_gt_offsets[GlobalGeometryTypeIndex::offset(dim-codim)];
#else
        return codim < dim ? _mapped_gt_offsets[GlobalGeometryTypeIndex::size(dim-codim)] -
          _mapped_gt_offsets[GlobalGeometryTypeIndex::size(dim-codim-1)] : 0;
#endif
      }

      template<typename Entity>
      bool contains(const Entity& e) const
      {
        return Partitions::contains(e.partitionType()) ? baseIndexSet().contains(e) : false;
      }

      bool contains(dim_type codim) const
      {
        return _active_codims.test(codim);
      }

      bool contains(const GeometryType& gt) const
      {
        return _active_geometry_types.test(GlobalGeometryTypeIndex::index(gt));
      }

      const BaseIndexSet& baseIndexSet() const
      {
        return _gv.indexSet();
      }

      Types types(dim_type codim) const
      {
        assert(!needsUpdate());
        return _per_codim_geometry_types[codim];
      }

      Types types() const
      {
        assert(!needsUpdate());
        return {_geometry_types.begin(),_geometry_types.end()};
      }

      PartitionViewEntitySetIndexSetBase(const GV& gv, CodimMask wanted_codims)
        : _gv(gv)
        , _needs_update(true)
        , _wanted_codims(wanted_codims)
      {}

      const GridView& gridView() const
      {
        return _gv;
      }

      bool needsUpdate() const
      {
        return _needs_update;
      }

    protected:

      void reset()
      {
        _needs_update = true;
        _wanted_codims.reset();
      }

      void addCodim(dim_type codim)
      {
        _wanted_codims.set(codim);
        _needs_update = _wanted_codims != _active_codims || _wanted_codims.none();
      }

      void removeCodim(dim_type codim)
      {
        _wanted_codims.reset(codim);
        _needs_update = _wanted_codims != _active_codims || _wanted_codims.none();
      }

      GV _gv;
      bool _needs_update;
      CodimMask _wanted_codims;
      std::bitset<GlobalGeometryTypeIndex::size(GV::dimension)> _active_geometry_types;
      CodimMask _active_codims;
      std::array<size_type,GlobalGeometryTypeIndex::size(GV::dimension) + 1> _gt_offsets;
      std::array<size_type,GlobalGeometryTypeIndex::size(GV::dimension) + 1> _mapped_gt_offsets;

    private:

      std::vector<GeometryType> _geometry_types;
      std::array<Types,GV::dimension + 1> _per_codim_geometry_types;

    };

    template<typename GV, typename P>
    class PartitionViewEntitySetIndexSet
      : public PartitionViewEntitySetIndexSetBase<GV,P>
    {

      using Base = PartitionViewEntitySetIndexSetBase<GV,P>;

      template<typename,typename>
      friend class PartitionViewEntitySet;

    public:

      using typename Base::Index;
      using typename Base::Partitions;
      using typename Base::size_type;
      using typename Base::dim_type;

      using typename Base::Grid;

      using Base::gridView;
      using Base::baseIndexSet;
      using Base::invalidIndex;
      using Base::contains;
      using typename Base::CodimMask;
      using Base::needsUpdate;

    private:

      static constexpr bool hasAllEntities(Dune::Dim<Grid::dimension + 1>)
      {
        return true;
      }

      template<dim_type dim = 0>
      static constexpr bool hasAllEntities(Dune::Dim<dim> = {})
      {
        return Capabilities::hasEntity<Grid,dim>::v && hasAllEntities(Dune::Dim<dim+1>{});
      }

      bool update(bool force)
      {
        if (!Base::update(force))
          return false;
        _indices.assign(_gt_offsets.back(),invalidIndex());
        _mapped_gt_offsets[0] = 0;
        update_codims(std::integral_constant<bool,hasAllEntities()>{});
        std::partial_sum(_mapped_gt_offsets.begin(),_mapped_gt_offsets.end(),_mapped_gt_offsets.begin());
        return true;
      }

      void update_codims(std::true_type)
      {
        update_codim(Dune::Codim<0>{});
      }

      void update_codim(Dune::Codim<GV::dimension+1>)
      {}

      template<dim_type cd>
      void update_codim(Dune::Codim<cd> codim)
      {
        if (_active_codims.test(codim))
          for (const auto& e : entities(gridView(),codim,Dune::Partitions::all))
            {
              auto gt = e.type();
              auto gt_index = GlobalGeometryTypeIndex::index(gt);
              if (Partitions::contains(e.partitionType()))
                _indices[_gt_offsets[gt_index] + baseIndexSet().index(e)] = _mapped_gt_offsets[gt_index + 1]++;
            }
        update_codim(Dune::Codim<cd+1>{});
      }


      void update_codims(std::false_type)
      {
        std::fill(_indices.begin(),_indices.end(),invalidIndex());

        auto& index_set = baseIndexSet();

        for (const auto& e : elements(gridView(),Dune::Partitions::all))
          {
            if (!Partitions::contains(e.partitionType()))
              continue;

            auto ref_el = ReferenceElements<typename Base::Traits::CoordinateField,GV::dimension>::general(e.type());
            for (dim_type codim = 0; codim <= Grid::dimension; ++codim)
              {
                if (!_active_codims.test(codim))
                  continue;

                size_type sub_entity_count = ref_el.size(codim);

                for(size_type i = 0; i < sub_entity_count; ++i)
                  {
                    auto gt = ref_el.type(i,codim);
                    auto gt_index = GlobalGeometryTypeIndex::index(gt);
                    auto index = index_set.subIndex(e,i,codim);
                    if (_indices[_gt_offsets[gt_index] + index] == invalidIndex())
                      _indices[_gt_offsets[gt_index] + index] = _mapped_gt_offsets[gt_index + 1]++;
                  }
              }
          }
      }


    public:

      template<typename E>
      Index index(const E& e) const
      {
        assert(!needsUpdate());
        assert(Partitions::contains(e.partitionType()));
        assert(contains(e.type()));
        auto gt_index = GlobalGeometryTypeIndex::index(e.type());
        return _indices[_gt_offsets[gt_index] + baseIndexSet().index(e)];
      }

      template<typename E>
      Index subIndex(const E& e, size_type i, dim_type codim) const
      {
        assert(!needsUpdate());
        assert(Partitions::contains(e.partitionType()));
        auto gt = ReferenceElements<typename Base::Traits::CoordinateField,GV::dimension>::general(e.type()).type(i,codim);
        assert(contains(gt));
        auto gt_index = GlobalGeometryTypeIndex::index(gt);
        return _indices[_gt_offsets[gt_index] + baseIndexSet().subIndex(e,i,codim)];
      }


      template<typename E>
      Index uniqueIndex(const E& e) const
      {
        assert(!needsUpdate());
        assert(Partitions::contains(e.partitionType()));
        assert(contains(e.type()));
        auto gt_index = GlobalGeometryTypeIndex::index(e.type());
        return _indices[_gt_offsets[gt_index] + baseIndexSet().index(e)] + _mapped_gt_offsets[gt_index];
      }

      template<typename E>
      Index uniqueSubIndex(const E& e, size_type i, dim_type codim) const
      {
        assert(!needsUpdate());
        assert(Partitions::contains(e.partitionType()));
        auto gt = ReferenceElements<typename Base::Traits::CoordinateField,GV::dimension>::general(e.type()).type(i,codim);
        assert(contains(gt));
        auto gt_index = GlobalGeometryTypeIndex::index(gt);
        return _indices[_gt_offsets[gt_index] + baseIndexSet().subIndex(e,i,codim)] + _mapped_gt_offsets[gt_index];
      }


      PartitionViewEntitySetIndexSet(const GV& gv, CodimMask wanted_codims, bool initialize)
        : Base(gv,wanted_codims)
      {
        if (initialize)
          update(true);
      }

    private:

      using Base::_active_codims;
      using Base::_gt_offsets;
      using Base::_mapped_gt_offsets;

      std::vector<Index> _indices;

    };

    template<typename GV>
    class PartitionViewEntitySetIndexSet<GV,Partitions::All>
      : public PartitionViewEntitySetIndexSetBase<GV,Partitions::All>
    {

      using Base = PartitionViewEntitySetIndexSetBase<GV,Dune::Partitions::All>;

      template<typename,typename>
      friend class Dune::PDELab::PartitionViewEntitySet;

    public:

      using typename Base::Index;
      using typename Base::Partitions;
      using typename Base::size_type;
      using typename Base::dim_type;
      using typename Base::CodimMask;

      using Base::baseIndexSet;
      using Base::contains;

    private:

      bool update(bool force)
      {
        if (!Base::update(force))
          return false;
        _mapped_gt_offsets[0] = 0;
        for (const auto& gt : Base::types())
          _mapped_gt_offsets[GlobalGeometryTypeIndex::index(gt) + 1] = baseIndexSet().size(gt);
        std::partial_sum(_mapped_gt_offsets.begin(),_mapped_gt_offsets.end(),_mapped_gt_offsets.begin());
        return true;
      }

    public:

      template<typename E>
      Index index(const E& e) const
      {
        assert(contains(e.type()));
        return baseIndexSet().index(e);
      }

      template<typename E>
      Index uniqueIndex(const E& e) const
      {
        assert(contains(e.type()));
        return baseIndexSet().index(e) + _mapped_gt_offsets[Dune::GlobalGeometryTypeIndex::index(e.type())];
      }

      template<typename E>
      Index subIndex(const E& e, size_type i, dim_type codim) const
      {
#ifndef NDEBUG
        auto gt = ReferenceElements<typename Base::Traits::CoordinateField,GV::dimension>::general(e.type()).type(i,codim);
        assert(contains(gt));
#endif
        return baseIndexSet().subIndex(e,i,codim);
      }

      template<typename E>
      Index uniqueSubIndex(const E& e, size_type i, dim_type codim) const
      {
        auto gt = ReferenceElements<typename Base::Traits::CoordinateField,GV::dimension>::general(e.type()).type(i,codim);
        assert(contains(gt));
        return baseIndexSet().subIndex(e,i,codim) + _mapped_gt_offsets[Dune::GlobalGeometryTypeIndex::index(gt)];
      }

      PartitionViewEntitySetIndexSet(const GV& gv, CodimMask wanted_codims, bool initialize = true)
        : Base(gv,wanted_codims)
      {
        if (initialize)
          update(true);
      }

    private:

      using Base::_mapped_gt_offsets;

    };

    template<typename GV>
    using AllEntitySet = PartitionViewEntitySet<GV,Partitions::All>;

    template<typename GV>
    using OverlappingEntitySet = PartitionViewEntitySet<GV,Partitions::InteriorBorderOverlapFront>;

    template<typename GV>
    using NonOverlappingEntitySet = PartitionViewEntitySet<GV,Partitions::InteriorBorder>;

#ifndef DOXYGEN

    namespace impl {

      template<typename T>
      struct _isEntitySet
      {
        using type = std::false_type;
      };

      template<typename GV,typename P>
      struct _isEntitySet<PartitionViewEntitySet<GV,P>>
      {
        using type = std::true_type;
      };

    }

#endif // DOXYGEN

    //! Type Trait to determine whether T is an EntitySet.
    template<typename T>
    using isEntitySet = typename impl::_isEntitySet<T>::type;

  } // namespace PDELab
} // namespace Dune


#endif // DUNE_PDELAB_COMMON_PARTITIONVIEWENTITYSET_HH
