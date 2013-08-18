// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_ORDERING_LEAFGRIDVIEWORDERING_HH
#define DUNE_PDELAB_ORDERING_LEAFGRIDVIEWORDERING_HH

#include <dune/pdelab/common/typetree.hh>
#include <dune/pdelab/ordering/utility.hh>
#include <dune/pdelab/ordering/orderingbase.hh>
#include <dune/pdelab/ordering/directleaflocalordering.hh>

namespace Dune {
  namespace PDELab {



    //! \addtogroup GridFunctionSpace
    //! \ingroup PDELab
    //! \{

    //! Gridview ordering for leaf spaces
    template<typename LocalOrdering>
    class LeafGridViewOrdering
      : public TypeTree::VariadicCompositeNode<LocalOrdering>
      , public VirtualOrderingBase<typename LocalOrdering::Traits::DOFIndex,
                                   typename LocalOrdering::Traits::GlobalDOFIndex,
                                   typename LocalOrdering::Traits::ContainerIndex>
      , public OrderingBase<typename LocalOrdering::Traits::DOFIndex,
                            typename LocalOrdering::Traits::GlobalDOFIndex,
                            typename LocalOrdering::Traits::ContainerIndex>
    {
    public:
      typedef typename LocalOrdering::Traits Traits;

      static const bool has_dynamic_ordering_children = false;

      static const bool consume_tree_index = false;

    private:

      typedef typename Traits::GridView GV;

      typedef TypeTree::VariadicCompositeNode<LocalOrdering> NodeT;

      typedef OrderingBase<typename LocalOrdering::Traits::DOFIndex,
                           typename LocalOrdering::Traits::GlobalDOFIndex,
                           typename LocalOrdering::Traits::ContainerIndex> BaseT;

    public:

      LocalOrdering& localOrdering()
      {
        return this->template child<0>();
      }

      const LocalOrdering& localOrdering() const
      {
        return this->template child<0>();
      }


      LeafGridViewOrdering(const typename NodeT::NodeStorage& local_ordering, bool container_blocked, typename BaseT::GFSData* gfs_data)
        : NodeT(local_ordering)
        , BaseT(*this,container_blocked,gfs_data,this)
        , _gv(this->template child<0>().gridView())
      {
        // copy grid partition information from local ordering.
        this->setPartitionSet(localOrdering());
      }

#ifndef DOXYGEN

// we need to override the default copy / move ctor to fix the delegate pointer, but that is
// hardly interesting to our users...

      LeafGridViewOrdering(const LeafGridViewOrdering& r)
        : NodeT(r.nodeStorage())
        , BaseT(r)
        , _gv(r._gv)
        , _gt_dof_offsets(r._gt_dof_offsets)
      {
        this->setDelegate(this);
      }

#if HAVE_RVALUE_REFERENCES

      LeafGridViewOrdering(LeafGridViewOrdering&& r)
        : NodeT(r.nodeStorage())
        , BaseT(std::move(r))
        , _gv(r._gv)
        , _gt_dof_offsets(std::move(r._gt_dof_offsets))
      {
        this->setDelegate(this);
      }

#endif // HAVE_RVALUE_REFERENCES

#endif // DOXYGEN

      virtual void map_index_dynamic(typename Traits::DOFIndexView di, typename Traits::ContainerIndex& ci) const
      {
        mapIndex(di,ci);
      }

      typename Traits::ContainerIndex mapIndex(const typename Traits::DOFIndex& di) const
      {
        typename Traits::ContainerIndex ci;
        mapIndex(di.view(),ci);
        return ci;
      }

      void mapIndex(typename Traits::DOFIndexView di, typename Traits::ContainerIndex& ci) const
      {

        const typename Traits::SizeType geometry_type_index = Traits::DOFIndexAccessor::geometryType(di);
        const typename Traits::SizeType entity_index = Traits::DOFIndexAccessor::entityIndex(di);
        assert (di.treeIndex().size() == 1);
        ci.push_back(di.treeIndex().back());

        if (localOrdering()._fixed_size)
          {
            if (_container_blocked)
              {
                const typename Traits::SizeType geometry_type_index = Traits::DOFIndexAccessor::geometryType(di);
                const typename Traits::SizeType entity_index = Traits::DOFIndexAccessor::entityIndex(di);
                ci.push_back(_gt_dof_offsets[geometry_type_index] + entity_index);
              }
            else
              {
                ci.back() += _gt_dof_offsets[geometry_type_index] + entity_index * localOrdering()._gt_dof_sizes[geometry_type_index];
              }
          }
        else
          {
            if (_container_blocked)
              {
                // This check is needed to avoid a horrid stream of compiler warnings about
                // exceeding array bounds in ReservedVector!
                if (ci.size() < ci.capacity())
                  ci.push_back(localOrdering()._gt_entity_offsets[geometry_type_index] + entity_index);
                else
                  {
                    DUNE_THROW(Dune::Exception,"Container blocking incompatible with backend structure");
                  }
              }
            else
              {
                ci.back() += localOrdering()._entity_dof_offsets[localOrdering()._gt_entity_offsets[geometry_type_index] + entity_index];
              }
          }
      }


      template<typename ItIn, typename ItOut>
      void map_lfs_indices(const ItIn begin, const ItIn end, ItOut out) const
      {
        typedef typename Traits::SizeType size_type;

        if (localOrdering()._fixed_size)
          {
            if (_container_blocked)
              {
                for (ItIn in = begin; in != end; ++in, ++out)
                  {
                    assert(in->treeIndex().size() == 1);
                    out->push_back(in->treeIndex().back());
                    const size_type geometry_type_index = Traits::DOFIndexAccessor::geometryType(*in);
                    const size_type entity_index = Traits::DOFIndexAccessor::entityIndex(*in);
                    out->push_back(_gt_dof_offsets[geometry_type_index] + entity_index);
                  }
              }
            else
              {
                for (ItIn in = begin; in != end; ++in, ++out)
                  {
                    assert(in->treeIndex().size() == 1);
                    out->push_back(in->treeIndex().back());
                    const size_type geometry_type_index = Traits::DOFIndexAccessor::geometryType(*in);
                    const size_type entity_index = Traits::DOFIndexAccessor::entityIndex(*in);
                    out->back() += _gt_dof_offsets[geometry_type_index] + entity_index * localOrdering()._gt_dof_sizes[geometry_type_index];
                  }
              }
          }
        else
          {
            if (_container_blocked)
              {
                for (ItIn in = begin; in != end; ++in, ++out)
                  {
                    assert(in->treeIndex().size() == 1);
                    out->push_back(in->treeIndex().back());
                    const size_type geometry_type_index = Traits::DOFIndexAccessor::geometryType(*in);
                    const size_type entity_index = Traits::DOFIndexAccessor::entityIndex(*in);
                    out->push_back(localOrdering()._gt_entity_offsets[geometry_type_index] + entity_index);
                  }
              }
            else
              {
                for (ItIn in = begin; in != end; ++in, ++out)
                  {
                    assert(in->treeIndex().size() == 1);
                    out->push_back(in->treeIndex().back());
                    const size_type geometry_type_index = Traits::DOFIndexAccessor::geometryType(*in);
                    const size_type entity_index = Traits::DOFIndexAccessor::entityIndex(*in);
                    out->back() += localOrdering()._entity_dof_offsets[localOrdering()._gt_entity_offsets[geometry_type_index] + entity_index];
                  }
              }
          }
      }

      template<typename CIOutIterator>
      typename Traits::SizeType
      extract_entity_indices(const typename Traits::DOFIndex::EntityIndex& ei,
                             typename Traits::SizeType child_index,
                             CIOutIterator ci_out, const CIOutIterator ci_end) const
      {
        typedef typename Traits::SizeType size_type;

        const size_type geometry_type_index = Traits::DOFIndexAccessor::GeometryIndex::geometryType(ei);
        const size_type entity_index = Traits::DOFIndexAccessor::GeometryIndex::entityIndex(ei);

        if (!localOrdering()._gt_used[geometry_type_index])
          return 0;

        if (localOrdering()._fixed_size)
          {
            size_type size = localOrdering()._gt_dof_sizes[geometry_type_index];
            if (_container_blocked)
              {
                for (size_type i = 0; i < size; ++i, ++ci_out)
                  {
                    ci_out->push_back(i);
                    ci_out->push_back(_gt_dof_offsets[geometry_type_index] + entity_index);
                  }
              }
            else
              {
                for (size_type i = 0; i < size; ++i, ++ci_out)
                  {
                    ci_out->push_back(i);
                    ci_out->back() += _gt_dof_offsets[geometry_type_index] + entity_index * localOrdering()._gt_dof_sizes[geometry_type_index];
                  }
              }
            return 0;
          }
        else
          {
            size_type index = localOrdering()._gt_entity_offsets[geometry_type_index] + entity_index;
            size_type size = localOrdering()._entity_dof_offsets[index+1] - localOrdering()._entity_dof_offsets[index];
            if (_container_blocked)
              {
                for (size_type i = 0; i < size; ++i, ++ci_out)
                  {
                    ci_out->push_back(i);
                    ci_out->push_back(index);
                  }
              }
            else
              {
                for (size_type i = 0; i < size; ++i, ++ci_out)
                  {
                    ci_out->push_back(i);
                    ci_out->back() += localOrdering()._entity_dof_offsets[index];
                  }
              }
            return 0;
          }
      }

      void update()
      {
        LocalOrdering& lo = localOrdering();
        lo.update_a_priori_fixed_size();

        const std::size_t dim = GV::dimension;

        typedef typename Traits::SizeType size_type;
        typedef std::vector<GeometryType> GTVector;
        GTVector geom_types;

        for (size_type cc = 0; cc <= dim; ++cc)
          {
            const GTVector& per_codim_geom_types = _gv.indexSet().geomTypes(cc);
            std::copy(per_codim_geom_types.begin(),per_codim_geom_types.end(),std::back_inserter(geom_types));
          }

        if (lo._fixed_size)
          {
            lo.update_fixed_size(geom_types.begin(),geom_types.end());
          }
        else
          {
            lo.pre_collect_used_geometry_types_from_cell();

            typedef typename GV::template Codim<0>::Iterator CellIterator;
            typedef typename GV::template Codim<0>::Entity Cell;

            const CellIterator end_it = _gv.template end<0>();
            for (CellIterator it = _gv.template begin<0>(); it != end_it; ++it)
              {
                lo.collect_used_geometry_types_from_cell(*it);
              }

            lo.allocate_entity_offset_vector(geom_types.begin(),geom_types.end());

            for (CellIterator it = _gv.template begin<0>(); it != end_it; ++it)
              {
                lo.extract_per_entity_sizes_from_cell(*it);
              }

            // FIXME: handling of blocked containers!
            lo.finalize_non_fixed_size_update();
          }

        // we need to re-test here, as the local ordering could have detected a fixed size ordering
        // and switched its implementation
        if (lo._fixed_size)
          {
            _gt_dof_offsets.assign(GlobalGeometryTypeIndex::size(dim) + 1,0);
            _size = 0;

            const GTVector::const_iterator end_it = geom_types.end();
            for (GTVector::const_iterator it = geom_types.begin(); it != end_it; ++it)
              {
                const size_type gt_index = GlobalGeometryTypeIndex::index(*it);
                size_type gt_size = lo.size(gt_index,0);
                size_type entity_count = _gv.indexSet().size(*it);
                _size += gt_size * entity_count;
                if (_container_blocked)
                  gt_size = gt_size > 0;
                _gt_dof_offsets[gt_index + 1] = gt_size * entity_count;
              }
            std::partial_sum(_gt_dof_offsets.begin(),_gt_dof_offsets.end(),_gt_dof_offsets.begin());
            _block_count = _gt_dof_offsets.back();
          }
        else
          {
            _block_count = _size = lo._entity_dof_offsets.back();
          }

        _fixed_size = lo._fixed_size;
        _max_local_size = lo.maxLocalSize();

        _codim_used = lo._codim_used;
        _codim_fixed_size = lo._codim_fixed_size;

      }

      using BaseT::fixedSize;

    private:

      using BaseT::_max_local_size;
      using BaseT::_size;
      using BaseT::_block_count;
      using BaseT::_container_blocked;
      using BaseT::_fixed_size;
      using BaseT::_codim_used;
      using BaseT::_codim_fixed_size;

      typename Traits::GridView _gv;
      std::vector<typename Traits::SizeType> _gt_dof_offsets;

    };


    template<typename GFS, typename Transformation>
    struct direct_leaf_gfs_to_gridview_ordering_descriptor
    {

      static const bool recursive = false;

      typedef DirectLeafLocalOrdering<typename GFS::Traits::OrderingTag,
                                      typename GFS::Traits::FiniteElementMap,
                                      typename GFS::Traits::GridView,
                                      typename Transformation::DOFIndex,
                                      typename Transformation::ContainerIndex
                                      > LocalOrdering;

      typedef LeafGridViewOrdering<LocalOrdering> GridViewOrdering;

      typedef GridViewOrdering transformed_type;
      typedef shared_ptr<transformed_type> transformed_storage_type;

      static transformed_type transform(const GFS& gfs, const Transformation& t)
      {
        return transformed_type(make_tuple(make_shared<LocalOrdering>(gfs.finiteElementMapStorage(),gfs.gridView())),gfs.backend().blocked(),const_cast<GFS*>(&gfs));
      }

      static transformed_storage_type transform_storage(shared_ptr<const GFS> gfs, const Transformation& t)
      {
        return make_shared<transformed_type>(make_tuple(make_shared<LocalOrdering>(gfs->finiteElementMapStorage(),gfs->gridView())),gfs->backend().blocked(),const_cast<GFS*>(gfs.get()));
      }

    };


    template<typename GFS, typename Transformation, typename Params>
    direct_leaf_gfs_to_gridview_ordering_descriptor<GFS,Transformation>
    register_leaf_gfs_to_ordering_descriptor(GFS*,Transformation*,LeafOrderingTag<Params>*);


   //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ORDERING_LEAFORDERING_HH
