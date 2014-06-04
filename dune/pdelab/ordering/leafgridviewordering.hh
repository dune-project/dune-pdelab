// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_ORDERING_LEAFGRIDVIEWORDERING_HH
#define DUNE_PDELAB_ORDERING_LEAFGRIDVIEWORDERING_HH

#include <dune/pdelab/ordering/directleaflocalordering.hh>
#include <dune/pdelab/ordering/leaforderingbase.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup Ordering
    //! \{

    //! Gridview ordering for leaf spaces
    template<typename LocalOrdering>
    class LeafGridViewOrdering
      : public LeafOrderingBase<LocalOrdering>
    {
    public:
      typedef typename LocalOrdering::Traits Traits;

    private:

      typedef typename Traits::GridView GV;

      typedef LeafOrderingBase<LocalOrdering> BaseT;
      typedef typename BaseT::NodeT NodeT;

    public:

      LeafGridViewOrdering(const typename NodeT::NodeStorage& local_ordering, bool container_blocked, typename BaseT::GFSData* gfs_data)
        : BaseT(local_ordering, container_blocked, gfs_data)
        , _gv(this->template child<0>().gridView())
      {}

#ifndef DOXYGEN

// we need to override the default copy / move ctor to fix the delegate pointer, but that is
// hardly interesting to our users...

      LeafGridViewOrdering(const LeafGridViewOrdering& r)
        : BaseT(r)
        , _gv(r._gv)
      {}

#if HAVE_RVALUE_REFERENCES

      LeafGridViewOrdering(LeafGridViewOrdering&& r)
        : BaseT(std::move(r))
        , _gv(r._gv)
      {}

#endif // HAVE_RVALUE_REFERENCES

#endif // DOXYGEN

      virtual void update()
      {
        LocalOrdering& lo = this->localOrdering();
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
            _codim_fixed_size.set();
          }
        else
          {
            _block_count = _size = lo._entity_dof_offsets.back();
            _codim_fixed_size.reset();
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
      using BaseT::_gt_dof_offsets;

      typename Traits::GridView _gv;
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
        return transformed_type(make_tuple(make_shared<LocalOrdering>(gfs.finiteElementMapStorage(),gfs.gridView())),gfs.backend().blocked(gfs),const_cast<GFS*>(&gfs));
      }

      static transformed_storage_type transform_storage(shared_ptr<const GFS> gfs, const Transformation& t)
      {
        return make_shared<transformed_type>(make_tuple(make_shared<LocalOrdering>(gfs->finiteElementMapStorage(),gfs->gridView())),gfs->backend().blocked(*gfs),const_cast<GFS*>(gfs.get()));
      }

    };


    template<typename GFS, typename Transformation, typename Params>
    direct_leaf_gfs_to_gridview_ordering_descriptor<GFS,Transformation>
    register_leaf_gfs_to_ordering_descriptor(GFS*,Transformation*,LeafOrderingTag<Params>*);

   //! \} group Ordering
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ORDERING_LEAFORDERING_HH
