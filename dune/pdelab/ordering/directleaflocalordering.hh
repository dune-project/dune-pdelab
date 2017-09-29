// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_ORDERING_DIRECTLEAFLOCALORDERING_HH
#define DUNE_PDELAB_ORDERING_DIRECTLEAFLOCALORDERING_HH

#include <dune/typetree/leafnode.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/localfunctions/common/interfaceswitch.hh>
#include <dune/localfunctions/common/localkey.hh>
#include <dune/pdelab/ordering/utility.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspacebase.hh>

#include <vector>
#include <numeric>

namespace Dune {
  namespace PDELab {

    //! \addtogroup Ordering
    //! \{

    template<typename OrderingTag, typename FEM, typename ES, typename DI, typename CI>
    class DirectLeafLocalOrdering
      : public TypeTree::LeafNode
    {

      template<typename>
      friend class LeafGridViewOrdering;

      template<typename>
      friend class LeafOrderingBase;

      template<typename size_type>
      friend struct ::Dune::PDELab::impl::update_ordering_data;

    public:

      typedef LocalOrderingTraits<ES,DI,CI> Traits;

    private:

      typedef impl::GridFunctionSpaceOrderingData<typename Traits::SizeType> GFSData;

    public:

      void map_local_index(const typename Traits::SizeType geometry_type_index,
                           const typename Traits::SizeType entity_index,
                           typename Traits::TreeIndexView mi,
                           typename Traits::ContainerIndex& ci) const
      {
        DUNE_THROW(NotImplemented,"not implemented");
      }

      template<typename ItIn, typename ItOut>
      void map_lfs_indices(const ItIn begin, const ItIn end, ItOut out) const
      {
        // don't do anything - this is handled by the specialized GridViewOrdering
      }

      template<typename CIOutIterator, typename DIOutIterator = DummyDOFIndexIterator>
      typename Traits::SizeType
      extract_entity_indices(const typename Traits::DOFIndex::EntityIndex& ei,
                             typename Traits::SizeType child_index,
                             CIOutIterator ci_out, const CIOutIterator ci_end,
                             DIOutIterator di_out = DIOutIterator()) const
      {
        const typename Traits::SizeType s = size(ei);

        // Handle DOF indices
        for (typename Traits::SizeType i = 0; i < s; ++i, ++di_out)
          di_out->treeIndex().push_back(i);

        // only return the size, as the tree visitor expects that from all leaf nodes.
        // The actual index processing is done by the specialized GridViewOrdering.
        return s;
      }

      typename Traits::SizeType size(const typename Traits::DOFIndex::EntityIndex& index) const
      {
        return size(
          Traits::DOFIndexAccessor::GeometryIndex::geometryType(index),
          Traits::DOFIndexAccessor::GeometryIndex::entityIndex(index)
        );
      }

      typename Traits::SizeType size(const typename Traits::SizeType geometry_type_index, const typename Traits::SizeType entity_index) const
      {
        typedef typename Traits::SizeType size_type;
        if (_fixed_size)
          return _gt_dof_sizes[geometry_type_index];
        else if (_gt_used[geometry_type_index])
          {
            const size_type index = _gt_entity_offsets[geometry_type_index] + entity_index;
            return _entity_dof_offsets[index+1] - _entity_dof_offsets[index];
          }
        else
          return 0;
      }

      typename Traits::SizeType size(const typename Traits::SizeType geometry_type_index, const typename Traits::SizeType entity_index, const typename Traits::SizeType child_index) const
      {
        DUNE_THROW(NotImplemented,"not implemented");
      }

      typename Traits::SizeType offset(const typename Traits::SizeType geometry_type_index, const typename Traits::SizeType entity_index, const typename Traits::SizeType child_index) const
      {
        assert(child_index == 0);
        return 0;
      }

      DirectLeafLocalOrdering(const shared_ptr<const FEM>& fem, const ES& es)
        : _fem(fem)
        , _es(es)
        , _fixed_size(false)
        , _container_blocked(false)
        , _gfs_data(nullptr)
      {}

      const typename Traits::EntitySet& entitySet() const
      {
        return _es;
      }

      const FEM& finiteElementMap() const
      {
        return *_fem;
      }

    private:

      typedef FiniteElementInterfaceSwitch<
      typename FEM::Traits::FiniteElement
      > FESwitch;


      void update_a_priori_fixed_size()
      {
        _fixed_size = _fem->fixedSize();
      }

      template<typename CodimMask>
      void collect_used_codims(CodimMask& codims) const
      {
        for (typename ES::dim_type codim = 0; codim <= ES::dimension; ++codim)
          if (_fem->hasDOFs(codim))
            codims.set(codim);
      }

      template<typename It>
      void update_fixed_size(It it, const It end)
      {
        assert(_fixed_size);

        _max_local_size = _fem->maxLocalSize();

        typedef typename Traits::SizeType size_type;
        const size_type dim = Traits::GridView::dimension;
        _codim_used.reset();
        _gt_used.assign(GlobalGeometryTypeIndex::size(dim),false);
        _gt_dof_sizes.assign(GlobalGeometryTypeIndex::size(dim),0);
        for (; it != end; ++it)
          {
            size_type size = _fem->size(*it);
            _gt_dof_sizes[GlobalGeometryTypeIndex::index(*it)] = size;
            _gt_used[GlobalGeometryTypeIndex::index(*it)] = size > 0;
            _codim_used[dim - it->dim()] = _codim_used[dim - it->dim()] || (size > 0);
          }

        _codim_fixed_size.set();
      }


      void pre_collect_used_geometry_types_from_cell()
      {
        typedef typename Traits::SizeType size_type;
        const size_type dim = Traits::GridView::dimension;

        _codim_used.reset();
        _gt_used.assign(GlobalGeometryTypeIndex::size(dim),false);
        _gt_dof_sizes.assign(GlobalGeometryTypeIndex::size(dim),0);
        _local_gt_dof_sizes.resize(GlobalGeometryTypeIndex::size(dim));
        _max_local_size = 0;
        _fixed_size_possible = true;
      }


      void collect_used_geometry_types_from_cell(const typename Traits::GridView::template Codim<0>::Entity& cell)
      {
        FESwitch::setStore(_fe_store,_fem->find(cell));

        const typename FESwitch::Coefficients& coeffs =
          FESwitch::coefficients(*_fe_store);

        _max_local_size = std::max(_max_local_size,coeffs.size());

        auto ref_el = ReferenceElements<typename Traits::GridView::ctype,Traits::GridView::dimension>::general(cell.type());

        for (std::size_t i = 0; i < coeffs.size(); ++i)
          {
            const LocalKey& key = coeffs.localKey(i);
            GeometryType gt = ref_el.type(key.subEntity(),key.codim());
            _gt_used[GlobalGeometryTypeIndex::index(gt)] = true;
            _codim_used.set(key.codim());
          }
      }


      template<typename It>
      void allocate_entity_offset_vector(It it, const It end)
      {
        _gt_entity_offsets.assign(GlobalGeometryTypeIndex::size(ES::dimension) + 1,0);
        for (; it != end; ++it)
          {
            if (_gt_used[GlobalGeometryTypeIndex::index(*it)])
              _gt_entity_offsets[GlobalGeometryTypeIndex::index(*it) + 1] = _es.indexSet().size(*it);
          }
        std::partial_sum(_gt_entity_offsets.begin(),_gt_entity_offsets.end(),_gt_entity_offsets.begin());
        _entity_dof_offsets.assign(_gt_entity_offsets.back() + 1,0);

        // Don't claim fixed size for any codim for now
        _codim_fixed_size.reset();
      }


      void extract_per_entity_sizes_from_cell(const typename Traits::GridView::template Codim<0>::Entity& cell)
      {
        if (this->_fixed_size_possible)
          std::fill(_local_gt_dof_sizes.begin(),_local_gt_dof_sizes.end(),0);

        FESwitch::setStore(_fe_store,_fem->find(cell));

        const typename FESwitch::Coefficients& coeffs =
          FESwitch::coefficients(*_fe_store);

        typedef typename Traits::SizeType size_type;

        auto ref_el = ReferenceElements<typename Traits::GridView::ctype,Traits::GridView::dimension>::general(cell.type());

        for (std::size_t i = 0; i < coeffs.size(); ++i)
          {
            const LocalKey& key = coeffs.localKey(i);
            GeometryType gt = ref_el.type(key.subEntity(),key.codim());
            const size_type geometry_type_index = GlobalGeometryTypeIndex::index(gt);

            const size_type entity_index = _es.indexSet().subIndex(cell,key.subEntity(),key.codim());
            const size_type index = _gt_entity_offsets[geometry_type_index] + entity_index;
            _local_gt_dof_sizes[geometry_type_index] = _entity_dof_offsets[index+1] = std::max(_entity_dof_offsets[index+1],static_cast<size_type>(key.index() + 1));
          }

        if (_fixed_size_possible)
          {
            for (size_type i = 0; i < _local_gt_dof_sizes.size(); ++i)
              {
                if (_local_gt_dof_sizes[i] > 0)
                  {
                    if (_gt_dof_sizes[i] == 0)
                      _gt_dof_sizes[i] = _local_gt_dof_sizes[i];
                    else if (_gt_dof_sizes[i] != _local_gt_dof_sizes[i])
                      {
                        _fixed_size_possible = false;
                        break;
                      }
                  }
                else
                  {
                    // catch special case where some grid entities sometimes don't have any DOFs associated with
                    // them, but are const size if there are any DOFs, i.e. a pattern of #(DOFs) == 0 || #(DOFs) == c > 0
                    _fixed_size_possible = _fixed_size_possible && (_gt_dof_sizes[i] == 0);
                  }
              }
          }
      }


      void finalize_non_fixed_size_update()
      {
        if (_fixed_size_possible)
          {
            // free per-entity offsets
            _entity_dof_offsets = std::vector<typename Traits::SizeType>();
            _fixed_size = true;
            _codim_fixed_size.set();
          }
        else
          {
            // convert per-entity sizes to offsets
            std::partial_sum(_entity_dof_offsets.begin(),_entity_dof_offsets.end(),_entity_dof_offsets.begin());
            _fixed_size = false;
            _codim_fixed_size.reset();
          }
      }


      typename Traits::SizeType maxLocalSize() const
      {
        return _max_local_size;
      }

    private:

      bool update_gfs_data_size(typename Traits::SizeType& size, typename Traits::SizeType& block_count) const
      {
        return false;
      }

    protected:

      shared_ptr<const FEM> _fem;
      typename FESwitch::Store _fe_store;

      ES _es;
      bool _fixed_size;
      bool _fixed_size_possible;
      typename Traits::SizeType _max_local_size;
      const bool _container_blocked;

      typename Traits::CodimFlag _codim_used;
      typename Traits::CodimFlag _codim_fixed_size;
      std::vector<bool> _gt_used;

      std::vector<typename Traits::SizeType> _gt_entity_offsets;
      std::vector<typename Traits::SizeType> _gt_dof_sizes;
      std::vector<typename Traits::SizeType> _entity_dof_offsets;
      std::vector<typename Traits::SizeType> _local_gt_dof_sizes;

      // This is only here to make the visitor happy that traverses all
      // Orderings to manipulate the contained GFSData
      GFSData* _gfs_data;

    };

    //! \} group ordering

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ORDERING_DIRECTLEAFLOCALORDERING_HH
