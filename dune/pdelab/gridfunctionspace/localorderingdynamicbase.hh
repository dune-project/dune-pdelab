// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_LOCALORDERINGDYNAMICBASE_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_LOCALORDERINGDYNAMICBASE_HH

#include <dune/pdelab/gridfunctionspace/orderingutility.hh>

#include <vector>

namespace Dune {
  namespace PDELab {

    template<typename GV, typename DI, typename CI>
    class LocalOrderingBase
    {

      friend struct collect_a_priori_fixed_size;

      template<typename>
      friend struct update_fixed_size;

      template<typename>
      friend struct post_collect_used_geometry_types;

      template<typename>
      friend struct post_extract_per_entity_sizes;

      friend struct pre_collect_used_geometry_types;

      template<typename>
      friend struct collect_used_geometry_types_from_cell;

      template<typename>
      friend struct extract_per_entity_sizes_from_cell;

      template<typename>
      friend class GridViewOrdering;

    public:

      static const bool has_dynamic_ordering_children = true;

      typedef LocalOrderingTraits<GV,DI,CI> Traits;

      void map_local_index(const typename Traits::SizeType geometry_type_index,
                           const typename Traits::SizeType entity_index,
                           typename Traits::TreeIndexView mi,
                           typename Traits::ContainerIndex& ci) const
      {
        if (_child_count == 0)
          {
            assert(mi.size() == 1 && "MultiIndex length must match GridFunctionSpace tree depth");
            ci.push_back(mi.back());
          }
        else
          {
            const typename Traits::SizeType child_index = mi.back();
            if (!mi.empty())
              _children[child_index]->map_local_index(geometry_type_index,entity_index,mi.back_popped(),ci);
            if (_container_blocked)
              {
                ci.push_back(child_index);
              }
            else if (child_index > 0)
              {
                if (_fixed_size)
                  {
                    const typename Traits::SizeType index = geometry_type_index * _child_count + child_index - 1;
                    ci.back() += _gt_dof_offsets[index];
                  }
                else
                  {
                    assert(_gt_used[geometry_type_index]);
                    const typename Traits::SizeType index = (_gt_entity_offsets[geometry_type_index] + entity_index) * _child_count + child_index - 1;
                    ci.back() += _entity_dof_offsets[index];
                  }
              }
          }
      }


      template<typename ItIn, typename ItOut>
      void map_indices(const ItIn begin, const ItIn end, ItOut out) const
      {
        if (_child_count == 0)
          {
            for (ItIn in = begin; in != end; ++in, ++out)
              {
                assert(in->size() == 1 && "MultiIndex length must match GridFunctionSpace tree depth");
                out->push_back(in->treeIndex().back());
              }
          }
        else if (_container_blocked)
          {
            for (ItIn in = begin; in != end; ++in, ++out)
              out->push_back(in->treeIndex().back());
          }
        else if (_fixed_size)
          {
            for (ItIn in = begin; in != end; ++in, ++out)
              {
                const typename Traits::SizeType child_index = in->treeIndex().back();
                if (child_index > 0)
                  {
                    const typename Traits::SizeType index = in->entityIndex()[0] * _child_count + child_index - 1;
                    out->back() += _gt_dof_offsets[index];
                  }
              }
          }
        else
          {
            for (ItIn in = begin; in != end; ++in, ++out)
              {
                const typename Traits::SizeType child_index = in->treeIndex().back();
                if (child_index > 0)
                  {
                    assert(_gt_used[in->entityIndex()[0]]);
                    const typename Traits::SizeType index = (_gt_entity_offsets[in->entityIndex()[0]] + in->entityIndex()[1]) * _child_count + child_index - 1;
                    out->back() += _entity_dof_offsets[index];
                  }
              }
          }
      }


      typename Traits::SizeType size(const typename Traits::SizeType geometry_type_index, const typename Traits::SizeType entity_index) const
      {
        if (_fixed_size)
          return _gt_dof_offsets[geometry_type_index * _child_count + _child_count - 1];
        else
          return _gt_used[geometry_type_index]
            ? _entity_dof_offsets[(_gt_entity_offsets[geometry_type_index] + entity_index) * _child_count + _child_count - 1]
            : 0;
      }

      typename Traits::SizeType size(const typename Traits::SizeType geometry_type_index, const typename Traits::SizeType entity_index, const typename Traits::SizeType child_index) const
      {
        assert(child_index < _child_count);
        if (_fixed_size)
          {
            const typename Traits::SizeType index = geometry_type_index * _child_count + child_index;
            return child_index > 0 ? _gt_dof_offsets[index] - _gt_dof_offsets[index-1] : _gt_dof_offsets[index];
          }
        else
          {
            if (_gt_used[geometry_type_index])
              {
                const typename Traits::SizeType index = (_gt_entity_offsets[geometry_type_index] + entity_index) * _child_count + child_index;
                return child_index > 0 ? _entity_dof_offsets[index] - _entity_dof_offsets[index-1] : _entity_dof_offsets[index];
              }
            else
              {
                return 0;
              }
          }
      }

      typename Traits::SizeType offset(const typename Traits::SizeType geometry_type_index, const typename Traits::SizeType entity_index, const typename Traits::SizeType child_index) const
      {
        assert(child_index < _child_count);
        assert(_gt_used[geometry_type_index]);
        if (_fixed_size)
          return child_index > 0 ? _gt_dof_offsets[geometry_type_index * _child_count + child_index - 1] : 0;
        else
          return child_index > 0 ? _entity_dof_offsets[(_gt_entity_offsets[geometry_type_index] + entity_index) * _child_count + child_index - 1] : 0;
      }

      template<typename Node>
      LocalOrderingBase(Node& node, bool container_blocked)
        : _container_blocked(container_blocked)
        , _child_count(Node::CHILDREN)
        , _children(Node::CHILDREN,nullptr)
      {
        TypeTree::applyToTree(node,extract_child_bases<LocalOrderingBase>(_children));
      }

      bool fixedSize() const
      {
        return _fixed_size;
      }

      bool contains(const GeometryType& gt) const
      {
        return _gt_used[GlobalGeometryTypeIndex::index(gt)];
      }

      bool contains(typename Traits::SizeType codim) const
      {
        return _codim_used[codim];
      }

      typename Traits::SizeType maxLocalSize() const
      {
        return _max_local_size;
      }

    protected:

      LocalOrderingBase& dynamic_child(typename Traits::SizeType i)
      {
        return *_children[i];
      }

      const LocalOrderingBase& dynamic_child(typename Traits::SizeType i) const
      {
        return *_children[i];
      }

      void disable_container_blocking()
      {
        _container_blocked = false;
      }


      bool _fixed_size;
      bool _fixed_size_possible;
      bool _container_blocked;
      std::size_t _max_local_size;

      const std::size_t _child_count;
      std::vector<LocalOrderingBase*> _children;

      std::vector<bool> _codim_used;
      std::vector<bool> _gt_used;

      std::vector<typename Traits::SizeType> _gt_entity_offsets;
      std::vector<typename Traits::SizeType> _gt_dof_offsets;
      std::vector<typename Traits::SizeType> _entity_dof_offsets;

    };


  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_LOCALORDERINGDYNAMICBASE_HH
