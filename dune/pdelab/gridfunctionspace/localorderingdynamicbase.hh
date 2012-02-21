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

    protected:
      bool _fixed_size;
      bool _container_blocked;
      bool _max_local_size;

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
