// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_LOCALORDERINGDYNAMICBASE_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_LOCALORDERINGDYNAMICBASE_HH

#include <dune/pdelab/gridfunctionspace/orderingutility.hh>

namespace Dune {
  namespace PDELab {

    template<typename MI, typename CI>
    class LocalOrderingBase
    {

    public:

      typedef OrderingTraits<MI,CI> Traits;

      void map_local_index(const typename Traits::SizeType geometry_type_index,
                           const typename Traits::SizeType entity_index,
                           const typename Traits::MultiIndex& mi,
                           typename Traits::ContainerIndex& ci) const
      {
        if (_child_count == 0)
          {
            assert(mi.size() == 1 && "MultiIndex length must match tree depth");
            ci.push_back(mi.back());
          }
        else
          {
            const typename Traits::SizeType child_index = mi.back();
            mi.pop_back();
            _children[child_index]->map_local_index(geometry_type_index,entity_index,mi,ci);
            mi.push_back(child_index);
            if (_container_blocked)
              {
                ci.push_back(child_index);
              }
            else if (_fixed_size)
              {
                ci.back() += _gt_dof_offsets[geometry_type_index * (_child_count + 1) + child_index];
              }
            else
              {
                ci.back() += _entity_dof_offsets[(_gt_entity_offsets[geometry_type_index] + entity_index) * (_child_count + 1) + child_index];
              }
          }
      }

      typename Traits::SizeType size(const typename Traits::SizeType geometry_type_index, const typename Traits::SizeType entity_index) const
      {
        if (_child_count == 0)
          return size(geometry_type_index,entity_index,0);
        else
          return offset(geometry_type_index,entity_index,_child_count);
      }

      typename Traits::SizeType size(const typename Traits::SizeType geometry_type_index, const typename Traits::SizeType entity_index, const typename Traits::SizeType child_index) const
      {
        if (_child_count == 0)
          {
            assert (child_index == 0);
            if (_fixed_size)
              {
                return _gt_dof_offsets[geometry_type_index];
              }
            else
              {
                return _entity_dof_offsets[_gt_entity_offsets[geometry_type_index + entity_index]];
              }
          }
        else
          {
            if (_fixed_size)
              {
                const typename Traits::SizeType index = geometry_type_index * (_child_count + 1) + child_index;
                return _gt_dof_offsets[index+1] - gt_dof_offsets[index];
              }
            else
              {
                const typename Traits::SizeType index = (_gt_entity_offsets[geometry_type_index] + entity_index) * (_child_count + 1) + child_index;
                return _entity_dof_offsets[index+1] - entity_dof_offsets[index];
              }
          }
      }

      typename Traits::SizeType offset(const typename Traits::SizeType geometry_type_index, const typename Traits::SizeType entity_index, const typename Traits::SizeType child_index) const
      {
        if (_child_count == 0)
          {
            assert(child_index < 2);
            return child_index == 0 ? 0 : size(geometry_type_index,entity_index,0);
          }
        else
          {
            if (_fixed_size)
              {
                return _gt_dof_offsets[geometry_type_index * (_child_count + 1) + child_index];
              }
            else
              {
                return _entity_dof_offsets[(_gt_entity_offsets[geometry_type_index] + entity_index) * (_child_count + 1) + child_index];
              }
          }
      }

      void update()
      {
      }

      template<typename Node>
      LocalOrderingBase(Node& node)
        : _container_blocked(node.container_blocked())
        , _child_count(Node::CHILDREN)
        , _children(Node::CHILDREN,nullptr)
      {
        TypeTree::applyToTree(node,extract_child_bases<LocalOrderingBase>(_children));
      }

    protected:
      bool _fixed_size;
      const bool _container_blocked;

      const std::size_t _child_count;
      std::vector<LocalOrderingBase*> _children;

      std::vector<SizeType> _gt_entity_offsets;
      std::vector<SizeType> _gt_dof_offsets;
      std::vector<SizeType> _entity_dof_offsets;

    };


  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_LOCALORDERINGDYNAMICBASE_HH
