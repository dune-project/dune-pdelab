// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_ORDERINGDYNAMICBASE_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_ORDERINGDYNAMICBASE_HH

#include <dune/common/shared_ptr.hh>
#include <dune/pdelab/gridfunctionspace/orderingutility.hh>

#include <vector>

namespace Dune {
  namespace PDELab {

    template<typename MI, typename CI>
    class OrderingBase
    {

    public:

      typedef OrderingTraits<MI,CI> Traits;

      void map_index(const typename Traits::MultiIndex& mi, typename Traits::ContainerIndex& ci) const
      {
        if (_delegate)
          _delegate->map_index_dynamic(mi,ci);
        else
          {
            typename Traits::SizeType child_index = mi.back();
            mi.pop_back();
            _children[child_index]->map_index(mi,ci);
            mi.push_back(child_index);
            if (_container_blocked)
              ci.push_back(child_index);
            else
              ci.back() += offset(child_index);
          }
      }

      typename Traits::SizeType size() const
      {
        return _child_offsets.back();
      }

      typename Traits::SizeType size(const typename Traits::SizeType child_index) const
      {
        return _child_offsets[child_index + 1] - _child_offsets[child_index];
      }

      typename Traits::SizeType offset(const typename Traits::SizeType child_index) const
      {
        return _child_offsets[child_index];
      }

      void update()
      {
        typename Traits::SizeType carry = 0;
        for (typename Traits::SizeType i = 0; i < _child_count; ++i)
          _child_offsets[i+1] = (carry += _children[i].size());
      }

      template<typename Node>
      OrderingBase(Node& node, VirtualOrderingBase<MI,CI>* delegate = nullptr)
        : _container_blocked(node.container_blocked())
        , _child_count(Node::CHILDREN)
        , _children(Node::CHILDREN,nullptr)
        , _child_offsets(Node::CHILDREN + 1,0)
        , _delegate(delegate)
      {
        TypeTree::applyToTree(node,extract_child_bases<OrderingBase>(_children));
      }

    protected:

      bool _fixed_size;
      const bool _container_blocked;

      const std::size_t _child_count;
      std::vector<OrderingBase*> _children;

      std::vector<SizeType> _child_offsets;

      const VirtualOrderingBase<MI,CI>* _delegate;

    };


  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_ORDERINGDYNAMICBASE_HH
