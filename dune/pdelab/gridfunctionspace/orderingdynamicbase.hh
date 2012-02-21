// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_ORDERINGDYNAMICBASE_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_ORDERINGDYNAMICBASE_HH

#include <dune/common/shared_ptr.hh>
#include <dune/pdelab/gridfunctionspace/orderingutility.hh>

#include <vector>

namespace Dune {
  namespace PDELab {

    template<typename DI, typename CI>
    class OrderingBase
    {

    public:

      typedef OrderingTraits<DI,CI> Traits;

      static const bool has_dynamic_ordering_children = true;

      void map_index(typename Traits::DOFIndexView di, typename Traits::ContainerIndex& ci) const
      {
        if (_delegate)
          _delegate->map_index_dynamic(di,ci);
        else
          {
            typename Traits::SizeType child_index = di.treeIndex().back();
            _children[child_index]->map_index(di.back_popped(),ci);
            if (_container_blocked)
              ci.push_back(child_index);
            else
              ci.back() += offset(child_index);
          }
      }

      typename Traits::SizeType size() const
      {
        return _size;
      }

      typename Traits::SizeType blockCount() const
      {
        return _block_count;
      }

      typename Traits::SizeType size(const typename Traits::SizeType child_index) const
      {
        return _child_offsets[child_index + 1] - _child_offsets[child_index];
      }

      typename Traits::SizeType offset(const typename Traits::SizeType child_index) const
      {
        return _child_offsets[child_index];
      }

      typename Traits::SizeType maxLocalSize() const
      {
        return _max_local_size;
      }

      void update()
      {
        std::fill(_child_offsets.begin(),_child_offsets.end(),0);
        typename Traits::SizeType carry = 0;
        _max_local_size = 0;
        _block_count = 0;
        for (typename Traits::SizeType i = 0; i < _child_count; ++i)
          {
            _child_offsets[i+1] = (carry += _children[i]->size());
            _block_count += _children[i]->blockCount();
            _max_local_size += _children[i]->maxLocalSize();
          }
        if (_container_blocked)
          _block_count = _child_count;
        _size = _child_offsets.back();
      }

      template<typename Node>
      OrderingBase(Node& node, bool container_blocked, VirtualOrderingBase<DI,CI>* delegate = nullptr)
        : _container_blocked(container_blocked)
        , _child_count(Node::has_dynamic_ordering_children ? Node::CHILDREN : 0)
        , _children(_child_count,nullptr)
        , _child_offsets(Node::CHILDREN + 1,0)
        , _max_local_size(0)
        , _size(0)
        , _block_count(0)
        , _delegate(delegate)
      {
        TypeTree::applyToTree(node,extract_child_bases<OrderingBase>(_children));
      }

      bool container_blocked() const
      {
        return _container_blocked;
      }

      std::size_t dynamic_child_count() const
      {
        return _child_count;
      }

      OrderingBase& dynamic_child(typename Traits::SizeType i)
      {
        return *_children[i];
      }

      const OrderingBase& dynamic_child(typename Traits::SizeType i) const
      {
        return *_children[i];
      }

    public:

      bool _fixed_size;
      const bool _container_blocked;

      const std::size_t _child_count;
      std::vector<OrderingBase*> _children;

      std::vector<typename Traits::SizeType> _child_offsets;

      std::size_t _max_local_size;
      std::size_t _size;
      std::size_t _block_count;

      const VirtualOrderingBase<DI,CI>* _delegate;

    };


  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_ORDERINGDYNAMICBASE_HH
