// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_ORDERING_ORDERINGBASE_HH
#define DUNE_PDELAB_ORDERING_ORDERINGBASE_HH

#include <dune/pdelab/common/exceptions.hh>
#include <dune/pdelab/ordering/utility.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspacebase.hh>

#include <vector>

namespace Dune {
  namespace PDELab {

    //! \addtogroup Ordering
    //! \{

    template<typename DI, typename CI>
    class OrderingBase
    {

      template<typename size_type>
      friend struct ::Dune::PDELab::impl::update_ordering_data;

    public:

      typedef OrderingTraits<DI,CI> Traits;

    protected:

      typedef Dune::PDELab::impl::GridFunctionSpaceOrderingData<typename Traits::SizeType> GFSData;

    public:

      typedef HierarchicContainerAllocationTag ContainerAllocationTag;

      typedef DefaultLFSCacheTag CacheTag;

      static const bool has_dynamic_ordering_children = true;

      static const bool consume_tree_index = true;

      typename Traits::ContainerIndex mapIndex(const typename Traits::DOFIndex& di) const
      {
        typename Traits::ContainerIndex ci;
        mapIndex(di.view(),ci);
        return ci;
      }

      void mapIndex(typename Traits::DOFIndexView di, typename Traits::ContainerIndex& ci) const
      {
        if (_delegate)
          _delegate->map_index_dynamic(di,ci);
        else
          {
            _mapIndex(di,ci);
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
        return _child_size_offsets[child_index + 1] - _child_size_offsets[child_index];
      }

      typename Traits::SizeType sizeOffset(const typename Traits::SizeType child_index) const
      {
        return _child_size_offsets[child_index];
      }

      typename Traits::SizeType blockOffset(const typename Traits::SizeType child_index) const
      {
        assert(_merge_mode == MergeMode::lexicographic);
        return _child_block_offsets[child_index];
      }

      typename Traits::SizeType maxLocalSize() const
      {
        return _max_local_size;
      }

      MergeMode::type mergeMode() const
      {
        return _merge_mode;
      }

      void update()
      {
        std::fill(_child_size_offsets.begin(),_child_size_offsets.end(),0);
        std::fill(_child_block_offsets.begin(),_child_block_offsets.end(),0);
        _codim_used.reset();
        _codim_fixed_size.set();
        _max_local_size = 0;
        _block_count = 0;
        typename Traits::SizeType block_carry = 0;
        typename Traits::SizeType size_carry = 0;
        for (typename Traits::SizeType i = 0; i < _child_count; ++i)
          {
            _child_block_offsets[i+1] = (block_carry += _children[i]->blockCount());
            _child_size_offsets[i+1] = (size_carry += _children[i]->size());
            _codim_used |= _children[i]->_codim_used;
            _codim_fixed_size &= _children[i]->_codim_fixed_size;
            _block_count += _children[i]->blockCount();
            _max_local_size += _children[i]->maxLocalSize();
          }
        if (_container_blocked)
          if (_merge_mode == MergeMode::lexicographic)
            {
              _block_count = _child_count;
            }
          else
            {
              if (_child_block_offsets.back() % _child_block_merge_offsets.back() != 0)
                DUNE_THROW(OrderingStructureError,
                           "Invalid ordering structure: "
                           << "total number of blocks ("
                           << _child_block_offsets.back()
                           << ") is not a multiple of the interleaved block size ("
                           << _child_block_merge_offsets.back()
                           << ")."
                           );
              _block_count = _child_block_offsets.back() / _child_block_merge_offsets.back();
            }
        else
          _block_count = _child_block_offsets.back();
        _size = _child_size_offsets.back();
      }

      template<typename Node>
      OrderingBase(Node& node,
                   bool container_blocked,
                   GFSData* gfs_data,
                   VirtualOrderingBase<DI,CI>* delegate = nullptr)
        : _fixed_size(false)
        , _container_blocked(container_blocked)
        , _merge_mode(MergeMode::lexicographic)
        , _child_count(Node::has_dynamic_ordering_children ? TypeTree::degree(node) : 0)
        , _children(_child_count,nullptr)
        , _child_size_offsets((_child_count > 0 ? _child_count + 1 : 0),0)
        , _child_block_offsets((_child_count > 0 ? _child_count + 1 : 0),0)
        , _max_local_size(0)
        , _size(0)
        , _block_count(0)
        , _delegate(delegate)
        , _gfs_data(gfs_data)
      {
        TypeTree::applyToTree(node,extract_child_bases<OrderingBase>(_children));
      }

      template<typename Node>
      OrderingBase(Node& node,
                   bool container_blocked,
                   const std::vector<std::size_t>& merge_offsets,
                   GFSData* gfs_data,
                   VirtualOrderingBase<DI,CI>* delegate = nullptr)
        : _fixed_size(false)
        , _container_blocked(container_blocked)
        , _merge_mode(MergeMode::interleaved)
        , _child_count(Node::has_dynamic_ordering_children ? TypeTree::degree(node) : 0)
        , _children(_child_count,nullptr)
        , _child_size_offsets((_child_count > 0 ? _child_count + 1 : 0),0)
        , _child_block_offsets((_child_count > 0 ? _child_count + 1 : 0),0)
        , _child_block_merge_offsets(merge_offsets)
        , _max_local_size(0)
        , _size(0)
        , _block_count(0)
        , _delegate(delegate)
        , _gfs_data(gfs_data)
      {
        TypeTree::applyToTree(node,extract_child_bases<OrderingBase>(_children));
      }


      bool containerBlocked() const
      {
        return _container_blocked;
      }

      std::size_t childOrderingCount() const
      {
        return _child_count;
      }

      OrderingBase& childOrdering(typename Traits::SizeType i)
      {
        return *_children[i];
      }

      const OrderingBase& childOrdering(typename Traits::SizeType i) const
      {
        return *_children[i];
      }

      bool contains(typename Traits::SizeType codim) const
      {
        return _codim_used.test(codim);
      }

      bool fixedSize(typename Traits::SizeType codim) const
      {
        return _codim_fixed_size.test(codim);
      }

    protected:

      //! Set the delegate called in mapIndex().
      /**
       * When copying an Ordering with a delegate, the derived Ordering
       * *must* call this method with 'this' as its argument in the copy
       * and the move constructors!
       */
      void setDelegate(const VirtualOrderingBase<DI,CI>* delegate)
      {
        _delegate = delegate;
      }

      void _mapIndex(typename Traits::DOFIndexView di, typename Traits::ContainerIndex& ci) const
      {
        typedef typename Traits::SizeType size_type;
        size_type child_index = di.treeIndex().back();
        _children[child_index]->mapIndex(di.back_popped(),ci);
        if (_merge_mode == MergeMode::lexicographic)
          {
            if (_container_blocked)
              ci.push_back(child_index);
            else
              ci.back() += blockOffset(child_index);
          }
        else
          {
            size_type child_block_offset = _child_block_merge_offsets[child_index];
            size_type child_block_size = _child_block_merge_offsets[child_index + 1] - child_block_offset;
            size_type block_index = ci.back() / child_block_size;
            size_type offset = ci.back() % child_block_size;
            if (_container_blocked)
              {
                ci.back() = child_block_offset + offset;
                ci.push_back(block_index);
              }
            else
              {
                size_type block_size = _child_block_merge_offsets.back();
                ci.back() = block_index * block_size + child_block_offset + offset;
              }
          }
      }

    private:

      bool update_gfs_data_size(typename Traits::SizeType& size, typename Traits::SizeType& block_count) const
      {
        size = _size;
        block_count = _block_count;
        return true;
      }

    public:

      bool _fixed_size;
      const bool _container_blocked;
      const MergeMode::type _merge_mode;

      const std::size_t _child_count;
      std::vector<OrderingBase*> _children;

      std::vector<typename Traits::SizeType> _child_size_offsets;
      std::vector<typename Traits::SizeType> _child_block_offsets;
      std::vector<typename Traits::SizeType> _child_block_merge_offsets;
      typename Traits::CodimFlag _codim_used;
      typename Traits::CodimFlag _codim_fixed_size;

      std::size_t _max_local_size;
      std::size_t _size;
      std::size_t _block_count;

      const VirtualOrderingBase<DI,CI>* _delegate;
      GFSData* _gfs_data;

    };

    //! \}

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ORDERING_ORDERINGBASE_HH
