// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_ORDERING_LOCALORDERINGBASE_HH
#define DUNE_PDELAB_ORDERING_LOCALORDERINGBASE_HH

#include <dune/pdelab/ordering/utility.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspacebase.hh>

#include <dune/common/rangeutilities.hh>

#include <vector>

namespace Dune {
  namespace PDELab {

    //! \addtogroup Ordering
    //! \{

    /**
     * @brief Entity-wise orderings.
     * @details A tree structure of this class is able to provide
     * (multi-index) indexations for each entity in an entity set.
     * @tparam ES Entity Set
     * @tparam DI DOFIndex
     * @tparam CI ContainerIndex
     */
    template<typename ES, typename DI, typename CI>
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

      template<typename size_type>
      friend struct ::Dune::PDELab::impl::update_ordering_data;

    public:

      static const bool has_dynamic_ordering_children = true;

      static const bool consume_tree_index = true;

      typedef LocalOrderingTraits<ES,DI,CI,MultiIndexOrder::Inner2Outer> Traits;

      static constexpr auto GT_UNUSED = ~std::size_t(0);

    protected:

      typedef impl::GridFunctionSpaceOrderingData<typename Traits::SizeType> GFSData;

    public:

      void map_local_index(const typename Traits::SizeType geometry_type_index,
                           const typename Traits::SizeType entity_index,
                           typename Traits::TreeIndexView mi,
                           typename Traits::ContainerIndex& ci) const
      {
        if (_child_count == 0) // leaf nodes
          {
            assert(mi.size() == 1 && "MultiIndex length must match GridFunctionSpace tree depth");
            ci.push_back(mi.back());
          }
        else // inner nodes
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


      /**
       * @brief Set last index of container indices
       * @details For a (view) list of dof indices, this method will fill or
       * append the last index of the same ammount of container indices.
       *
       * @tparam ItIn       DOFIndexViewIterator<DOFIterator>
       * @tparam ItOut      std::vector<CI>::iterator
       * @param begin       begin iterator for (view) list of dof indices
       * @param end         begin iterator for (view) list of dof indices
       * @param out         begin iterator of container indices to be filled out
       */
      template<typename ItIn, typename ItOut>
      void map_lfs_indices(const ItIn begin, const ItIn end, ItOut out) const
      {
        if (_child_count == 0) // leaf nodes
        {
          for (ItIn in = begin; in != end; ++in, ++out) {
            assert(in->size() == 1 &&
                   "MultiIndex length must match GridFunctionSpace tree depth");
            out->push_back(in->treeIndex().back());
          }
        } else if (_container_blocked) // blocked inner nodes
        {
          for (ItIn in = begin; in != end; ++in, ++out)
            out->push_back(in->treeIndex().back());
        } else if (_fixed_size) // non-blocked inner nodes with fixed sizes
        {
          for (ItIn in = begin; in != end; ++in, ++out) {
            const typename Traits::SizeType child_index =
                in->treeIndex().back();
            const typename Traits::SizeType gt_index =
                Traits::DOFIndexAccessor::geometryType(*in);
            if (child_index > 0) {
              const typename Traits::SizeType index =
                  gt_index * _child_count + child_index - 1;
              out->back() += _gt_dof_offsets[index];
            }
          }
        } else // non-blocked inner nodes with variable sizes
        {
          for (ItIn in = begin; in != end; ++in, ++out) {
            const typename Traits::SizeType child_index =
                in->treeIndex().back();
            if (child_index > 0) {
              const typename Traits::SizeType gt_index =
                  Traits::DOFIndexAccessor::geometryType(*in);
              const typename Traits::SizeType entity_index =
                  Traits::DOFIndexAccessor::entityIndex(*in);

              assert(_gt_used[gt_index]);

              const typename Traits::SizeType index =
                  (_gt_entity_offsets[gt_index] + entity_index) *
                      _child_count +
                  child_index - 1;
              out->back() += _entity_dof_offsets[index];
            }
          }
        }
      }

      template<typename CIOutIterator, typename DIOutIterator = DummyDOFIndexIterator>
      typename Traits::SizeType
      extract_entity_indices(const typename Traits::DOFIndex::EntityIndex& ei,
                             typename Traits::SizeType child_index,
                             CIOutIterator ci_out, const CIOutIterator ci_end,
                             DIOutIterator di_out = DIOutIterator()) const
      {
        typedef typename Traits::SizeType size_type;

        const size_type geometry_type_index = Traits::DOFIndexAccessor::GeometryIndex::geometryType(ei);
        const size_type entity_index = Traits::DOFIndexAccessor::GeometryIndex::entityIndex(ei);

        if (!_gt_used[geometry_type_index])
          return 0;

        if (_child_count == 0)
          {
            const size_type size = _fixed_size
              ? _gt_dof_offsets[geometry_type_index]
              : _entity_dof_offsets[(_gt_entity_offsets[geometry_type_index] + entity_index)];

            for (size_type i = 0; i < size; ++i, ++ci_out, ++di_out)
              {
                ci_out->push_back(i);
                di_out->treeIndex().push_back(i);
              }
            return size;
          }
        else
          {
            if (_container_blocked)
              {
                for (; ci_out != ci_end; ++ci_out)
                  {
                    ci_out->push_back(child_index);
                  }
              }
            else if (child_index > 0)
              {
                if (_fixed_size)
                  for (; ci_out != ci_end; ++ci_out)
                    {
                      const typename Traits::SizeType index = geometry_type_index * _child_count + child_index - 1;
                      ci_out->back() += _gt_dof_offsets[index];
                    }
                else
                  for (; ci_out != ci_end; ++ci_out)
                    {
                      const typename Traits::SizeType index = (_gt_entity_offsets[geometry_type_index] + entity_index) * _child_count + child_index - 1;
                      ci_out->back() += _entity_dof_offsets[index];
                    }
              }

            // The return value is not used for non-leaf orderings.
            return 0;
          }
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
        if (_fixed_size)
          return _child_count > 0
            ? _gt_dof_offsets[geometry_type_index * _child_count + _child_count - 1]
            : _gt_dof_offsets[geometry_type_index];

        if (!_gt_used[geometry_type_index])
          return 0;

        return _child_count > 0
          ? _entity_dof_offsets[(_gt_entity_offsets[geometry_type_index] + entity_index) * _child_count + _child_count - 1]
          : _entity_dof_offsets[(_gt_entity_offsets[geometry_type_index] + entity_index)];
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

    protected:

      /**
       * @brief Gives the size for a given entity and suffix
       * @details This method is used by typetree class derived from this class
       *
       * @param node TypeTree version of this class
       * @param suffix  MultiIndex with a partial path to a container
       * @param index Entity index to compute the size
       * @return Traits::SizeType  The size required for such a path.
       */
      template<class Node>
      typename Traits::SizeType
      node_size(const Node& node, typename Traits::ContainerIndex suffix,
           const typename Traits::DOFIndex::EntityIndex &index) const {
        using size_type = typename Traits::size_type;

        // suffix wants the size for this node
        if (suffix.size() == 0)
          return node.size(index);

        if constexpr (Node::isLeaf) {
          return 0; // Assume leaf local orderings are always field vectors
        } else {
          // the next index to find out its size
          auto back_index = suffix.back();
          // task: find child local ordering because it should know its own size
          std::size_t _child;

          if (node.containerBlocked()) {
            // in this case back index is the child ordering itself
            _child = back_index;
            suffix.pop_back();
          } else {
            // here we need to find the child that describes the back_index (solve child in map_lfs_indices)
            const size_type gt_index = Traits::DOFIndexAccessor::GeometryIndex::geometryType(index);
            const size_type entity_index = Traits::DOFIndexAccessor::GeometryIndex::entityIndex(index);
            auto dof_begin = node._fixed_size ? node._gt_dof_offsets.begin() : node._entity_dof_offsets.begin();
            auto dof_end = node._fixed_size ? node._gt_dof_offsets.end() : node._entity_dof_offsets.end();
            auto dof_it = std::prev(std::upper_bound(dof_begin, dof_end, back_index));
            size_type dof_dist = std::distance(dof_begin, dof_it);
            if (node._fixed_size)
              _child = dof_dist - gt_index * node._child_count + 1;
            else
              _child = dof_dist - (node._gt_entity_offsets[gt_index] + entity_index) * node._child_count + 1;
          }

          assert(node.degree() > _child);
          typename Traits::SizeType _size;
          // create a dynamic or static index range
          auto indices = Dune::range(node.degree());
          // get size for required child
          Hybrid::forEach(indices, [&](auto i){
            if (i == _child)
              _size = node.child(i).size(suffix, index);
          });

          return _size;
        }
      }

    public:
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
      LocalOrderingBase(Node& node, bool container_blocked, GFSData* gfs_data)
        : _fixed_size(false)
        , _fixed_size_possible(false)
        , _container_blocked(container_blocked)
        , _max_local_size(0)
        , _child_count(TypeTree::degree(node))
        , _children(TypeTree::degree(node),nullptr)
        , _gfs_data(gfs_data)
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

      bool contains_geometry_type(typename Traits::SizeType gt_index) const
      {
        return _gt_used[gt_index];
      }

      bool contains(typename Traits::SizeType codim) const
      {
        return _codim_used.test(codim);
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

      bool containerBlocked() const
      {
        return _container_blocked;
      }

      std::size_t childOrderingCount() const
      {
        return _child_count;
      }

      LocalOrderingBase& childOrdering(typename Traits::SizeType i)
      {
        return *_children[i];
      }

      const LocalOrderingBase& childOrdering(typename Traits::SizeType i) const
      {
        return *_children[i];
      }

      void disable_container_blocking()
      {
        _container_blocked = false;
      }

      //! Initial setup of the flag indicating whether a fixed size ordering is possible.
      /**
       * For a non-leaf ordering, a fixed size ordering is possible if all children can
       * support it, so we implement that logic here.
       *
       * \note Leaf orderings will usually want to extract this a priori information from somewhere
       * else, so they should override this method (the correct method will get called even
       * without a virtual call, as the call happens from a TypeTree visitor that is aware of
       * the precise type of the ordering).
       */
      void setup_fixed_size_possible()
      {
        _fixed_size_possible = true;
        for (const auto& child : _children)
          _fixed_size_possible &= child->_fixed_size_possible;
      }



      bool _fixed_size;
      bool _fixed_size_possible;
      bool _container_blocked;
      std::size_t _max_local_size;

      const std::size_t _child_count;
      std::vector<LocalOrderingBase*> _children;

      typename Traits::CodimFlag _codim_used;
      std::vector<bool> _gt_used;

      std::vector<typename Traits::SizeType> _gt_entity_offsets;
      std::vector<typename Traits::SizeType> _gt_dof_offsets;
      std::vector<typename Traits::SizeType> _entity_dof_offsets;

      GFSData* _gfs_data;

    };

    //! \}

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ORDERING_LOCALORDERINGBASE_HH
