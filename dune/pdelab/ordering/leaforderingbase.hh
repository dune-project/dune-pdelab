// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_ORDERING_LEAFORDERINGBASE_HH
#define DUNE_PDELAB_ORDERING_LEAFORDERINGBASE_HH

#include <dune/typetree/typetree.hh>

#include <dune/pdelab/ordering/utility.hh>
#include <dune/pdelab/ordering/orderingbase.hh>
#include <dune/pdelab/ordering/directleaflocalordering.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup Ordering
    //! \{

    //! Generic infrastructure for orderings for leaf spaces
    template<typename LocalOrdering>
    class LeafOrderingBase
      : public TypeTree::CompositeNode<LocalOrdering>
      , public VirtualOrderingBase<typename LocalOrdering::Traits::DOFIndex,
                                   typename LocalOrdering::Traits::ContainerIndex>
      , public OrderingBase<typename LocalOrdering::Traits::DOFIndex,
                            typename LocalOrdering::Traits::ContainerIndex>
    {
    public:
      typedef typename LocalOrdering::Traits Traits;

      static const bool has_dynamic_ordering_children = false;

      static const bool consume_tree_index = false;

    protected:

      typedef TypeTree::CompositeNode<LocalOrdering> NodeT;

      typedef OrderingBase<typename LocalOrdering::Traits::DOFIndex,
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


      LeafOrderingBase(const typename NodeT::NodeStorage& local_ordering, bool container_blocked, typename BaseT::GFSData* gfs_data)
        : NodeT(local_ordering)
        , BaseT(*this,container_blocked,gfs_data,this)
      {}

#ifndef DOXYGEN

// we need to override the default copy / move ctor to fix the delegate pointer, but that is
// hardly interesting to our users...

      LeafOrderingBase(const LeafOrderingBase& r)
        : NodeT(r.nodeStorage())
        , BaseT(r)
        , _gt_dof_offsets(r._gt_dof_offsets)
      {
        this->setDelegate(this);
      }

      LeafOrderingBase(LeafOrderingBase&& r)
        : NodeT(r.nodeStorage())
        , BaseT(std::move(r))
        , _gt_dof_offsets(std::move(r._gt_dof_offsets))
      {
        this->setDelegate(this);
      }

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
                // This check is needed to avoid a horrid stream of compiler warnings about
                // exceeding array bounds in ReservedVector!
                if (ci.size() < ci.capacity())
                  ci.push_back(_gt_dof_offsets[geometry_type_index] + entity_index);
                else
                  {
                    DUNE_THROW(Dune::Exception,"Container blocking incompatible with backend structure");
                  }
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
                    assert(localOrdering()._gt_used[geometry_type_index]);
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
                    assert(localOrdering()._gt_used[geometry_type_index]);
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
                    assert(localOrdering()._gt_used[geometry_type_index]);
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
                    assert(localOrdering()._gt_used[geometry_type_index]);
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

      /**
       The actual implementation has to implement an update() method to fill ...
      */
      virtual void update() = 0;

      using BaseT::fixedSize;

    protected:

      using BaseT::_max_local_size;
      using BaseT::_size;
      using BaseT::_block_count;
      using BaseT::_container_blocked;
      using BaseT::_fixed_size;
      using BaseT::_codim_used;
      using BaseT::_codim_fixed_size;

      std::vector<typename Traits::SizeType> _gt_dof_offsets;

    };

   //! \} group Ordering
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ORDERING_LEAFORDERINGBASE_HH
