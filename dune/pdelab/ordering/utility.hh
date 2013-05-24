// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_ORDERING_UTILITY_HH
#define DUNE_PDELAB_ORDERING_UTILITY_HH

#include <vector>
#include <bitset>

#include <dune/pdelab/common/dofindex.hh>
#include <dune/pdelab/common/globaldofindex.hh>
#include <dune/pdelab/common/typetree/traversal.hh>
#include <dune/pdelab/common/typetree/accumulate_static.hh>

#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>

namespace Dune {
  namespace PDELab {

#ifndef DOXYGEN

    namespace ordering {

      // This is an implementation detail of the composite orderings, no need to confuse our users!
      struct update_direct_children
        : public TypeTree::DirectChildrenVisitor
        , public TypeTree::DynamicTraversal
      {

        template<typename GFS, typename Child, typename TreePath, typename ChildIndex>
        void afterChild(const GFS& gfs, Child& child, TreePath tp, ChildIndex childIndex) const
        {
          child.update();
        }

      };

    }

#endif // DOXYGEN

    struct extract_max_container_depth
    {

      typedef std::size_t result_type;

      template<typename Node, typename TreePath>
      struct doVisit
      {
        static const bool value = true;
      };

      template<typename Node, typename TreePath>
      struct visit
      {
        static const std::size_t result = Node::Traits::Backend::Traits::max_blocking_depth;
      };

    };

    template<typename RootGFS>
    struct gfs_to_ordering
    {
      static const std::size_t ci_depth =
        TypeTree::AccumulateValue<RootGFS,
                                  extract_max_container_depth,
                                  TypeTree::max<std::size_t>,
                                  0,
                                  TypeTree::plus<std::size_t>
                                  >::result + 1;

      typedef typename gfs_to_lfs<RootGFS>::DOFIndex DOFIndex;
      typedef MultiIndex<std::size_t,ci_depth> ContainerIndex;
    };

    template<typename GlobalTransformation>
    struct gfs_to_local_ordering
    {
      typedef typename GlobalTransformation::DOFIndex DOFIndex;
      typedef typename GlobalTransformation::ContainerIndex ContainerIndex;
    };


    struct DefaultDOFIndexAccessor
    {

      template<typename DOFIndex, typename SizeType, typename IndexType>
      static typename enable_if<
        std::is_integral<IndexType>::value
        >::type
      store(DOFIndex& dof_index, const GeometryType& gt, SizeType entity_index, IndexType tree_index)
      {
        dof_index.clear();
        dof_index.entityIndex()[0] = GlobalGeometryTypeIndex::index(gt);
        dof_index.entityIndex()[1] = entity_index;
        dof_index.treeIndex().push_back(tree_index);
      }

      template<typename DOFIndex, typename SizeType, typename IndexType>
      static typename enable_if<
        !std::is_integral<IndexType>::value
        >::type
      store(DOFIndex& dof_index, const GeometryType& gt, SizeType entity_index, IndexType tree_index)
      {
        dof_index.entityIndex()[0] = GlobalGeometryTypeIndex::index(gt);
        dof_index.entityIndex()[1] = entity_index;
        dof_index.treeIndex() = tree_index;
      }

      template<typename DOFIndex, typename SizeType, typename IndexType>
      static typename enable_if<
        std::is_integral<IndexType>::value
        >::type
      store(DOFIndex& dof_index, SizeType gt_index, SizeType entity_index, IndexType tree_index)
      {
        dof_index.clear();
        dof_index.entityIndex()[0] = gt_index;
        dof_index.entityIndex()[1] = entity_index;
        dof_index.treeIndex().push_back(tree_index);
      }

      template<typename DOFIndex, typename SizeType, typename IndexType>
      static typename enable_if<
        !std::is_integral<IndexType>::value
        >::type
      store(DOFIndex& dof_index, SizeType gt_index, SizeType entity_index, IndexType tree_index)
      {
        dof_index.entityIndex()[0] = gt_index;
        dof_index.entityIndex()[1] = entity_index;
        dof_index.treeIndex() = tree_index;
      }


      struct GeometryIndex
      {

        template<typename Index>
        static std::size_t geometryType(const Index& geometry_index)
        {
          return geometry_index[0];
        }

        template<typename Index>
        static std::size_t entityIndex(const Index& geometry_index)
        {
          return geometry_index[1];
        }

        template<typename Index, typename SizeType>
        static void store(Index& index, const GeometryType& gt, SizeType entity_index)
        {
          index[0] = GlobalGeometryTypeIndex::index(gt);
          index[1] = entity_index;
        }

      };

      template<typename DOFIndex>
      static std::size_t geometryType(const DOFIndex& dof_index)
      {
        return GeometryIndex::geometryType(dof_index.entityIndex());
      }

      template<typename DOFIndex>
      static std::size_t entityIndex(const DOFIndex& dof_index)
      {
        return GeometryIndex::entityIndex(dof_index.entityIndex());
      }

    };

    struct SimpleDOFIndexAccessor
    {

      template<typename DOFIndex, typename SizeType>
      static void store(DOFIndex& dof_index, const GeometryType& gt, SizeType entity_index, SizeType tree_index)
      {
        dof_index = entity_index;
      }

    };


    template<typename DI, typename CI>
    struct SimpleOrderingTraits
    {

      typedef DI DOFIndex;

      typedef CI ContainerIndex;

      typedef std::size_t SizeType;

      typedef DefaultDOFIndexAccessor DOFIndexAccessor;

    };


    template<typename SizeType_, typename CI>
    struct SimpleOrderingTraits<SimpleDOFIndex<SizeType_>,CI>
    {

      typedef SimpleDOFIndex<SizeType_> DOFIndex;

      typedef CI ContainerIndex;

      typedef SizeType_ SizeType;

      typedef SimpleDOFIndexAccessor DOFIndexAccessor;

    };



    template<typename DI, typename GDI, typename CI>
    struct OrderingTraits
      : public SimpleOrderingTraits<DI,CI>
    {

      // The maximum dimension supported (length of bitsets)
      // 32 dimensions should probably be fine for now... ;-)
      static const std::size_t max_dim = 32;

      typedef std::bitset<max_dim> CodimFlag;

      typedef typename DI::TreeIndex TreeIndex;

      typedef typename DI::View DOFIndexView;
      typedef typename DI::View::TreeIndex TreeIndexView;

      typedef typename DI::size_type SizeType;

      typedef GDI GlobalDOFIndex;

    };


    template<typename GV, typename DI, typename CI>
    struct LocalOrderingTraits
      : public OrderingTraits<DI,
                              Dune::PDELab::GlobalDOFIndex<
                                typename DI::value_type,
                                DI::max_depth,
                                typename GV::Grid::GlobalIdSet::IdType
                                >,
                              CI
                              >
    {

      typedef GV GridView;

    };

    template<typename GV, typename DI, typename CI>
    struct GridViewOrderingTraits
      : public LocalOrderingTraits<GV,DI,CI>
    {

      typedef typename DI::EntityIndex EntityIndex;
      typedef typename DI::View::EntityIndex EntityIndexView;

    };


    template<typename DI, typename GDI, typename CI>
    class VirtualOrderingBase
    {
    public:

      typedef OrderingTraits<DI,GDI,CI> Traits;

      virtual void map_index_dynamic(typename Traits::DOFIndexView di, typename Traits::ContainerIndex& ci) const = 0;
    };


    template<typename child_type>
    struct extract_child_bases
      : public TypeTree::DirectChildrenVisitor
      , public TypeTree::DynamicTraversal
    {

      template<typename Node, typename Child, typename TreePath, typename ChildIndex>
      void afterChild(const Node& node, Child& child, TreePath tp, ChildIndex child_index)
      {
        extract_child(node,child,child_index);
      }

      template<typename Node, typename Child, typename ChildIndex>
      typename enable_if<Node::has_dynamic_ordering_children>::type
      extract_child(const Node& node, Child& child, ChildIndex child_index)
      {
        _children[child_index] = &child;
      }

      template<typename Node, typename Child, typename ChildIndex>
      typename enable_if<!Node::has_dynamic_ordering_children>::type
      extract_child(const Node& node, Child& child, ChildIndex child_index)
      {
      }

      extract_child_bases(std::vector<child_type*>& children)
        : _children(children)
      {}

    private:
      std::vector<child_type*>& _children;

    };


    //! Dummy iterator type over DOF indices.
    /**
     * This dummy iterator is used to support omitting the calculation
     * of DOFIndex values in the per-entity index lookup methods of
     * orderings. By defining all operations performed on the DOFIndex
     * iterator and its value by this methods as no-ops, we can reuse the
     * combined implementation mapping both DOFIndex and ContainerIndex for
     * the (much more common) case of only having to map the ContainerIndex
     * values.
     */
    struct DummyDOFIndexIterator
    {

      typedef std::size_t size_type;

      DummyDOFIndexIterator& operator++()
      {
        return *this;
      }

      DummyDOFIndexIterator& operator+=(size_type i)
      {
        return *this;
      }

      DummyDOFIndexIterator& operator*()
      {
        return *this;
      }

      DummyDOFIndexIterator* operator->()
      {
        return this;
      }

      DummyDOFIndexIterator& treeIndex()
      {
        return *this;
      }

      bool operator==(const DummyDOFIndexIterator& r) const
      {
        return true;
      }

      bool operator!=(const DummyDOFIndexIterator& r) const
      {
        return !operator==(r);
      }

      void push_back(size_type i)
      {}

    };


    //! Mixin class for providing information about contained grid partitions.
    /**
     * This is a mixin class for orderings providing the common implementation
     * of the Dune::PartitionType query interface. As the number of partition types
     * is fixed, we can easily move the complete implementation into this mixin,
     * only requiring the ordering to update the contained information using the
     * protected API.
     */
    class PartitionInfoProvider
    {

    public:

      //! Returns whether this ordering contains entities with PartitionType partition.
      bool containsPartition(PartitionType partition) const
      {
        return _contained_partitions.test(static_cast<unsigned char>(partition));
      }

      //! Returns the internal representation of the set of contained entities.
      std::bitset<6> containedPartitions() const
      {
        return _contained_partitions;
      }

    protected:

      //! Empties the set of contained partitions.
      void clearPartitionSet()
      {
        _contained_partitions.reset();
      }

      //! Adds all partitions contained in r the set of contained partitions.
      void mergePartitionSet(const PartitionInfoProvider& r)
      {
        _contained_partitions |= r._contained_partitions;
      }

      //! Sets the set of contained partitions to the passed-in value.
      /**
       * \warning This is an internal interface that relies on the internal implementation
       *          of the partition set and may change without notice!
       */
      void setPartitionSet(const std::bitset<6>& partitions)
      {
        _contained_partitions = partitions;
      }

      //! Copies the set of contained partitions from r.
      void setPartitionSet(const PartitionInfoProvider& r)
      {
        _contained_partitions = r._contained_partitions;
      }

      //! Adds the partitions from all PartitionInfoProviders in the range [begin,end).
      /**
       * \note The passed-in iterators may yield both references and pointers to the
       *       PartitionInfoProviders in the range. This feature exists mostly to simplify
       *       implementation of the dynamic ordering base classes, which hold pointers
       *       to their children.
       */
      template<typename It>
      void mergePartitionSets(It begin, It end)
      {
        clearPartitionSet();
        for (; begin != end; ++begin)
          mergePartitionSet(reference(*begin));
      }

    private:

      // ********************************************************************************
      // The following two function are here to make mergePartitionSets() work with both
      // normal iterators and iterators of pointers. It would be nice to just use
      // boost::indirect_iterator, but alas...
      // ********************************************************************************

      static const PartitionInfoProvider& reference(const PartitionInfoProvider& provider)
      {
        return provider;
      }

      static const PartitionInfoProvider& reference(const PartitionInfoProvider* provider)
      {
        return *provider;
      }

      std::bitset<6> _contained_partitions;

    };

   //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ORDERING_UTILITY_HH
