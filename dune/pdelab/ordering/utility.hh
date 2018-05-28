// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_ORDERING_UTILITY_HH
#define DUNE_PDELAB_ORDERING_UTILITY_HH

#include <vector>
#include <bitset>

#include <dune/pdelab/common/dofindex.hh>
#include <dune/pdelab/common/globaldofindex.hh>
#include <dune/pdelab/ordering/transformations.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup Ordering
    //! \{

    //! Index merging algorithm for global orderings.
    struct MergeMode
    {

      enum type {
        lexicographic, //!< Lexicographically ordered ([i1,i2],[j1,j2] -> [i1,i2,j1,j2]).
        interleaved  //!< Indices are interleaved according to a user-supplied pattern ([i1,i2],[j1,j2] -> [i1,j1,i2,j2]).
      };

    };

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

    } // end namespace ordering

#endif // DOXYGEN


    struct DefaultDOFIndexAccessor
    {

      template<typename DOFIndex, typename SizeType, typename IndexType>
      static typename std::enable_if<
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
      static typename std::enable_if<
        !std::is_integral<IndexType>::value
        >::type
      store(DOFIndex& dof_index, const GeometryType& gt, SizeType entity_index, IndexType tree_index)
      {
        dof_index.entityIndex()[0] = GlobalGeometryTypeIndex::index(gt);
        dof_index.entityIndex()[1] = entity_index;
        dof_index.treeIndex() = tree_index;
      }

      template<typename DOFIndex, typename SizeType, typename IndexType>
      static typename std::enable_if<
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
      static typename std::enable_if<
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



    template<typename DI, typename CI>
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
      typedef typename DI::size_type size_type;

    };


    template<typename ES, typename DI, typename CI>
    struct LocalOrderingTraits
      : public OrderingTraits<DI,
                              CI
                              >
    {

      using EntitySet = ES;
      using GridView = typename ES::GridView;

    };

    template<typename ES, typename DI, typename CI>
    struct GridViewOrderingTraits
      : public LocalOrderingTraits<ES,DI,CI>
    {

      typedef typename DI::EntityIndex EntityIndex;
      typedef typename DI::View::EntityIndex EntityIndexView;

    };


    template<typename DI, typename CI>
    class VirtualOrderingBase
    {
    public:

      typedef OrderingTraits<DI,CI> Traits;

      VirtualOrderingBase() {}

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
      typename std::enable_if<Node::has_dynamic_ordering_children>::type
      extract_child(const Node& node, Child& child, ChildIndex child_index)
      {
        _children[child_index] = &child;
      }

      template<typename Node, typename Child, typename ChildIndex>
      typename std::enable_if<!Node::has_dynamic_ordering_children>::type
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

   //! \} group Ordering
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ORDERING_UTILITY_HH
