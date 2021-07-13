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

        template<typename Index, typename SizeType>
        static void store(Index& index, SizeType geometry_index, SizeType entity_index)
        {
          index[0] = geometry_index;
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

      using SizePrefix = MultiIndex<SizeType,ContainerIndex::max_depth>;

      typedef DefaultDOFIndexAccessor DOFIndexAccessor;

    };


    template<typename SizeType_, typename CI>
    struct SimpleOrderingTraits<SimpleDOFIndex<SizeType_>,CI>
    {

      typedef SimpleDOFIndex<SizeType_> DOFIndex;

      typedef CI ContainerIndex;

      typedef SizeType_ SizeType;

      using SizePrefix = MultiIndex<SizeType,ContainerIndex::max_depth>;

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
      virtual ~VirtualOrderingBase() = default;

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

    //! Information about order semantics on multi-indices
    enum class MultiIndexOrder {
      //! indices are ordered from inner to outer container: {inner,...,outer}
      Inner2Outer,
      //! indices are ordered from outer to inner container: {outer,...,inner}
      Outer2Inner,
    };

    /**
     * @brief Adapter to create a size provider from an ordering
     * @details This adapter is meant to be used in allocation and
     *   resizing of vectors containers.
     *   In particular, this adapter is needed because the ordering library give
     *   sizes for multi-indices ordered with Inner2Outer semantis, while
     *   resizing algorithms are faster and easier when using Outer2Inner
     *   semantics.
     *
     *   - This class makes type erasure on the ordering type.
     *   - This class has value semantics.
     *
     * @tparam SizeType_ return type of the size method
     * @tparam SizePrefix_ argument type of the size method
     * @tparam Order enum with MultiIndexOrder semantics of the adapter
     *
     * @tparam Ordering  Ordering tree from grid function space
     */
    template <class SizeType_, class SizePrefix_, MultiIndexOrder Order = MultiIndexOrder::Outer2Inner>
    struct SizeProviderAdapter {

      template <class Ordering>
      SizeProviderAdapter(const std::shared_ptr<const Ordering>& ordering)
          : _size_provider(
                [=](const SizePrefix &size) { return ordering->size(size); })
      {
        static_assert(Ordering::size_prefix_order == MultiIndexOrder::Inner2Outer);
      }

      //! Partial MultiIndex of a ContainerIndex
      using SizePrefix = SizePrefix_;

      //! Type that refers to the size of containers
      using SizeType = SizeType_;

      //! Inform about SizePrefix multi-index order semantics
      static constexpr MultiIndexOrder size_prefix_order = Order;

      /**
       * @brief Gives the size for a given prefix
       * @param prefix  MultiIndex with a partial path to a container
       * @return Traits::SizeType  The size required for such a path
       */
      auto size(const SizePrefix &prefix) const {
        if constexpr (size_prefix_order == MultiIndexOrder::Outer2Inner) {
          // reversing Outer2Inner prefix into a Inner2Outer prefix
          SizePrefix i2o_prefix;
          i2o_prefix.resize(prefix.size());
          std::reverse_copy(prefix.begin(), prefix.end(),i2o_prefix.begin());
          // forward size request to ordering with new Inner2Outer prefix
          return _size_provider(i2o_prefix);
        } else {
          // prefix is already Inner2Outer, forward to size provider
          return _size_provider(prefix);
        }
      }

    private:

      const std::function<SizeType(const SizePrefix &)> _size_provider;
    };

    //! template deduction guide for orderings
    template<class Ordering>
    SizeProviderAdapter(const std::shared_ptr<const Ordering>& ordering)
        -> SizeProviderAdapter<typename Ordering::Traits::SizeType,
                              typename Ordering::Traits::SizePrefix>;


   //! \} group Ordering
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ORDERING_UTILITY_HH
