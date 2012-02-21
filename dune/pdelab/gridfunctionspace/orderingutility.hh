// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_ORDERINGUTILITY_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_ORDERINGUTILITY_HH

#include <vector>

#include <dune/pdelab/common/multiindex.hh>
#include <dune/pdelab/common/typetree/traversal.hh>
#include <dune/pdelab/common/typetree/accumulate_static.hh>

namespace Dune {
  namespace PDELab {

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


    template<typename DI, typename CI>
    struct OrderingTraits
    {
      typedef DI DOFIndex;

      typedef typename DI::TreeIndex TreeIndex;

      typedef typename DI::View DOFIndexView;
      typedef typename DI::View::TreeIndex TreeIndexView;

      typedef CI ContainerIndex;

      typedef typename DI::size_type SizeType;
    };


    template<typename GV, typename DI, typename CI>
    struct LocalOrderingTraits
      : public OrderingTraits<DI,CI>
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


    template<typename DI, typename CI>
    class VirtualOrderingBase
    {
    public:

      typedef OrderingTraits<DI,CI> Traits;

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

   //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_ORDERINGUTILITY_HH
