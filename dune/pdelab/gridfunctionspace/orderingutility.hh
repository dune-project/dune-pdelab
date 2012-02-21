// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_ORDERINGUTILITY_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_ORDERINGUTILITY_HH

#include <vector>

namespace Dune {
  namespace PDELab {

    template<typename MI, typename CI>
    struct OrderingTraits
    {
      typedef MI MultiIndex;
      typedef CI ContainerIndex;

      typedef typename MI::Traits::SizeType SizeType;
    };

    template<typename MI, typename CI>
    class VirtualOrderingBase
    {
    public:

      typedef OrderingTraits<MI,CI> Traits;

      virtual void map_index_dynamic(const typename Traits::MultiIndex& mi, typename Traits::ContainerIndex& ci) const = 0;
    };


    template<typename child_type>
    struct extract_child_bases
      : public TypeTree::DirectChildrenVisitor
      , public TypeTree::DynamicTraversal
    {

      template<typename Node, typename Child, typename TreePath, typename ChildIndex>
      void afterChild(const Node& node, Child& child, TreePath tp, ChildIndex child_index);
      {
        _children[child_index] = &child;
      }

      extract_child_bases(std::vector<child_type*>& children)
        : _children(children)
      {}

    private:
      std::vector<child_type*> _children;

    };

   //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_ORDERINGUTILITY_HH
