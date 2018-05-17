//
// Created by marckoch on 14/05/18.
//

#ifndef DUNE_PDELAB_VISITORS_HH
#define DUNE_PDELAB_VISITORS_HH

#include <dune/typetree/visitor.hh>

namespace Dune{
  namespace Blockstructured{

    // the bogus template parameter is necessary to make GCC honor the friend declaration
    // in the LocalFunctionSpace (probably a GCC bug)
    template<typename = int>
    struct PropagateGlobalStorageVisitor
        : public TypeTree::TreeVisitor
            , public TypeTree::DynamicTraversal
    {

      template<typename LFS, typename Child, typename TreePath, typename ChildIndex>
      void beforeChild(const LFS& lfs, Child& child, TreePath treePath, ChildIndex childIndex) const
      {
        child._dof_index_storage_subentity_wise_ptr = lfs._dof_index_storage_subentity_wise_ptr;
      }
    };

    template<typename Entity, bool fast>
    struct ComputeSizeVisitor
        : public Dune::TypeTree::TreeVisitor
            , public TypeTree::DynamicTraversal
    {

      template<typename Node, typename TreePath>
      void pre(Node& node, TreePath treePath)
      {
        node.offset = offset;
        node.offsetLeafs = leafOffset;
      }

      template<typename Node, typename TreePath>
      void post(Node& node, TreePath treePath)
      {
          node.n = offset - node.offset;
          node.nLeafs = leafOffset - node.offsetLeafs;
      }

      template<typename Node, typename TreePath>
      void leaf(Node& node, TreePath treePath)
      {
          node.offset = offset;
          node.offsetLeafs = leafOffset;
          node.nLeafs = 1;
          if (fast)
            {
              node.pfe = nullptr;
              node.n = node.pgfs->finiteElementMap().maxLocalSize();
              Node::FESwitch::setStore(node.pfe, node.pgfs->finiteElementMap().find(e));
            }
          else
            {
              Node::FESwitch::setStore(node.pfe, node.pgfs->finiteElementMap().find(e));
              node.n = Node::FESwitch::basis(*node.pfe).size();
            }
          offset += node.n;
          leafOffset++;
      }

      ComputeSizeVisitor(const Entity& entity, std::size_t offset_ = 0)
          : e(entity)
          , offset(offset_)
          , leafOffset(0)
      {}

      const Entity& e;
      std::size_t offset;
      std::size_t leafOffset;

    };


    template<typename Entity, bool fast>
    struct FillIndicesVisitor
        : public TypeTree::TreeVisitor
            , public TypeTree::DynamicTraversal
    {

      template<typename Node, typename TreePath>
      void leaf(Node& node, TreePath treePath)
      {
        // setup DOFIndices for this finite element
        node.dofIndices(e, node.offsetLeafs);
      }

      template<typename Node, typename Child, typename TreePath, typename ChildIndex>
      void afterChild(const Node& node, const Child& child, TreePath treePath, ChildIndex childIndex)
      {
        // Just skip the entire function space structure handling in fast mode
        // This **really** breaks any attempt at using the DOFIndex for anything other
        // than as input to the FastDGGridOperator machine.
        // You have been warned!
        for (int i = 0; i < child.nLeafs; ++i)
          for (auto &codim: (*node._dof_index_storage_subentity_wise_ptr)[child.offsetLeafs + i])
            for (auto &subentity: codim)
              subentity.treeIndex().push_back(childIndex);
      }

      FillIndicesVisitor(const Entity& entity)
          : e(entity)
      {}

      const Entity& e;
    };


  }
}

#endif //DUNE_PDELAB_VISITORS_HH
