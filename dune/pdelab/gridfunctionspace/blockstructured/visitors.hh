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
        : public Dune::PDELab::PropagateGlobalStorageVisitor<int>
    {
      template<typename LFS, typename Child, typename TreePath, typename ChildIndex>
      void beforeChild(const LFS& lfs, Child& child, TreePath treePath, ChildIndex childIndex) const
      {
        Dune::PDELab::PropagateGlobalStorageVisitor<int>::beforeChild(lfs, child, treePath, childIndex);
        child._subentityWiseDOFs_ptr = lfs._subentityWiseDOFs_ptr;
      }
    };

    template<typename Entity, bool fast>
    struct ComputeSizeVisitor
        : Dune::PDELab::ComputeSizeVisitor<Entity, false>
    {
      using Base = Dune::PDELab::ComputeSizeVisitor<Entity, false>;

      template<typename Node, typename TreePath>
      void pre(Node& node, TreePath treePath)
      {
        Base::pre(node, treePath);
        node.offsetLeafs = leafOffset;
      }

      template<typename Node, typename TreePath>
      void post(Node& node, TreePath treePath)
      {
        Base::post(node, treePath);
        node.nLeafs = leafOffset - node.offsetLeafs;
      }

      template<typename Node, typename TreePath>
      void leaf(Node& node, TreePath treePath)
      {
        Base::leaf(node, treePath);
        node.offsetLeafs = leafOffset;
        node.nLeafs = 1;
        leafOffset++;
      }

      ComputeSizeVisitor(const Entity& entity, std::size_t offset_ = 0)
          : Base(entity, offset_),
            leafOffset(0)
      {}

      std::size_t leafOffset;
    };


    template<typename Entity, bool fast>
    struct FillIndicesVisitor
        : public Dune::PDELab::FillIndicesVisitor<Entity, false>
    {
      template<typename Node, typename Child, typename TreePath, typename ChildIndex>
      void afterChild(const Node& node, const Child& child, TreePath treePath, ChildIndex childIndex)
      {
        for (int i = 0; i < child.nLeafs; ++i)
          for (auto &codim: (*node._subentityWiseDOFs_ptr)[child.offsetLeafs + i])
            for (auto &subentity: codim)
              subentity.treeIndex().push_back(childIndex);
      }

      FillIndicesVisitor(const Entity& entity)
          : Dune::PDELab::FillIndicesVisitor<Entity,false>(entity)
      {}
    };

  }
}

#endif //DUNE_PDELAB_VISITORS_HH
