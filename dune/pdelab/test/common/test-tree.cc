#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/pdelab/concepts/treenode.hh>
#include <dune/pdelab/concepts/tree.hh>

#include <dune/pdelab/common/tree_traversal.hh>

#include <dune/typetree/leafnode.hh>
#include <dune/typetree/powernode.hh>
#include <dune/typetree/dynamicpowernode.hh>
#include <dune/typetree/compositenode.hh>

struct LeafNode : public Dune::TypeTree::LeafNode {
  LeafNode() = default;
};

struct ArrayNode : public Dune::TypeTree::PowerNode<LeafNode,2> {
  ArrayNode() : Dune::TypeTree::PowerNode<LeafNode,2>{LeafNode{}, LeafNode{}} {}
};

struct VectorNode : public Dune::TypeTree::DynamicPowerNode<ArrayNode> {
  VectorNode() : Dune::TypeTree::DynamicPowerNode<ArrayNode>{ArrayNode{}, ArrayNode{}} {}
};

struct TupleNode : public Dune::TypeTree::CompositeNode<ArrayNode, VectorNode, LeafNode> {
  TupleNode() : Dune::TypeTree::CompositeNode<ArrayNode, VectorNode, LeafNode>{ArrayNode{}, VectorNode{}, LeafNode{}} {}
};

int main()
{
  static_assert(Dune::PDELab::Concept::LeafTreeNode<LeafNode>);
  static_assert(Dune::PDELab::Concept::LeafTreeNode<LeafNode&>);
  static_assert(Dune::PDELab::Concept::LeafTreeNode<const LeafNode&>);
  static_assert(Dune::PDELab::Concept::TreeNode<LeafNode>);

  static_assert(not Dune::PDELab::Concept::ArrayTreeNode<LeafNode>);
  static_assert(not Dune::PDELab::Concept::VectorTreeNode<LeafNode>);

  static_assert(Dune::PDELab::Concept::ArrayTreeNode<ArrayNode>);
  static_assert(Dune::PDELab::Concept::TreeNode<ArrayNode>);
  static_assert(Dune::PDELab::Concept::VectorTreeNode<VectorNode>);
  static_assert(Dune::PDELab::Concept::TreeNode<VectorNode>);
  static_assert(Dune::PDELab::Concept::TupleTreeNode<TupleNode>);
  static_assert(Dune::PDELab::Concept::TreeNode<TupleNode>);

  static_assert(Dune::PDELab::Concept::TupleTreeNode<TupleNode&>);
  static_assert(Dune::PDELab::Concept::TupleTreeNode<const TupleNode&>);
  static_assert(Dune::PDELab::Concept::TupleTreeNode<TupleNode&&>);

  Dune::PDELab::forEach(ArrayNode{}, [](auto&& node){});
  Dune::PDELab::forEach(TupleNode{}, [](auto&& node){});
  Dune::PDELab::forEachNode(TupleNode{}, [](auto&& node){});
  Dune::PDELab::forEachLeafNode(TupleNode{}, [](auto&& node){});

  static_assert(Dune::PDELab::Concept::Tree<ArrayNode>);
  static_assert(Dune::PDELab::Concept::Tree<VectorNode>);
  static_assert(Dune::PDELab::Concept::Tree<TupleNode>);
}
