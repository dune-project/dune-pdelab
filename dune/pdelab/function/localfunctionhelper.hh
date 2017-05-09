//-*- tab-width: 4; c-basic-offset: 2; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_FUNCTION_LOCALFUNCTIONHELPER_HH
#define DUNE_PDELAB_FUNCTION_LOCALFUNCTIONHELPER_HH

#include <dune/pdelab/function/tags.hh>
#include <dune/typetree/visitor.hh>

namespace Dune {
namespace PDELab {

namespace Imp
{

  template<typename Entity>
  struct PowerCompositeBindVisitor
    : public TypeTree::TreeVisitor, public TypeTree::DynamicTraversal
  {
    PowerCompositeBindVisitor(const Entity & e) : e_(e) {}
    template<typename LeafNode, typename TreePath>
    void leaf(LeafNode& node, TreePath treePath) const
    {
      node.bind(e_);
    }
    const Entity & e_;
  };

  struct PowerCompositeUnbindVisitor
    : public TypeTree::TreeVisitor, public TypeTree::DynamicTraversal
  {
    template<typename LeafNode, typename TreePath>
    void leaf(LeafNode& node, TreePath treePath) const
    {
      node.unbind();
    }
  };

  template<typename F>
  class LocalFunctionLeafNodeWrapper
    : public TypeTree::LeafNode
    , public F
  {
  public:
    typedef DifferentiableFunctionLocalViewTag ImplementationTag;
    LocalFunctionLeafNodeWrapper(const F& f) :
      F(f)
    {}
  };

} // end namespace Imp
} // end namespace PDELab
} // end namespace Dune

#endif // DUNE_PDELAB_FUNCTION_LOCALFUNCTIONHELPER_HH
