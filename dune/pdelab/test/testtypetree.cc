#include "config.h"
#include <dune/pdelab/common/leafnode.hh>
#include <dune/pdelab/common/powernode.hh>
//#include <dune/pdelab/common/compositenode.hh>
#include <dune/pdelab/common/traversal.hh>

#include <iostream>

struct SimpleLeaf
  : public Dune::PDELab::TypeTree::LeafNode
{

  static const char* name()
  {
    return "SimpleLeaf";
  }
};

template<typename T, std::size_t k>
struct SimplePower
  : public Dune::PDELab::TypeTree::PowerNode<T,k>
{
  static const char* name()
  {
    return "SimplePower";
  }
};

struct TreePrinter
  : public Dune::PDELab::TypeTree::TypeTreeVisitor
{

  template<typename T, typename TreePath>
  void leaf(const T& t, TreePath treePath) const
  {
    pre(t,treePath);
  }

  template<typename T, typename TreePath>
  void pre(const T& t, TreePath treePath) const
  {
    for (std::size_t i = 0; i < Dune::PDELab::TypeTree::TreePathSize<TreePath>::value; ++i)
      std::cout << "  ";
    std::cout << t.name() << std::endl;
  }
};

int main(int argc, char** argv)
{
  SimpleLeaf sl1;
  SimplePower<SimpleLeaf,3> sp1;
  sp1.setChild(0,sl1);
  sp1.setChild(1,sl1);
  sp1.setChild(2,sl1);
  sp1.getChild<0>();
  //sp1.getChild(1);
  Dune::PDELab::TypeTree::applyToTree(sp1,TreePrinter());
  return 0;
}
