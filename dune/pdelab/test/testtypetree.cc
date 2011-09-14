#include "config.h"

#include "typetreetestswitch.hh"

#if TEST_TYPETREE_INVALID

int main()
{
  return 0;
}

#else

#include "typetreetestutility.hh"


struct NodeCountingFunctor
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
    static const result_type result = 1;
  };

};

struct LeafCountingFunctor
{

  typedef std::size_t result_type;

  template<typename Node, typename TreePath>
  struct doVisit
  {
    static const bool value = Node::isLeaf;
  };

  template<typename Node, typename TreePath>
  struct visit
  {
    static const result_type result = 1;
  };

};


struct DepthFunctor
{

  typedef std::size_t result_type;

  template<typename Node, typename TreePath>
  struct doVisit
  {
    static const bool value = Node::isLeaf;
  };

  template<typename Node, typename TreePath>
  struct visit
  {
    // the TreePath is always one entry shorter than the actual depth of the tree
    static const result_type result = Dune::PDELab::TypeTree::TreePathSize<TreePath>::value + 1;
  };

};


int main(int argc, char** argv)
{

  // basic tests

  // leaf node
  TreePrinter treePrinter;
  SimpleLeaf sl1;

  Dune::PDELab::TypeTree::applyToTree(sl1,treePrinter);

  typedef SimplePower<SimpleLeaf,3> SP1;
  SP1 sp1_1;
  sp1_1.setChild(0,sl1);
  sp1_1.setChild(1,sl1);
  sp1_1.setChild(2,sl1);

  Dune::PDELab::TypeTree::applyToTree(sp1_1,TreePrinter());

  SimpleLeaf sl2;
  SP1 sp1_2(sl2,false);

  Dune::PDELab::TypeTree::applyToTree(sp1_2,TreePrinter());

  SP1 sp1_2a(sl2,true);

  Dune::PDELab::TypeTree::applyToTree(sp1_2a,TreePrinter());


  typedef SimpleComposite<SimpleLeaf,SP1,SimpleLeaf> SC1;
  SC1 sc1_1(sl1,sp1_2,sl2);
  Dune::PDELab::TypeTree::applyToTree(const_cast<const SC1&>(sc1_1),treePrinter);

#if HAVE_VARIADIC_TEMPLATES

#if HAVE_RVALUE_REFERENCES

  typedef SimpleComposite<SimpleLeaf,SimpleLeaf,SimpleLeaf> SC2;
  SC2 sc2(sl1,sl1,sl1);

  typedef SimpleVariadicComposite<SimpleLeaf,SP1,SimpleLeaf,SC1> SVC1;
  SVC1 svc1_1(sl1,sp1_1,sl2,sc1_1);
  Dune::PDELab::TypeTree::applyToTree(svc1_1,treePrinter);

  SP1 sp1_3(SimpleLeaf(),SimpleLeaf(),sl1);
  Dune::PDELab::TypeTree::applyToTree(sp1_3,TreePrinter());

#if HAVE_VARIADIC_CONSTRUCTOR_SFINAE

  SVC1 svc1_2(SimpleLeaf(),SP1(sp1_2),sl2,const_cast<const SC1&>(sc1_1));
  Dune::PDELab::TypeTree::applyToTree(svc1_2,TreePrinter());

  typedef SimpleComposite<SimpleLeaf,SC2,SimpleLeaf,SC1> SVC2;
  SVC2 svc2_1(sl1,sc2,sl2,sc1_1);

  Dune::PDELab::TypeTree::applyToTreePair(svc1_2,svc2_1,PairPrinter());

  typedef Dune::PDELab::TypeTree::TreeInfo<SVC2> TI;

  // test TreeInfo
  dune_static_assert(TI::depth == 4 && TI::nodeCount == 14 && TI::leafCount == 10,
                     "TreeInfo yields wrong information");

  std::cout << "depth: " << TI::depth << std::endl
            << "nodes: " << TI::nodeCount << std::endl
            << "leafs: " << TI::leafCount << std::endl;

  dune_static_assert((Dune::PDELab::TypeTree::AccumulateValue<SVC2,
                                                              NodeCountingFunctor,
                                                              Dune::PDELab::TypeTree::plus<std::size_t>,
                                                              0>::result == TI::nodeCount),
                     "Error in AccumulateValue");

  dune_static_assert((Dune::PDELab::TypeTree::AccumulateValue<SVC2,
                                                              LeafCountingFunctor,
                                                              Dune::PDELab::TypeTree::plus<std::size_t>,
                                                              0>::result == TI::leafCount),
                     "Error in AccumulateValue");

  dune_static_assert((Dune::PDELab::TypeTree::AccumulateValue<SVC2,
                                                              DepthFunctor,
                                                              Dune::PDELab::TypeTree::max<std::size_t>,
                                                              0>::result == TI::depth),
                     "Error in AccumulateValue");

#endif

#endif

#endif

  return 0;
}

#endif
