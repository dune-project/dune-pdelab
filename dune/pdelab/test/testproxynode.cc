#include "config.h"

#include <dune/pdelab/common/typetree/proxynode.hh>


#include "typetreetestutility.hh"

template<typename Node>
class SimpleProxy
  : public Dune::PDELab::TypeTree::ProxyNode<Node>
{

  typedef Dune::PDELab::TypeTree::ProxyNode<Node> BaseT;

public:

  static const char* name()
  {
    return "SimpleProxy";
  }

  int id() const
  {
    return this->proxiedNode().id();
  }

  SimpleProxy(Node& node)
    : BaseT(node)
  {}

};

template<typename Node>
void testProxyNode(Node& node)
{
  typedef SimpleProxy<Node> ProxyNode;
  ProxyNode proxyNode(node);
  Dune::PDELab::TypeTree::applyToTree(proxyNode,TreePrinter());
  typedef Dune::PDELab::TypeTree::TreeInfo<Node> TI_Node;
  typedef Dune::PDELab::TypeTree::TreeInfo<ProxyNode> TI_ProxyNode;
  dune_static_assert(TI_Node::depth == TI_ProxyNode::depth, "Proxy node has wrong depth");
  dune_static_assert(TI_Node::nodeCount == TI_ProxyNode::nodeCount, "Proxy node has wrong node count");
  dune_static_assert(TI_Node::leafCount == TI_ProxyNode::leafCount, "Proxy node has wrong leaf count");
}


int main(int argc, char** argv)
{

  // basic tests

  // leaf node
  SimpleLeaf sl1;

  typedef SimplePower<SimpleLeaf,3> SP1;
  SP1 sp1_1;
  sp1_1.setChild(0,sl1);
  sp1_1.setChild(1,sl1);
  sp1_1.setChild(2,sl1);

  SimpleLeaf sl2;
  SP1 sp1_2(sl2,false);

  Dune::PDELab::TypeTree::applyToTree(sp1_1,TreePrinter());

  typedef SimpleComposite<SimpleLeaf,SP1,SimpleLeaf> SC1;
  SC1 sc1_1(sl1,sp1_2,sl2);

  typedef SimpleComposite<SimpleLeaf,SimpleLeaf,SimpleLeaf> SC2;
  SC2 sc2(sl1,sl1,sl1);

  testProxyNode(sl1);
  testProxyNode(sp1_1);
  testProxyNode(sc1_1);

#if HAVE_VARIADIC_TEMPLATES && HAVE_VARIADIC_CONSTRUCTOR_SFINAE

  typedef SimpleVariadicComposite<SimpleLeaf,SP1,SimpleLeaf,SC1> SVC1;
  SVC1 svc1_1(sl1,sp1_1,sl2,sc1_1);

  SP1 sp1_3(SimpleLeaf(),SimpleLeaf(),sl1);

  SVC1 svc1_2(SimpleLeaf(),SP1(sp1_2),sl2,const_cast<const SC1&>(sc1_1));

  typedef SimpleVariadicComposite<SimpleLeaf,SVC1,SimpleLeaf,SP1,SC1> SVC2;
  SVC2 svc2_1(sl1,svc1_2,sl2,sp1_3,sc1_1);

  testProxyNode(svc2_1);

#endif // HAVE_VARIADIC_TEMPLATES

  return 0;
}
