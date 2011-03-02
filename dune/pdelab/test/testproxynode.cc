#include "config.h"

#include "typetreetestswitch.hh"

#if TEST_TYPETREE_INVALID

int main()
{
  return 0;
}

#else

#include <dune/pdelab/common/typetree/proxynode.hh>


#include "typetreetestutility.hh"
#include "typetreetargetnodes.hh"

#include <type_traits>


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

  typedef SimpleVariadicComposite<SimpleLeaf,SP1,SimpleLeaf,SC1> SVC1;
  SVC1 svc1_1(sl1,sp1_1,sl2,sc1_1);

  SP1 sp1_3(SimpleLeaf(),SimpleLeaf(),sl1);

  SVC1 svc1_2(SimpleLeaf(),SP1(sp1_2),sl2,const_cast<const SC1&>(sc1_1));

  typedef SimpleVariadicComposite<SimpleLeaf,SVC1,SimpleLeaf,SP1,SC1> SVC2;
  SVC2 svc2_1(sl1,svc1_2,sl2,sp1_3,sc1_1);

  typedef SimpleProxy<SimpleLeaf> PSL;
  PSL psl(sl1);
  Dune::PDELab::TypeTree::applyToTree(psl,TreePrinter());
  typedef Dune::PDELab::TypeTree::TreeInfo<SimpleLeaf> TI_SL;
  typedef Dune::PDELab::TypeTree::TreeInfo<PSL> TI_PSL;
  dune_static_assert(TI_SL::depth == TI_PSL::depth, "Proxy node has wrong depth");
  dune_static_assert(TI_SL::nodeCount == TI_PSL::nodeCount, "Proxy node has wrong node count");
  dune_static_assert(TI_SL::leafCount == TI_PSL::leafCount, "Proxy node has wrong leaf count");

  typedef SimpleProxy<SP1> PSP1;
  PSP1 psp1(sp1_1);
  Dune::PDELab::TypeTree::applyToTree(psp1,TreePrinter());
  typedef Dune::PDELab::TypeTree::TreeInfo<SP1> TI_SP1;
  typedef Dune::PDELab::TypeTree::TreeInfo<PSP1> TI_PSP1;
  dune_static_assert(TI_SP1::depth == TI_PSP1::depth, "Proxy node has wrong depth");
  dune_static_assert(TI_SP1::nodeCount == TI_PSP1::nodeCount, "Proxy node has wrong node count");
  dune_static_assert(TI_SP1::leafCount == TI_PSP1::leafCount, "Proxy node has wrong leaf count");

  typedef SimpleProxy<SC1> PSC1;
  PSC1 psc1(sc1_1);
  Dune::PDELab::TypeTree::applyToTree(psc1,TreePrinter());
  typedef Dune::PDELab::TypeTree::TreeInfo<SC1> TI_SC1;
  typedef Dune::PDELab::TypeTree::TreeInfo<PSC1> TI_PSC1;
  dune_static_assert(TI_SC1::depth == TI_PSC1::depth, "Proxy node has wrong depth");
  dune_static_assert(TI_SC1::nodeCount == TI_PSC1::nodeCount, "Proxy node has wrong node count");
  dune_static_assert(TI_SC1::leafCount == TI_PSC1::leafCount, "Proxy node has wrong leaf count");

  typedef SimpleProxy<SVC2> PSVC2;
  PSVC2 psvc2(svc2_1);
  Dune::PDELab::TypeTree::applyToTree(psvc2,TreePrinter());
  typedef Dune::PDELab::TypeTree::TreeInfo<SVC2> TI_SVC2;
  typedef Dune::PDELab::TypeTree::TreeInfo<PSVC2> TI_PSVC2;
  dune_static_assert(TI_SVC2::depth == TI_PSVC2::depth, "Proxy node has wrong depth");
  dune_static_assert(TI_SVC2::nodeCount == TI_PSVC2::nodeCount, "Proxy node has wrong node count");
  dune_static_assert(TI_SVC2::leafCount == TI_PSVC2::leafCount, "Proxy node has wrong leaf count");

  return 0;
}

#endif
