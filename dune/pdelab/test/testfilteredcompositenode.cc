#include "config.h"

#include "typetreetestswitch.hh"

#if TEST_TYPETREE_INVALID

int main()
{
  return 0;
}

#else

#include <dune/pdelab/common/typetree/filteredcompositenode.hh>

#include "typetreetestutility.hh"
#include "typetreetargetnodes.hh"

#include <type_traits>

struct LeafFilter
  : public Dune::PDELab::TypeTree::DefaultFilter
{
  template<typename T, std::size_t new_k, std::size_t old_k>
  struct apply
  {
    static const bool value = std::is_same<typename T::NodeTag,Dune::PDELab::TypeTree::LeafNodeTag>::value;
  };
};


template<typename Node, typename Filter>
class SimpleFilteredNode
  : public Dune::PDELab::TypeTree::FilteredCompositeNode<Node,Filter>
  , public Counter
{

  typedef Dune::PDELab::TypeTree::FilteredCompositeNode<Node,Filter> BaseT;

public:

  typedef SimpleVariadicCompositeTag ImplementationTag;

  static const char* name()
  {
    return "SimpleFilteredNode";
  }

  SimpleFilteredNode(Node& node)
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

  typedef Dune::PDELab::TypeTree::TreeInfo<SVC2> TI;

  typedef SimpleFilteredNode<SVC2,LeafFilter> FCN1;
  typedef Dune::PDELab::TypeTree::TreeInfo<FCN1> FCN1_TI;

  FCN1 fcn1_1(svc2_1);

  typedef SimpleFilteredNode<SVC2,Dune::PDELab::TypeTree::IndexFilter<1,2,3> > FCN2;
  typedef Dune::PDELab::TypeTree::TreeInfo<FCN2> FCN2_TI;

  FCN2 fcn2_1(svc2_1);

  typedef SimpleFilteredNode<SVC2,Dune::PDELab::TypeTree::IndexFilter<3,1,0,4,1,3,3,3,3> > FCN3;
  typedef Dune::PDELab::TypeTree::TreeInfo<FCN3> FCN3_TI;

  FCN3 fcn3_1(svc2_1);


  // test TreeInfo
  //dune_static_assert(TI::depth == 4 && TI::nodeCount == 14 && TI::leafCount == 10,
  //                   "TreeInfo yields wrong information");

  Dune::PDELab::TypeTree::applyToTree(svc2_1,TreePrinter());
  std::cout << "depth: " << TI::depth << std::endl
            << "nodes: " << TI::nodeCount << std::endl
            << "leafs: " << TI::leafCount << std::endl;

  std::cout << std::endl << "leaf filter:" << std::endl;
  Dune::PDELab::TypeTree::applyToTree(fcn1_1,TreePrinter());
  std::cout << "depth: " << FCN1_TI::depth << std::endl
            << "nodes: " << FCN1_TI::nodeCount << std::endl
            << "leafs: " << FCN1_TI::leafCount << std::endl;

  std::cout << std::endl << "index filter (1,2,3):" << std::endl;
  Dune::PDELab::TypeTree::applyToTree(fcn2_1,TreePrinter());
  std::cout << "depth: " << FCN2_TI::depth << std::endl
            << "nodes: " << FCN2_TI::nodeCount << std::endl
            << "leafs: " << FCN2_TI::leafCount << std::endl;

  std::cout << std::endl << "index filter (3,1,0,4):" << std::endl;
  Dune::PDELab::TypeTree::applyToTree(fcn3_1,TreePrinter());
  std::cout << "depth: " << FCN3_TI::depth << std::endl
            << "nodes: " << FCN3_TI::nodeCount << std::endl
            << "leafs: " << FCN3_TI::leafCount << std::endl;

  typedef Dune::PDELab::TypeTree::TransformTree<FCN3,TestTransformation> Transformation;

  typedef typename Transformation::Type TFCN3;

  TFCN3 tfcn3_1 = Transformation::transform(fcn3_1);

  typedef Dune::PDELab::TypeTree::TreeInfo<TFCN3> TFCN3_TI;

  dune_static_assert(TFCN3_TI::depth == FCN3_TI::depth, "error in transformation with filtered node");

  dune_static_assert(TFCN3_TI::nodeCount == FCN3_TI::nodeCount, "error in transformation with filtered node");

  dune_static_assert(TFCN3_TI::leafCount == FCN3_TI::leafCount, "error in transformation with filtered node");

  Dune::PDELab::TypeTree::applyToTree(tfcn3_1,TreePrinter());

  return 0;
}

#endif
