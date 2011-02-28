#include "config.h"

#include "typetreetestswitch.hh"

#if TEST_TYPETREE_INVALID

int main()
{
  return 0;
}

#else

#include "typetreetestutility.hh"
#include "typetreetargetnodes.hh"

int main(int argc, char** argv)
{

  // basic tests

  // leaf node
  TreePrinter treePrinter;
  SimpleLeaf sl1;

  Dune::PDELab::TypeTree::applyToTree(sl1,treePrinter);

  Dune::PDELab::TypeTree::TransformTree<SimpleLeaf,TestTransformation>::transformed_type tl1 =
    Dune::PDELab::TypeTree::TransformTree<SimpleLeaf,TestTransformation>::transform(sl1,TestTransformation());

  typedef SimplePower<SimpleLeaf,3> SP1;
  SP1 sp1_1;
  sp1_1.setChild(0,sl1);
  sp1_1.setChild(1,sl1);
  sp1_1.setChild(2,sl1);

  SimpleLeaf sl2;
  SP1 sp1_2(sl2,false);

  typedef SimpleVariadicComposite<SimpleLeafDerived,SP1,SimpleLeaf> SVC1;

  SVC1 svc1_1(SimpleLeafDerived(),sp1_2,sl1);

  Dune::PDELab::TypeTree::applyToTree(sp1_1,TreePrinter());

  TestTransformation trafo;

  Dune::PDELab::TypeTree::TransformTree<SP1,TestTransformation>::transformed_type tp1_1 =
    Dune::PDELab::TypeTree::TransformTree<SP1,TestTransformation>::transform(sp1_1,trafo);

  Dune::PDELab::TypeTree::TransformTree<SVC1,TestTransformation>::transformed_type tvc1_1 =
    Dune::PDELab::TypeTree::TransformTree<SVC1,TestTransformation>::transform(svc1_1,TestTransformation());

  Dune::PDELab::TypeTree::applyToTree(tvc1_1,TreePrinter());

  /*
  typedef SimpleComposite<SimpleLeaf,SP1,SimpleLeaf> SC1;
  SC1 sc1_1(sl1,sp1_2,sl2);
  Dune::PDELab::TypeTree::applyToTree(const_cast<const SC1&>(sc1_1),treePrinter);

#if HAVE_VARIADIC_TEMPLATES

#if HAVE_RVALUE_REFERENCES

  typedef SimpleVariadicComposite<SimpleLeaf,SP1,SimpleLeaf,SC1> SVC1;
  SVC1 svc1_1(sl1,sp1_1,sl2,sc1_1);
  Dune::PDELab::TypeTree::applyToTree(svc1_1,treePrinter);

  SP1 sp1_3(SimpleLeaf(),SimpleLeaf(),sl1);
  Dune::PDELab::TypeTree::applyToTree(sp1_3,TreePrinter());

#if HAVE_VARIADIC_CONSTRUCTOR_SFINAE

  SVC1 svc1_2(SimpleLeaf(),SP1(sp1_2),sl2,const_cast<const SC1&>(sc1_1));
  Dune::PDELab::TypeTree::applyToTree(svc1_2,TreePrinter());

#endif

#endif

#endif
  */
  return 0;
}

#endif
