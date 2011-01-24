#include "config.h"

#include "typetreetestswitch.hh"

#if TEST_TYPETREE_INVALID

int main()
{
  return 0;
}

#else

#include <dune/pdelab/common/leafnode.hh>
#include <dune/pdelab/common/powernode.hh>
#include <dune/pdelab/common/compositenode.hh>
#include <dune/pdelab/common/traversal.hh>

#if HAVE_VARIADIC_TEMPLATES
#include <dune/pdelab/common/variadiccompositenode.hh>
#endif

#include <iostream>

struct Counter
{
  Counter()
    : _id(_ids++)
  {
    std::cout << "Constructed id = " << id() << std::endl;
  }

  Counter(const Counter& rhs)
    : _id(_ids++)
  {
    std::cout << "Copy-Constructed id = " << id() << " from id = " << rhs.id() << std::endl;
  }

  Counter(Counter&& rhs)
    : _id(rhs._id)
  {
    rhs._id = -1;
    std::cout << "Move-Constructed id = " << id() << " from id = " << rhs.id() << std::endl;
  }

  ~Counter()
  {
    std::cout << "Destructed id = " << id() << std::endl;
  }

  Counter& operator=(const Counter& rhs)
  {
    std::cout << "Assigned id = " << id() << " from id = " << rhs.id() << std::endl;
    return *this;
  }

  Counter& operator=(Counter&& rhs)
  {
    std::cout << "Move-Assigned id = " << id() << " from id = " << rhs.id() << std::endl;
    return *this;
  }

  int id() const
  {
    return _id;
  }

  int _id;
  static int _ids;
};

int Counter::_ids = 0;

struct SimpleLeaf
  : public Dune::PDELab::TypeTree::LeafNode
  , public Counter
{

  static const char* name()
  {
    return "SimpleLeaf";
  }

  SimpleLeaf() {}

  SimpleLeaf(SimpleLeaf&& rhs)
    : Dune::PDELab::TypeTree::LeafNode(std::move(rhs))
    , Counter(std::move(rhs))
  {
    std::cout << "move ctor" << std::endl;
  }
};

template<typename T, typename U>
struct check_args
{
  typedef typename std::enable_if<std::is_same<typename std::remove_reference<T>::type,U>::value,T>::type type;
};

template<typename T, std::size_t k>
struct SimplePower
  : public Dune::PDELab::TypeTree::PowerNode<T,k>
  , public Counter
{
  static const char* name()
  {
    return "SimplePower";
  }

  typedef Dune::PDELab::TypeTree::PowerNode<T,k> BaseT;

  SimplePower() {}

  SimplePower(T& c, bool copy)
    : BaseT(c,copy)
  {}

#if HAVE_VARIADIC_TEMPLATES && HAVE_RVALUE_REFERENCES

  template<typename C1, typename C2, typename... Children>
  SimplePower(C1&& c1, C2&& c2, Children&&... children)
    : BaseT(std::forward<C1>(c1),std::forward<C2>(c2),std::forward<Children>(children)...)
  {}

#else

  SimplePower(T& c0, T& c1)
    : BaseT(c0,c1)
  {}

  SimplePower(T& c0, T& c1, T& c2)
    : BaseT(c0,c1,c2)
  {}

#endif

};

#if HAVE_VARIADIC_TEMPLATES && HAVE_RVALUE_REFERENCES

template<typename... Children>
struct SimpleVariadicComposite
  : public Dune::PDELab::TypeTree::VariadicCompositeNode<Children...>
  , public Counter
{

  static const char* name()
  {
    return "SimpleVariadicComposite";
  }

  typedef Dune::PDELab::TypeTree::VariadicCompositeNode<Children...> BaseT;

#if HAVE_VARIADIC_CONSTRUCTOR_SFINAE

  template<typename... Args, typename = typename std::enable_if<(sizeof...(Args) == BaseT::CHILDREN)>::type>
  SimpleVariadicComposite(Args&&... args)
    : BaseT(std::forward<Args>(args)...)
  {}

#else

  SimpleVariadicComposite(Children&... children)
    : BaseT(children...)
  {}

#endif

};

#endif

template<typename C1, typename C2 = Dune::PDELab::TypeTree::EmptyNode, typename C3 = Dune::PDELab::TypeTree::EmptyNode, typename C4 = Dune::PDELab::TypeTree::EmptyNode>
struct SimpleComposite
  : public Dune::PDELab::TypeTree::CompositeNode<C1,C2,C3,C4>
  , public Counter
{

  static const char* name()
  {
    return "SimpleComposite";
  }

  typedef Dune::PDELab::TypeTree::CompositeNode<C1,C2,C3,C4> BaseT;

  SimpleComposite(C1& c1,
                  typename Dune::PDELab::TypeTree::OptionalChild<C2>::type c2 = Dune::PDELab::TypeTree::OptionalChild<C2>::default_value(),
                  typename Dune::PDELab::TypeTree::OptionalChild<C3>::type c3 = Dune::PDELab::TypeTree::OptionalChild<C3>::default_value(),
                  typename Dune::PDELab::TypeTree::OptionalChild<C4>::type c4 = Dune::PDELab::TypeTree::OptionalChild<C4>::default_value())
    : BaseT(c1,c2,c3,c4)
  {}

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
    std::cout << t.name() << " " << t.id() << std::endl;
  }
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

  SimpleLeaf sl2;
  SP1 sp1_2(sl2,false);

  Dune::PDELab::TypeTree::applyToTree(sp1_1,TreePrinter());

  typedef SimpleComposite<SimpleLeaf,SP1,SimpleLeaf> SC1;
  SC1 sc1_1(sl1,sp1_2);
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

  return 0;
}

#endif
