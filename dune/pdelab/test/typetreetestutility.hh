
#include <dune/pdelab/common/typetree.hh>
#include <dune/pdelab/common/typetree/pairtraversal.hh>

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

struct SimpleLeafTag {};

struct SimpleLeaf
  : public Dune::PDELab::TypeTree::LeafNode
  , public Counter
{

  typedef SimpleLeafTag ImplementationTag;

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

struct SimpleLeafDerived
  : public SimpleLeaf
{

  static const char* name()
  {
    return "SimpleLeafDerived";
  }

};

struct SimplePowerTag {};

template<typename T, std::size_t k>
struct SimplePower
  : public Dune::PDELab::TypeTree::PowerNode<T,k>
  , public Counter
{

  typedef SimplePowerTag ImplementationTag;

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

struct SimpleVariadicCompositeTag {};

template<typename... Children>
struct SimpleVariadicComposite
  : public Dune::PDELab::TypeTree::VariadicCompositeNode<Children...>
  , public Counter
{

  typedef SimpleVariadicCompositeTag ImplementationTag;

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
  : public Dune::PDELab::TypeTree::TreeVisitor
  , public Dune::PDELab::TypeTree::DynamicTraversal
{

  template<typename T, typename TreePath>
  void leaf(const T& t, TreePath treePath) const
  {
    pre(t,treePath);
  }

  template<typename T, typename TreePath>
  void pre(const T& t, TreePath treePath) const
  {
    for (std::size_t i = 0; i < treePath.size(); ++i)
      std::cout << "  ";
    std::cout << t.name() << " " << t.id() << std::endl;
  }
};




struct PairPrinter
  : public Dune::PDELab::TypeTree::TreePairVisitor
  , public Dune::PDELab::TypeTree::DynamicTraversal
{

  template<typename T1, typename T2, typename TreePath>
  void leaf(const T1& t1, const T2& t2, TreePath treePath) const
  {
    pre(t1,t2,treePath);
  }

  template<typename T1, typename T2, typename TreePath>
  void pre(const T1& t1, const T2& t2, TreePath treePath) const
  {
    for (std::size_t i = 0; i < treePath.size(); ++i)
      std::cout << "  ";
    std::cout << t1.name() << " " << t1.id() << "      " << t2.name() << " " << t2.id() << std::endl;
  }
};
