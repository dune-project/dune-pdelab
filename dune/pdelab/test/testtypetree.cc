#include "config.h"
#include <dune/pdelab/common/leafnode.hh>
#include <dune/pdelab/common/powernode.hh>
//#include <dune/pdelab/common/compositenode.hh>
#include <dune/pdelab/common/traversal.hh>

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
    : _id(_ids++)
  {
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

  const int _id;
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

private:
  //SimplePower(const SimplePower&);
public:

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
  TreePrinter treePrinter;
  SimpleLeaf sl1;
  typedef SimplePower<SimpleLeaf,3> SP1;
  SimplePower<SimpleLeaf,3> sp1;
  sp1.setChild(0,sl1);
  sp1.setChild(1,sl1);
  sp1.setChild(2,sl1);
  Dune::PDELab::TypeTree::applyToTree(sp1,TreePrinter());
#if HAVE_RVALUE_REFERENCES
  SimplePower<SimplePower<SimpleLeaf,3>,2> sp2(sp1,const_cast<const SP1&>(sp1));
  SimplePower<SimplePower<SimpleLeaf,3>,2> sp3(SP1(SimpleLeaf(),SimpleLeaf(),sl1),SP1(sl1,sl1,SimpleLeaf()));
  Dune::PDELab::TypeTree::applyToTree(sp3,treePrinter);
#else
  SimplePower<SimplePower<SimpleLeaf,3>,2> sp2(sp1,sp1);
  Dune::PDELab::TypeTree::applyToTree(sp2,treePrinter);
#endif
  Dune::PDELab::TypeTree::applyToTree(const_cast<const SimplePower<SimplePower<SimpleLeaf,3>,2>&>(sp2),TreePrinter());
  Dune::PDELab::TypeTree::applyToTree(const_cast<const SP1&>(sp1),treePrinter);
  return 0;
}
