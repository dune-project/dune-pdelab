#include "config.h"

#include "typetreetestswitch.hh"

#if TEST_TYPETREE_INVALID

int main()
{
  return 0;
}

#else

#include <dune/pdelab/common/typetree.hh>

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
};

template<typename T, typename U>
struct check_args
{
  typedef typename std::enable_if<std::is_same<typename std::remove_reference<T>::type,U>::value,T>::type type;
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


namespace Dune {
  namespace PDELab {
    namespace TypeTree {

      template<typename S, typename T, typename transformed>
      struct WrappingLeafTransformation
      {

        typedef transformed transformed_type;
        typedef Dune::shared_ptr<transformed> transformed_storage_type;

        static transformed_type transform(const S& s, const T& t)
        {
          return transformed_type(s);
        }

        static transformed_storage_type transform_storage(Dune::shared_ptr<const S> s, const T& t)
        {
          return Dune::make_shared<transformed_type>(*s);
        }

      };


      template<typename S, typename T, template<typename U> class transformation>
      struct WrappingPowerTransformation
      {

        template<typename TC>
        struct result
        {
          typedef typename transformation<TC>::type transformed_type;
          typedef Dune::shared_ptr<transformed_type> transformed_storage_type;
        };

        template<typename TC>
        static typename result<TC>::transformed_type transform(const S& s, const T& t, const Dune::array<Dune::shared_ptr<TC>,result<TC>::transformed_type::CHILDREN>& children)
        {
          return typename result<TC>::transformed_type(s,children);
        }

        template<typename TC>
        static typename result<TC>::transformed_storage_type transform_storage(Dune::shared_ptr<const S> s, const T& t, const Dune::array<Dune::shared_ptr<TC>,result<TC>::transformed_type::CHILDREN>& children)
        {
          return Dune::make_shared<typename result<TC>::transformed_type>(*s,children);
        }

      };


      template<typename S, typename T, template<typename... U> class transformation>
      struct WrappingVariadicCompositeTransformation
      {

        template<typename... TC>
        struct result
        {
          typedef typename transformation<TC...>::type transformed_type;
          typedef Dune::shared_ptr<transformed_type> transformed_storage_type;
        };

        template<typename... TC>
        static typename result<TC...>::transformed_type transform(const S& s, const T& t, std::shared_ptr<TC>... children)
        {
          return typename result<TC...>::transformed_type(s,children...);
        }

        template<typename... TC>
        static typename result<TC...>::transformed_storage_type transform_storage(Dune::shared_ptr<const S> s, const T& t, std::shared_ptr<TC>... children)
        {
          return std::make_shared<typename result<TC...>::transformed_type>(*s,children...);
        }

      };

      template<typename S, template<typename...> class TransformedNode>
      struct DefaultWrappingVariadicCompositeTransformationHelper
      {
        template<typename... TC>
        struct result
        {
          typedef TransformedNode<S,TC...> type;
        };
      };

      template<typename S, typename T, template<typename...> class TransformedNode>
      struct DefaultWrappingVariadicCompositeTransformation
        : public WrappingVariadicCompositeTransformation<S,T,DefaultWrappingVariadicCompositeTransformationHelper<S,TransformedNode>::template result>
      {};

    }
  }
}


struct SimpleCompositeTag {};

template<typename C1, typename C2 = Dune::PDELab::TypeTree::EmptyNode, typename C3 = Dune::PDELab::TypeTree::EmptyNode, typename C4 = Dune::PDELab::TypeTree::EmptyNode>
struct SimpleComposite
  : public Dune::PDELab::TypeTree::CompositeNode<C1,C2,C3,C4>
  , public Counter
{

  typedef SimpleCompositeTag ImplementationTag;

  static const char* name()
  {
    return "SimpleComposite";
  }

  typedef Dune::PDELab::TypeTree::CompositeNode<C1,C2,C3,C4> BaseT;

  SimpleComposite(C1& c1,
                  typename Dune::PDELab::TypeTree::OptionalChild<C2>::type c2 = typename Dune::PDELab::TypeTree::OptionalChild<C2>::type(),
                  typename Dune::PDELab::TypeTree::OptionalChild<C3>::type c3 = typename Dune::PDELab::TypeTree::OptionalChild<C3>::type(),
                  typename Dune::PDELab::TypeTree::OptionalChild<C4>::type c4 = typename Dune::PDELab::TypeTree::OptionalChild<C4>::type())
    : BaseT(c1,c2,c3,c4)
  {}

};


struct TargetLeaf
  : public Dune::PDELab::TypeTree::LeafNode
{

  TargetLeaf(const SimpleLeaf& sl)
    : s(sl)
  {}

  const SimpleLeaf& s;

  const char* name() const
  {
    return "TargetLeaf";
  }

  int id() const
  {
    return s.id();
  }

};

template<typename S, typename T, std::size_t k>
struct TargetPower
  : public Dune::PDELab::TypeTree::PowerNode<T,k>
{

  TargetPower(const S& sc, const Dune::array<Dune::shared_ptr<T>,k>& children)
    : Dune::PDELab::TypeTree::PowerNode<T,k>(children)
    , s(sc)
  {}

  const S& s;

  const char* name() const
  {
    return "TargetPower";
  }

  int id() const
  {
    return s.id();
  }


};

template<typename S, typename... Children>
struct TargetVariadicComposite
  : public Dune::PDELab::TypeTree::VariadicCompositeNode<Children...>
{

  TargetVariadicComposite(const S& sc, Dune::shared_ptr<Children>... children)
    : Dune::PDELab::TypeTree::VariadicCompositeNode<Children...>(children...)
    , s(sc)
  {}

  const S& s;

  const char* name() const
  {
    return "TargetVariadicComposite";
  }

  int id() const
  {
    return s.id();
  }


};

struct TestTransformation {};

// register leaf node
template<typename SL>
Dune::PDELab::TypeTree::WrappingLeafTransformation<SimpleLeaf,TestTransformation,TargetLeaf>
transformNode(const SL& sl, const TestTransformation& t, SimpleLeafTag tag)
{
  return Dune::PDELab::TypeTree::WrappingLeafTransformation<SimpleLeaf,TestTransformation,TargetLeaf>();
}


// register power node
template<typename S>
struct TransformSimplePower
{
  template<typename TC>
  struct result
  {
    typedef TargetPower<S,TC,S::CHILDREN> type;
  };
};

template<typename SP>
Dune::PDELab::TypeTree::WrappingPowerTransformation<SP,TestTransformation,TransformSimplePower<SP>::template result>
transformNode(const SP& sp, const TestTransformation& t, SimplePowerTag tag)
{
  return Dune::PDELab::TypeTree::WrappingPowerTransformation<SP,TestTransformation,TransformSimplePower<SP>::template result>();
}


// register variadic composite node
/*
template<typename S>
struct TransformSimpleVariadicComposite
{
  template<typename... TC>
  struct result
  {
    typedef TargetVariadicComposite<S,TC...> type;
  };
};

template<typename SVC>
Dune::PDELab::TypeTree::WrappingVariadicCompositeTransformation<SVC,TestTransformation,TransformSimpleVariadicComposite<SVC>::template result>
transformNode(const SVC& svc, const TestTransformation& t, SimpleVariadicCompositeTag tag)
{
  return Dune::PDELab::TypeTree::WrappingVariadicCompositeTransformation<SVC,TestTransformation,TransformSimpleVariadicComposite<SVC>::template result>();
}
*/

template<typename SVC>
Dune::PDELab::TypeTree::DefaultWrappingVariadicCompositeTransformation<SVC,TestTransformation,TargetVariadicComposite>
transformNode(const SVC& svc, const TestTransformation& t, SimpleVariadicCompositeTag tag)
{
  return Dune::PDELab::TypeTree::DefaultWrappingVariadicCompositeTransformation<SVC,TestTransformation,TargetVariadicComposite>();
}


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
