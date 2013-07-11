struct TargetLeaf
  : public Dune::PDELab::TypeTree::LeafNode
{

  template<typename Transformation>
  TargetLeaf(const SimpleLeaf& sl, const Transformation& t)
    : s(Dune::stackobject_to_shared_ptr(sl))
  {}

  template<typename Transformation>
  TargetLeaf(Dune::shared_ptr<const SimpleLeaf> sl, const Transformation& t)
    : s(sl)
  {}

  Dune::shared_ptr<const SimpleLeaf> s;

  const char* name() const
  {
    return "TargetLeaf";
  }

  int id() const
  {
    return s->id();
  }

};

template<typename S, typename T, std::size_t k>
struct TargetPower
  : public Dune::PDELab::TypeTree::PowerNode<T,k>
{

  template<typename Transformation>
  TargetPower(const S& sc, const Transformation& t, const Dune::array<Dune::shared_ptr<T>,k>& children)
    : Dune::PDELab::TypeTree::PowerNode<T,k>(children)
    , s(Dune::stackobject_to_shared_ptr(sc))
  {}

  template<typename Transformation>
  TargetPower(Dune::shared_ptr<const S> sc, const Transformation& t, const Dune::array<Dune::shared_ptr<T>,k>& children)
    : Dune::PDELab::TypeTree::PowerNode<T,k>(children)
    , s(sc)
  {}

  Dune::shared_ptr<const S> s;

  const char* name() const
  {
    return "TargetPower";
  }

  int id() const
  {
    return s->id();
  }


};

template<typename S, typename... Children>
struct TargetVariadicComposite
  : public Dune::PDELab::TypeTree::VariadicCompositeNode<Children...>
{

  template<typename Transformation>
  TargetVariadicComposite(const S& sc, const Transformation& t, Dune::shared_ptr<Children>... children)
    : Dune::PDELab::TypeTree::VariadicCompositeNode<Children...>(children...)
    , s(Dune::stackobject_to_shared_ptr(sc))
  {}

  template<typename Transformation>
  TargetVariadicComposite(Dune::shared_ptr<const S> sc, const Transformation& t, Dune::shared_ptr<Children>... children)
    : Dune::PDELab::TypeTree::VariadicCompositeNode<Children...>(children...)
    , s(sc)
  {}

  Dune::shared_ptr<const S> s;

  const char* name() const
  {
    return "TargetVariadicComposite";
  }

  int id() const
  {
    return s->id();
  }


};


template<typename S,
         typename C0,
         typename C1 = Dune::PDELab::TypeTree::EmptyNode,
         typename C2 = Dune::PDELab::TypeTree::EmptyNode,
         typename C3 = Dune::PDELab::TypeTree::EmptyNode,
         typename C4 = Dune::PDELab::TypeTree::EmptyNode,
         typename C5 = Dune::PDELab::TypeTree::EmptyNode,
         typename C6 = Dune::PDELab::TypeTree::EmptyNode,
         typename C7 = Dune::PDELab::TypeTree::EmptyNode,
         typename C8 = Dune::PDELab::TypeTree::EmptyNode,
         typename C9 = Dune::PDELab::TypeTree::EmptyNode>
struct TargetComposite
  : public Dune::PDELab::TypeTree::CompositeNode<C0,C1,C2,C3,C4,C5,C6,C7,C8,C9>
{

  typedef Dune::PDELab::TypeTree::CompositeNode<C0,C1,C2,C3,C4,C5,C6,C7,C8,C9> BaseT;

  template<typename Transformation>
  TargetComposite(const S& sc, const Transformation& t,
                  Dune::shared_ptr<C0> c0,
                  Dune::shared_ptr< C1> c1,
                  Dune::shared_ptr< C2> c2,
                  Dune::shared_ptr< C3> c3,
                  Dune::shared_ptr< C4> c4,
                  Dune::shared_ptr< C5> c5,
                  Dune::shared_ptr< C6> c6,
                  Dune::shared_ptr< C7> c7,
                  Dune::shared_ptr< C8> c8,
                  Dune::shared_ptr< C9> c9)
    : BaseT(c0,c1,c2,c3,c4,c5,c6,c7,c8,c9)
    , s(Dune::stackobject_to_shared_ptr(sc))
  {}

  template<typename Transformation>
  TargetComposite(Dune::shared_ptr<const S> sc, const Transformation& t,
                  Dune::shared_ptr< C0> c0,
                  Dune::shared_ptr< C1> c1,
                  Dune::shared_ptr< C2> c2,
                  Dune::shared_ptr< C3> c3,
                  Dune::shared_ptr< C4> c4,
                  Dune::shared_ptr< C5> c5,
                  Dune::shared_ptr< C6> c6,
                  Dune::shared_ptr< C7> c7,
                  Dune::shared_ptr< C8> c8,
                  Dune::shared_ptr< C9> c9)
    : BaseT(c0,c1,c2,c3,c4,c5,c6,c7,c8,c9)
    , s(sc)
  {}

  Dune::shared_ptr<const S> s;

  const char* name() const
  {
    return "TargetComposite";
  }

  int id() const
  {
    return s->id();
  }


};

struct TestTransformation {};

// register leaf node
template<typename SL>
Dune::PDELab::TypeTree::GenericLeafNodeTransformation<SimpleLeaf,TestTransformation,TargetLeaf>
registerNodeTransformation(SL* sl, TestTransformation* t, SimpleLeafTag* tag);

template<typename SP>
Dune::PDELab::TypeTree::GenericPowerNodeTransformation<SP,TestTransformation,TargetPower>
registerNodeTransformation(SP* sp, TestTransformation* t, SimplePowerTag* tag);

template<typename SVC>
Dune::PDELab::TypeTree::GenericVariadicCompositeNodeTransformation<SVC,TestTransformation,TargetVariadicComposite>
registerNodeTransformation(SVC* svc, TestTransformation* t, SimpleVariadicCompositeTag* tag);

template<typename SC>
Dune::PDELab::TypeTree::GenericCompositeNodeTransformation<SC,TestTransformation,TargetComposite>
registerNodeTransformation(SC* svc, TestTransformation* t, SimpleCompositeTag* tag);
