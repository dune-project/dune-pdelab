//-*- tab-width: 4; c-basic-offset: 2; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_FUNCTION_OLADINERFACEADAPTER_HH
#define DUNE_PDELAB_FUNCTION_OLADINERFACEADAPTER_HH

#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/common/typetraits.hh>
#include <dune/functions/common/signature.hh>
#include <dune/functions/common/defaultderivativetraits.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>

namespace Dune {
namespace PDELab {

#warning should this go to dune-functions?
  struct DifferentiableFunctionBaseTag {};
  struct DifferentiableFunctionLocalViewBaseTag {};

  struct DifferentiableFunctionTag : public DifferentiableFunctionBaseTag {};
  struct DifferentiableFunctionLocalViewTag : public DifferentiableFunctionLocalViewBaseTag {};

  struct PowerDifferentiableFunctionTag : public DifferentiableFunctionBaseTag {};
  struct PowerDifferentiableFunctionLocalViewTag : public DifferentiableFunctionLocalViewBaseTag {};

  struct CompositeDifferentiableFunctionTag : public DifferentiableFunctionBaseTag {};
  struct CompositeDifferentiableFunctionLocalViewTag : public DifferentiableFunctionLocalViewBaseTag {};

  namespace Imp
  {

    template<typename Entity>
    struct PowerCompositeBindVisitor
      : public TypeTree::TreeVisitor, public TypeTree::DynamicTraversal
    {
      PowerCompositeBindVisitor(const Entity & e) : e_(e) {}
      template<typename LeafNode, typename TreePath>
      void leaf(LeafNode& node, TreePath treePath) const
      {
        node.bind(e_);
      }
      const Entity & e_;
    };

    struct PowerCompositeUnbindVisitor
      : public TypeTree::TreeVisitor, public TypeTree::DynamicTraversal
    {
      template<typename LeafNode, typename TreePath>
      void leaf(LeafNode& node, TreePath treePath) const
      {
        node.unbind();
      }
    };

    template<typename F>
    class LocalFunctionLeafNodeWrapper
      : public TypeTree::LeafNode
      , public F
    {
    public:
      typedef DifferentiableFunctionTag ImplementationTag;
      // LocalFunctionLeafNodeWrapper(F&& f) :
      //   F(std::forward(f))
      // {}
      LocalFunctionLeafNodeWrapper(const F& f) :
        F(f)
      {}
    };

  } // end namespace Imp

  template<class F, std::size_t k>
  class PowerLocalFunction
    : public TypeTree::PowerNode<F,k>
  {
    typedef TypeTree::PowerNode<F,k> NodeType;
  public:
    typedef DifferentiableFunctionTag ImplementationTag;

    #warning we should make this a free function
    //! Set the time in all leaf nodes of this function tree
    template <typename TT>
    void setTime(TT time){
      PowerCompositeSetTimeVisitor<TT> visitor(time);
      TypeTree::applyToTree(*this,visitor);
    }

    template <typename Entity>
    void bind(const Entity & e){
      Imp::PowerCompositeBindVisitor<Entity> visitor(e);
      TypeTree::applyToTree(*this,visitor);
    }

    void unbind(){
      Imp::PowerCompositeUnbindVisitor visitor;
      TypeTree::applyToTree(*this,visitor);
    }

    //! Deafult Constructor
    PowerLocalFunction()
    {}

    //! Construct a PowerGridFunction with k clones of the function t
    PowerLocalFunction (F& f)
      : NodeType(f) {}

    /** \brief Initialize all children with different function objects
     *
     *  @param t0 The initializer for the first child.
     *  @param t1 The initializer for the second child.
     *  @param ... more initializers
     */
    template<typename C0, typename C1, typename... Children>
    PowerLocalFunction (C0&& c0, C1&& c1, Children&&... children)
      : NodeType(std::forward(c0), std::forward(c1), std::forward(children)...)
    {
    }

    //! Transformation Constructor, taking the set of new children
    PowerLocalFunction(const Dune::array<Dune::shared_ptr<F>,k>& children)
      : NodeType(children)
    {}

  };

  template<typename... Children>
  class CompositeLocalFunction
    : public TypeTree::CompositeNode<Children...>
  {
    typedef TypeTree::CompositeNode<Children...> NodeType;
  public:
    typedef DifferentiableFunctionTag ImplementationTag;

    #warning we should make this a free function
    //! Set the time in all leaf nodes of this function tree
    template <typename TT>
    void setTime(TT time){
      PowerCompositeSetTimeVisitor<TT> visitor(time);
      TypeTree::applyToTree(*this,visitor);
    }

    template <typename Entity>
    void bind(const Entity & e){
      Imp::PowerCompositeBindVisitor<Entity> visitor(e);
      TypeTree::applyToTree(*this,visitor);
    }

    void unbind(){
      Imp::PowerCompositeUnbindVisitor visitor;
      TypeTree::applyToTree(*this,visitor);
    }

    //! Default Constructor
    CompositeLocalFunction()
    {}

    //! Initialize all children with the passed-in objects.
    template<typename... Args, typename = typename std::enable_if<(sizeof...(Args) == sizeof...(Children))>::type>
    CompositeLocalFunction(Args&&... args)
      : NodeType(std::forward<Args>(args)...)
    {}

  };

  template<class F, template<class> class DerivativeTraits = Functions::DefaultDerivativeTraits>
  class LocalGridViewFunctionAdapter
    : public TypeTree::LeafNode
  {
  public:
    using Range = typename F::Traits::RangeType;
    using LocalDomain = typename F::Traits::DomainType;
    using GridView = typename F::Traits::GridViewType;

    using Signature = Range(LocalDomain);
    using RawSignature =
      typename Functions::SignatureTraits<Signature>::RawSignature;
    using DerivativeSignature =
      typename DerivativeTraits<RawSignature>::Range(LocalDomain);

    using EntitySet = Functions::GridViewEntitySet<GridView, 0>;
    using Element = typename EntitySet::Element;
    using Geometry = typename std::decay<typename Element::Geometry>::type;

    typedef DifferentiableFunctionLocalViewTag ImplementationTag;

    // // Use the inderiction via derivativeIfImplemented to also support
    // // function types F that do not implement derivative. In this case
    // // the interface type DifferentiableFunction is used a dummy for
    // // the derivative type
    // using DerivativeDummy = DifferentiableFunction<DerivativeSignature>;
    // using GlobalRawDerivative = decltype(Imp::derivativeIfImplemented<DerivativeDummy, F>(std::declval<F>()));
    // using Derivative = AnalyticGridViewFunction<DerivativeSignature, GridView, GlobalRawDerivative, DerivativeTraits>;

    void bind(const Element& element)
    {
      element_ = element;
    }

    void unbind()
    {}

    Range operator()(const LocalDomain& x) const
    {
      Range v;
      f_->evaluate(element_, x, v);
      return v;
    }

    const Element& localContext() const
    {
      return element_;
    }

    // friend LocalDerivative derivative(const LocalAnalyticGridViewFunction& t)
    // {
    //     return LocalDerivative(Imp::derivativeIfImplemented<DerivativeDummy, F>(t.f_));
    // }

    LocalGridViewFunctionAdapter(const F & f) : f_(stackobject_to_shared_ptr(f)) {};

    // transforming constructor
    template<typename Transformation>
    LocalGridViewFunctionAdapter(shared_ptr<const F> f, const Transformation & t) : f_(f) {};

  private:
    Element element_;
    shared_ptr<const F> f_;
  };

  template<class F, template<class> class DerivativeTraits = Functions::DefaultDerivativeTraits>
  class GridViewFunctionAdapter
  {
  public:
    using Range = typename F::Traits::RangeFieldType;
    using Domain = typename F::Traits::DomainFieldType;
    using GridView = typename F::Traits::GridViewType;

    using Signature = Range(Domain);
    using RawSignature =
      typename Functions::SignatureTraits<Signature>::RawSignature;
    using DerivativeSignature =
      typename DerivativeTraits<RawSignature>::Range(Domain);

    using EntitySet = Functions::GridViewEntitySet<GridView, 0>;
    using Element = typename EntitySet::Element;
    using Geometry = typename Element::Geometry;

    // // Use the inderiction via derivativeIfImplemented to also support
    // // function types F that do not implement derivative. In this case
    // // the interface type DifferentiableFunction is used a dummy for
    // // the derivative type
    // using DerivativeDummy = DifferentiableFunction<DerivativeSignature>;
    // using GlobalRawDerivative = decltype(Imp::derivativeIfImplemented<DerivativeDummy, F>(std::declval<F>()));
    // using Derivative = AnalyticGridViewFunction<DerivativeSignature, GridView, GlobalRawDerivative, DerivativeTraits>;

    using LocalDomain = typename EntitySet::LocalCoordinate;
    using LocalFunction = LocalGridViewFunctionAdapter<F>; // , LocalDerivativeTraits<EntitySet, DerivativeTraits>::template Traits>;

    template<class FT>
    GridViewFunctionAdapter(FT&& f) :
      f_(std::forward<FT>(f))
    {}

    Range operator()(const Domain& x) const
    {
      Range v;
      f_.evaluate(x,v);
      return v;
    }

    // friend Derivative derivative(const AnalyticGridViewFunction& t)
    // {
    // }

    friend LocalFunction localFunction(const F& f)
    {
      return LocalFunction(f.f_);
    }

    const EntitySet& entitySet() const
    {
      return EntitySet(f_.getGridView());
    }

  private:
    F f_;
  };

  /*
    To create a tree with a local view we have to do destinguish
    between the following cases:

    a) HasFreeLocalFunction, but no NodeTag -> call localView(f) -> LocalViewLeafNodeWrapper
    b) operator(), but no Tag -> assert(false) // makeAnalyticGridViewFunction -> LocalViewLeafNodeWrapper
    c) HasFreeLocalFunction and NodeTag -> call localView(f)
    d) GridFunctionTag -> create a LocalGridViewFunctionAdapter
  */
  template<class F,
           typename std::enable_if<
             // case (a)
             Dune::Functions::Concept::models< Dune::Functions::Imp::HasFreeLocalFunction, F>()
             and
             not(TypeTree::has_node_tag<typename std::decay<F>::type>::value), int>::type = 0>
  auto makeLocalFunctionTree(F&& f)
    -> Imp::LocalFunctionLeafNodeWrapper< decltype(localView(f)) >
  {
    return Imp::LocalFunctionLeafNodeWrapper< decltype(localView(f)) >(localView(std::forward(f)));
  }

// template<class F>
// void makeLocalFunctionTree(F&& f)
// {
//     static_assert(AlwaysFalse<F>::value, "can't create a local function from this");
// }

  template<class F,
           typename std::enable_if<
             // case (c)
             Dune::Functions::Concept::models< Dune::Functions::Imp::HasFreeLocalFunction, F>()
             and
             TypeTree::has_node_tag<typename std::decay<F>::type>::value, int>::type = 0>
  auto makeLocalFunctionTree(F&& f)
    -> decltype(localView(f))
  {
    return localView(std::forward(f));
  }

  struct GridFunctionToLocalViewTransformation {};

  template<typename LeafNode>
  Dune::TypeTree::GenericLeafNodeTransformation<LeafNode,GridFunctionToLocalViewTransformation, LocalGridViewFunctionAdapter<LeafNode> >
  registerNodeTransformation(LeafNode* l, GridFunctionToLocalViewTransformation* t, GridFunctionTag* tag);

  template<typename PowerNode>
  Dune::TypeTree::SimplePowerNodeTransformation<PowerNode,GridFunctionToLocalViewTransformation,PowerLocalFunction>
  registerNodeTransformation(PowerNode* p, GridFunctionToLocalViewTransformation* t, PowerGridFunctionTag* tag);

  template<typename CompositeNode>
  Dune::TypeTree::SimpleCompositeNodeTransformation<CompositeNode,GridFunctionToLocalViewTransformation,CompositeLocalFunction>
  registerNodeTransformation(CompositeNode* c, GridFunctionToLocalViewTransformation* t, CompositeGridFunctionTag* tag);

  template<class F,
           typename std::enable_if<
             // case (d)
             IsGridFunction<F>::value, int>::type = 0>
  auto makeLocalFunctionTree(const F& f)
  -> typename Dune::TypeTree::TransformTree<typename std::decay<F>::type,
                                            GridFunctionToLocalViewTransformation>::transformed_type
  {
    // call the transformation
    return Dune::TypeTree::TransformTree<typename std::decay<F>::type,
                                         GridFunctionToLocalViewTransformation>::transform(f);
  }

} // end namespace PDELab
} // end namespace Dune

#endif // DUNE_PDELAB_FUNCTION_OLADINERFACEADAPTER_HH
