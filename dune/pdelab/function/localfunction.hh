//-*- tab-width: 4; c-basic-offset: 2; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_FUNCTION_LOCALFUNCTION_HH
#define DUNE_PDELAB_FUNCTION_LOCALFUNCTION_HH

#include <type_traits>

#include <dune/common/array.hh>
#include <dune/pdelab/function/tags.hh>
#include <dune/pdelab/function/localfunctionhelper.hh>
#include <dune/pdelab/function/oldinterfaceadapter.hh>
#include <dune/functions/common/signature.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>

namespace Dune {
namespace PDELab {

  template<class F, std::size_t k>
  class PowerLocalFunction
    : public TypeTree::PowerNode<F,k>
  {
    typedef TypeTree::PowerNode<F,k> NodeType;
  public:
    typedef PowerDifferentiableFunctionLocalViewTag ImplementationTag;

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
    typedef CompositeDifferentiableFunctionLocalViewTag ImplementationTag;

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

  /** \brief turn a given function tree into a local function tree

    To create a tree with a local view we have to do destinguish
    between the following cases:

    a) HasFreeLocalFunction, but no NodeTag -> call localView(f) -> LocalViewLeafNodeWrapper
    b) operator(), but no Tag -> makeGridViewFunction -> LocalViewLeafNodeWrapper
    c) HasFreeLocalFunction and NodeTag -> call localView(f)
    d) GridFunctionTag -> create a LocalGridViewFunctionAdapter
    e) LeafNodeTag, derived from GridFunctionInterface -> create a LocalGridViewFunctionAdapter (this is always a leaf node)
  */
  template<class F, class GV,
           typename std::enable_if<
             // case (a)
             models< Dune::Functions::Imp::HasFreeLocalFunction, F>()
             and
             not(TypeTree::has_node_tag<typename std::decay<F>::type>::value), int>::type = 0>
  auto makeLocalFunctionTree(const F& f, const GV & gv)
    -> Imp::LocalFunctionLeafNodeWrapper< decltype(localFunction(f)) >
  {
    return Imp::LocalFunctionLeafNodeWrapper< decltype(localFunction(f)) >(localFunction(f));
  }

  template<class F, class GV,
           // case (b)
           typename std::enable_if<
             Dune::Functions::Concept::isCallable<F, typename GV::template Codim<0>::Entity::Geometry::GlobalCoordinate>()
             and
             not(models< Dune::Functions::Imp::HasFreeLocalFunction, F>())
             and
             not(TypeTree::has_node_tag<typename std::decay<F>::type>::value), int>::type = 0>
  auto makeLocalFunctionTree(const F& f, const GV & gv)
    -> decltype(makeLocalFunctionTree(Functions::makeGridViewFunction(f,gv), gv))
  {
    return makeLocalFunctionTree(Functions::makeGridViewFunction(f,gv), gv);
  }

  template<class F, class GV,
           typename std::enable_if<
             // case (c)
             models< Dune::Functions::Imp::HasFreeLocalFunction, F>()
             and
             TypeTree::has_node_tag<typename std::decay<F>::type>::value, int>::type = 0>
  auto makeLocalFunctionTree(F&& f, const GV & gv)
    -> decltype(localView(f))
  {
    return localView(std::forward(f));
  }

  struct GridFunctionToLocalViewTransformation {};

  template<typename LeafNode>
  Dune::TypeTree::GenericLeafNodeTransformation<LeafNode,GridFunctionToLocalViewTransformation, Imp::LocalGridViewFunctionAdapter<LeafNode> >
  registerNodeTransformation(LeafNode* l, GridFunctionToLocalViewTransformation* t, GridFunctionTag* tag);

  template<typename PowerNode>
  Dune::TypeTree::SimplePowerNodeTransformation<PowerNode,GridFunctionToLocalViewTransformation,PowerLocalFunction>
  registerNodeTransformation(PowerNode* p, GridFunctionToLocalViewTransformation* t, PowerGridFunctionTag* tag);

  template<typename CompositeNode>
  Dune::TypeTree::SimpleCompositeNodeTransformation<CompositeNode,GridFunctionToLocalViewTransformation,CompositeLocalFunction>
  registerNodeTransformation(CompositeNode* c, GridFunctionToLocalViewTransformation* t, CompositeGridFunctionTag* tag);

  template<class F, class GV,
           typename std::enable_if<
             // case (d)
             IsGridFunction<F>::value, int>::type = 0>
  auto makeLocalFunctionTree(const F& f, const GV & gv)
  -> typename Dune::TypeTree::TransformTree<typename std::decay<F>::type,
                                            GridFunctionToLocalViewTransformation>::transformed_type
  {
    // call the transformation
    return Dune::TypeTree::TransformTree<typename std::decay<F>::type,
                                         GridFunctionToLocalViewTransformation>::transform(f);
  }

  template<class F, class GV,
           typename std::enable_if<
             // case (e)
             not(IsGridFunction<F>::value)
             &&
             std::is_same<TypeTree::NodeTag<F>,TypeTree::LeafNodeTag>::value, int>::type = 0>
  auto makeLocalFunctionTree(const GridFunctionInterface<typename F::Traits,F>& f, const GV & gv)
    -> Imp::LocalGridViewFunctionAdapter<F>
  {
    // call the transformation
    return Imp::LocalGridViewFunctionAdapter<F>(static_cast<const F&>(f));
  }

} // end namespace PDELab
} // end namespace Dune

#endif // DUNE_PDELAB_FUNCTION_LOCALFUNCTION_HH
