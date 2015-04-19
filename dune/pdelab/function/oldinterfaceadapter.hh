//-*- tab-width: 4; c-basic-offset: 2; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_FUNCTION_OLDINTERFACEADAPTER_HH
#define DUNE_PDELAB_FUNCTION_OLDINTERFACEADAPTER_HH

#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/common/typetraits.hh>
#include <dune/functions/common/signature.hh>
#include <dune/functions/common/defaultderivativetraits.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>

namespace Dune {
namespace PDELab {
namespace Imp {

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

    // transforming constructor
    template<typename Transformation>
    LocalGridViewFunctionAdapter(const F & f, const Transformation & t) : f_(stackobject_to_shared_ptr(f)) {};

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

} // end namespace Imp
} // end namespace PDELab
} // end namespace Dune

#endif // DUNE_PDELAB_FUNCTION_OLDINTERFACEADAPTER_HH
