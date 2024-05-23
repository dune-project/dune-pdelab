#ifndef DUNE_PDELAB_BASIS_DISCRETE_GLOBAL_FUNCTION_HH
#define DUNE_PDELAB_BASIS_DISCRETE_GLOBAL_FUNCTION_HH

#include <dune/pdelab/concepts/basis.hh>

#include <dune/pdelab/common/local_container.hh>
#include <dune/pdelab/common/container_entry.hh>

#include <dune/functions/gridfunctions/gridviewentityset.hh>

#include <dune/grid/utility/hierarchicsearch.hh>

#include <memory>

namespace Dune::PDELab::inline Experimenal {


// models a dune-functions grid function
template<PDELab::Concept::Basis Basis, class Container>
class ScalarDiscreteFunction {

public:
  static_assert(PDELab::Concept::LeafTreeNode<typename Basis::LocalView::Tree>);

  using EntitySet = Functions::GridViewEntitySet<typename Basis::EntitySet, 0>;

private:
  using Range = std::remove_cvref_t<decltype(containerEntry(std::declval<Container>(), std::declval<typename Basis::LocalView::Tree::MultiIndex>()))>;
  using Domain = typename EntitySet::GlobalCoordinate;

  class LocalFunction
  {
    using Domain = typename EntitySet::LocalCoordinate;
  public:

    using LocalContext = typename EntitySet::Element;

    LocalFunction(const Basis& basis, const std::shared_ptr<const Container>& container)
      : _container{container}
      , _lbasis{basis.localView()}
      , _lcontainer{basis, _container.get()}
    {}

    LocalFunction(const LocalFunction& other) = default;
    LocalFunction(LocalFunction&& other) = default;

    LocalFunction& operator=(const LocalFunction& other) = default;
    LocalFunction& operator=(LocalFunction&& other) = default;

    void bind(const LocalContext& element)
    {
      _lbasis->bind(element);
      _lcontainer.load(*_lbasis);
      _bound = true;
    }

    Range operator()(const Domain& xlocal, std::optional<Range> storage = {}) const
    {
      const auto& node = _lbasis->tree();
      if (node.size() == 0) DUNE_THROW(RangeError, "LocalFunction has no support in requested domain");
      node.finiteElement().localBasis().evaluateFunction(xlocal, _fe_range);
      Range value = std::move(storage).value_or(0.);
      for (std::size_t dof = 0; dof < node.size(); ++dof)
        value += _lcontainer(node, dof) * _fe_range[dof];
      return value;
    }

    void unbind()
    {
      _lcontainer.clear(*_lbasis);
      _lbasis->unbind();
      _bound = true;
    }

    bool bound() const {
      return _bound;
    }

    const LocalContext& localContext() const {
      return _lbasis->element();
    }

  private:
    using FiniteElement = typename Basis::LocalView::Tree::FiniteElement;
    using FERange = typename FiniteElement::Traits::LocalBasisType::Traits::RangeType;

    std::shared_ptr<const Container> _container;
    std::optional<typename Basis::LocalView> _lbasis;
    PDELab::LocalContainerBuffer<Basis,const Container> _lcontainer;
    mutable std::vector<FERange> _fe_range;
    bool _bound = false;
  };

public:
  friend LocalFunction localFunction(const ScalarDiscreteFunction& discrete_function) {
    return LocalFunction{discrete_function._basis, discrete_function._container};
  }


  ScalarDiscreteFunction(const Basis& basis, const std::shared_ptr<const Container>& container)
    : _basis{basis}
    , _container{container}
    , _entity_set{_basis.entitySet()}
  {}

  // very slow!
  Range operator()(const Domain& xglobal, std::optional<Range> storage = {}) const {
    HierarchicSearch search{_basis.entitySet().grid(), _basis.entitySet().indexSet()};
    const auto e = search.findEntity(xglobal);
    auto lfunc = localFunction(*this);
    lfunc.bind(e);
    return lfunc(e.geometry().local(xglobal), std::move(storage));
  }

  const EntitySet& entitySet() const {
    return _entity_set;
  }

private:
  Basis _basis;
  std::shared_ptr<const Container> _container;
  EntitySet _entity_set;
};

template <Concept::Basis Basis, class Container>
requires PDELab::Concept::LeafTreeNode<typename Basis::LocalView::Tree>
auto makeDiscreteGlobalBasisFunction(const Basis &basis, std::shared_ptr<Container> container)
{
  return ScalarDiscreteFunction<Basis, Container>{ basis,
                                                   std::move(container) };
}

} // namespace Dune::PDELab::inline Experimental

#endif // DUNE_PDELAB_BASIS_DISCRETE_GLOBAL_FUNCTION_HH
