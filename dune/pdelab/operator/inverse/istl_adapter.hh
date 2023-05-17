#ifndef DUNE_PDELAB_OPERATOR_ISTL_ADAPTER_HH
#define DUNE_PDELAB_OPERATOR_ISTL_ADAPTER_HH

#include <dune/istl/operators.hh>

#include <dune/pdelab/operator/operator.hh>

#include <dune/pdelab/common/algebra.hh>
#include <dune/pdelab/common/container_entry.hh>

#include <functional>

namespace Dune::PDELab::inline Experimental::ISTL {

// converts a dune-assembler forward operator into a dune-istl linear operator
template<class Domain, class Range>
class LinearAdapter : public Dune::LinearOperator<Domain,Range> {
public:
  LinearAdapter(std::shared_ptr<Operator<Domain,Range>> forward_op)
    : _forward{ std::move(forward_op) }
  {}

  virtual void apply (const Domain& x, Range& y) const override {
    forEachContainerEntry(std::execution::par_unseq, y, []<class T>(T& v){v = T{0};});
    _forward->apply(x,y).or_throw();
  }

  using typename Dune::LinearOperator<Domain,Range>::field_type;

  //! apply operator to x, scale and add:  \f$ y = y + \alpha A(x) \f$
  virtual void applyscaleadd (field_type alpha, const Domain& x, Range& y) const override {
    tmp = y;
    this->apply(x,tmp);
    Dune::PDELab::axpy(std::execution::par_unseq, y, alpha, tmp);
  }

  virtual SolverCategory::Category category() const override {
    return SolverCategory::Category::sequential;
  }

private:
  std::shared_ptr<Operator<Domain,Range>> _forward;
  mutable Range tmp;
};


} // namespace Dune::PDELab::inline Experimental::ISTL

#endif // DUNE_PDELAB_OPERATOR_ISTL_ADAPTER_HH
