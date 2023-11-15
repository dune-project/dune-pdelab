#ifndef DUNE_PDELAB_OPERATOR_LINE_SEARCH_HH
#define DUNE_PDELAB_OPERATOR_LINE_SEARCH_HH

#include <dune/pdelab/operator/operator.hh>

#include <dune/pdelab/common/error_condition.hh>
#include <dune/pdelab/common/algebra.hh>
#include <dune/pdelab/common/trace.hh>
#include <dune/pdelab/common/property_tree.hh>
#include <dune/pdelab/common/execution.hh>

#include <memory>
#include <functional>

namespace Dune::PDELab::inline Experimental {

  template<class Domain, class Range, class RangeField>
  class LineSearch : public PropertyTree {
  public:
    using Norm = std::function<RangeField(const Range&)>;

    virtual ~LineSearch () {}

    virtual ErrorCondition apply(
      Operator<Domain,Range>& op,
      const Norm& norm_op,
      Domain& point,
      const Domain& direction,
      Range& residual,
      RangeField& defect
    ) const = 0;
  };

  //! Class for simply updating the solution without line search
  template <class Domain, class Range, class RangeField>
  class LineSearchNone : public LineSearch<Domain,Range,RangeField>
  {
  public:
    using Norm = std::function<RangeField(const Range&)>;

    //! Do line search (in this case just update the solution)
    virtual ErrorCondition apply(
      Operator<Domain,Range>& op,
      const Norm& norm_op,
      Domain& point,
      const Domain& direction,
      Range& residual,
      RangeField& defect) const override
    {
      TRACE_EVENT("dune", "Newton::LineSearchNone");

      // point -= direction
      PDELab::axpy(PDELab::Execution::par_unseq, point, -1.0, direction);

      // evaluate || f(x) + z ||
      op.apply(point, residual).or_throw();
      defect = norm_op(residual);
      return {};
    }
  };

} // namespace Dune::PDELab::inline Experimental

#endif // DUNE_PDELAB_OPERATOR_LINE_SEARCH_HH
