#ifndef DUNE_PDELAB_OPERATOR_ISTL_ADAPTER_HH
#define DUNE_PDELAB_OPERATOR_ISTL_ADAPTER_HH


#include <dune/pdelab/operator/operator.hh>

#include <dune/pdelab/common/algebra.hh>
#include <dune/pdelab/common/container_entry.hh>

#include <dune/istl/operators.hh>
#include <dune/istl/scalarproducts.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/superlu.hh>
#include <dune/istl/umfpack.hh>
#include <dune/istl/preconditioners.hh>

#include <functional>
#include <concepts>

namespace Dune::PDELab::ISTL::inline Experimenal {


template<class Domain, std::copy_constructible Range>
class LinearSolver : public Operator<Range,Domain>
{

    // converts a dune-assembler forward operator into a dune-istl linear operator
    class ISTLLinearAdapter : public Dune::LinearOperator<Domain,Range> {
    public:
      ISTLLinearAdapter(Operator<Domain,Range>& forward_op)
        : _forward{ forward_op }
      {}

      void apply (const Domain& x, Range& y) const override {
        forEachContainerEntry(std::execution::par_unseq, y, []<class T>(T& v){v = T{0};});
        _forward.apply(x,y).or_throw();
      }

      using typename Dune::LinearOperator<Domain,Range>::field_type;

      //! apply operator to x, scale and add:  \f$ y = y + \alpha A(x) \f$
      void applyscaleadd (field_type alpha, const Domain& x, Range& y) const override {
        tmp = y;
        this->apply(x,tmp);
        Dune::PDELab::axpy(std::execution::par_unseq, y, alpha, tmp);
      }

      SolverCategory::Category category() const override {
        return SolverCategory::Category::sequential;
      }

    private:
      Operator<Domain,Range>& _forward;
      mutable Range tmp;
    };
public:

  LinearSolver() {}
  virtual ErrorCondition apply(Range& range, Domain& domain) override {
    uint64_t solver_timestamp = perfetto::TrackEvent::GetTraceTimeNs();
    TRACE_EVENT("dune", "LinearSolver", solver_timestamp);
    static_assert(std::is_same_v<Range,Domain>);
    auto& forward = this->template get<PDELab::Operator<Domain, Range>>("forward");

    auto istl_forward_op = std::make_shared<ISTLLinearAdapter>(forward);

    auto scalar_product_op = std::make_shared<Dune::ScalarProduct<Range>>();
    std::shared_ptr<Preconditioner<Domain,Range>> pre_op;
    pre_op = std::make_shared<Dune::Richardson<Domain,Range>>(0.1);

    // // try to get the jacobian...
    // using JacobianOperator = PDELab::AssembledLinearOperator<Jacobian,Domain,Range>;
    // auto assembled_derivative = dynamic_cast<JacobianOperator const *>(&derivative);
    // if (assembled_derivative) {
    //   // pre_op = std::make_shared<Dune::SeqSSOR<Jacobian,Domain,Range,2>>(assembled_derivative->matrix(), 5, 1);
    // }

    int verbosity = 4;
    int max_it = 20000;
    auto rel_tol = this->template get<double>("convergence_condition.relative_tolerance");
    Dune::InverseOperatorResult result;
    Dune::PDELab::ErrorCondition ec{};
    // // if (solver == "BiCGSTAB") {
    //   auto istl_solver = Dune::BiCGSTABSolver<Domain>{istl_forward_op, scalar_product_op, pre_op, rel_tol, int(max_it), verbosity};
    //   istl_solver.apply(x, b, result);
    //   DUNE_THROW(NotImplemented, "");
    // } else if (istl_solver == "RestartedGMRes") {
      int restart = 40;
      // auto istl_solver = Dune::CGSolver<Domain>{istl_forward_op, scalar_product_op, pre_op, rel_tol, max_it, verbosity};
      auto istl_solver = Dune::RestartedGMResSolver<Domain, Range>{istl_forward_op, scalar_product_op, pre_op, rel_tol, restart, int(max_it), verbosity};

      istl_solver.apply(domain, range, result);
    // } else {
    //   DUNE_THROW(IOError, "Not known linear solver");
    // }

    TRACE_COUNTER("dune", "LinearSolver::Iterations", solver_timestamp, result.iterations);
    TRACE_COUNTER("dune", "LinearSolver::Reduction", solver_timestamp, result.reduction);
    TRACE_COUNTER("dune", "LinearSolver::Converged", solver_timestamp, result.converged);

    if (result.converged)
      return ec;
    else
      return make_error_condition(Dune::PDELab::Convergence::Reason::DivergedNull);
  }

  virtual ErrorCondition apply(const Range& range, Domain& domain) override {
    Domain tmp_range = range;
    return this->apply(tmp_range, domain);
  }
};

} // namespace Dune::PDELab::inline Experimental::ISTL

#endif // DUNE_PDELAB_OPERATOR_ISTL_ADAPTER_HH
