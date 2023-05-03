#ifndef DUNE_PDELAB_COMMON_CONVERGENCE_REASON_HH
#define DUNE_PDELAB_COMMON_CONVERGENCE_REASON_HH

#include <dune/pdelab/common/error_condition.hh>

#include <string>

namespace Dune::PDELab::inline Experimental::Convergence {

/**
 * @brief Enumerator describing the reason of convergence
 * @details Negative numbers encode a non-converged reason and positive number a success.
 * These numbers roughly follow the PETSc convergence codes.
 */
enum class Reason {
  /* converged */
  ConvergedByRelativeTolerance        =  2,
  ConvergedByAbsoluteTolerance        =  3,
  ConvergedByIterations               =  4,
  ConvergedByStepLength               =  7,
  ConvergedByHappyBreakdown           =  8,
  /* diverged */
  DivergedNull                        = -2,
  DivergedByIterations                = -3,
  DivergedByDivergenceTolarance       = -4,
  DivergedByBreakdown                 = -5,
  DivergedByNonSymmetry               = -7,
  DivergedByIndefinitePreconditioner  = -8,
  DivergedByNanOrInf                  = -9,
  DivergedByIndefiniteMatrix          = -10,
  DivergedByFailedPreconditioner      = -11,

  Iterating                           = -1
};

namespace {

//! Category for error conditions arasing from iterative operators
class IterativeOperatorErrorCategory : public std::error_category {
public:

  //! Obtains the name of the category
  [[nodiscard]] const char* name() const noexcept override {
    return "iterative operator";
  }

  //! Obtains the explanatory string
  [[nodiscard]] virtual std::string message(int ec) const override {
    switch (static_cast<Convergence::Reason>(ec)) {
      case Convergence::Reason::DivergedNull:
        return "diverged by unkonw reason";
      case Convergence::Reason::DivergedByIterations:
        return "diverged by the maximum number of iterations";
      case Convergence::Reason::DivergedByDivergenceTolarance:
        return "diverged by a divergence tolerace";
      case Convergence::Reason::DivergedByBreakdown:
        return "diverged by solver breakdown";
      case Convergence::Reason::DivergedByNonSymmetry:
        return "diverged by detected non-symmetry";
      case Convergence::Reason::DivergedByIndefinitePreconditioner:
        return "diverged by detected indefinite preconditioner";
      case Convergence::Reason::DivergedByNanOrInf:
        return "diverged by resulting NaN or Inf values";
      case Convergence::Reason::DivergedByIndefiniteMatrix:
        return "diverged by detected indefinite matrix";
      case Convergence::Reason::DivergedByFailedPreconditioner:
        return "diverged by failure in the preconditioner";
      case Convergence::Reason::Iterating:
        return "still iterating";
      default:
        return "unkwnown reason";
    }
  }
};

//! Category for error conditions arasing from iterative operators
const IterativeOperatorErrorCategory iterative_operator_error_category {};
}

//! Creates an error condition for a convergence reason value cr.
ErrorCondition make_error_condition(Convergence::Reason cr) {
  // only negative values are errors
  int ecr = std::min(static_cast<int>(cr), int{0});
  return {ecr, iterative_operator_error_category};
}

} // namespace Dune::PDELab::inline Experimental::Convergence

namespace std {

//! Enable automatic conversions on Dune::PDELab::Convergence::Reason to error condition
template <>
struct is_error_code_enum<Dune::PDELab::Convergence::Reason> : public true_type {};

} // namespace std

#endif // DUNE_PDELAB_COMMON_CONVERGENCE_REASON_HH
