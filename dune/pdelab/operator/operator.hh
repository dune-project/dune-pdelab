#ifndef DUNE_PDELAB_OPERATOR_OPERATOR_HH
#define DUNE_PDELAB_OPERATOR_OPERATOR_HH

#include <dune/pdelab/common/error_condition.hh>
#include <dune/pdelab/common/property_tree.hh>

#include <utility>
#include <map>
#include <any>
#include <memory>
#include <string>
#include <cassert>

namespace Dune::PDELab::inline Experimental {

/**
 * @brief Operator @f$f@f$ from a domain @f$X@f$ to a range @f$X@f$ (i.e. @f$f:X\to Y@f$)
 * @details This object providers two versions of the operator, one where the
 * domain is constant and one where the domain is mutable.
 * @note The operator application is accumulative (i.e. @f$y_{out} := y_{in} + f(x)@f$).
 * This implies that range objects @f$y_{in}@f$ shall always be initialized and,
 * in case of vectors, be correctly sized.
 *
 * @tparam DomainType  Domain of @f$f@f$
 * @tparam RangeType   Range of @f$f@f$
 */
template<class DomainType, class RangeType>
class Operator : public PropertyTree
{
public:

  Operator() = default;
  Operator(const Operator&) = delete;
  Operator(Operator&&) = delete;

  virtual ~Operator() {}

  //! The domain of @f$f@f$
  using Domain = DomainType;

  //! The range of @f$f@f$
  using Range = RangeType;

  /**
   * @brief Application of @f$y_{out} := y_{in} + f(x)@f$
   *
   * @param domain   Value @f$x@f$ to evaluate @f$f@f$
   * @param range    Input/Output @f$y@f$ of the operator evaluation
   * @return         Value expresssing success/failure of the operation
   */
  virtual ErrorCondition apply(const Domain& domain, Range& range) = 0;

  /**
   * @brief Application of @f$y_{out} := y_{in} + f(x)@f$
   * @details This version of the operator should be prefered when the domain @f$x@f$
   * is not needed after evaluating the operator. This is because some implementations
   * may choose to perfom some computations in-place to be faster.
   * @warning The domain object may be mutated! Avoid this by passing a constant domain.
   *
   * @param domain   Mutable value @f$x@f$ to evaluate @f$f@f$
   * @param range    Input/Output  @f$y@f$ of the operator evaluation
   * @return         Value expresssing success/failure of the operation
   */
  virtual ErrorCondition apply(Domain& domain, Range& range) {
    return this->apply(std::as_const(domain), range);
  }

  /**
   * @brief Obtain a pointer to a derivative @f$Df(a)@f$ of @f$f@f$ at a point @f$a\in X@f$
   * @details The resulting operator is intended to be used where the derivatives
   * of this object are needed (e.g., newton's method), thus, operators are allowed
   * to return an approximate version of the derivative if the original one is not
   * possible or is too expensive to compute.
   *
   * @param x            Position of the derivative
   * @param reuse_dx     If not empty, the contests may be reused to build the new derivative
   * @return             operator holding the derivative
   */
  virtual std::shared_ptr<Operator<Domain,Range>> derivative(const Domain& x, std::shared_ptr<Operator<Domain,Range>> reuse_dx = nullptr) const {
    DUNE_THROW(NotImplemented, "Not known derivative");
    return nullptr;
  }
};


/**
 * @brief Forward time-evolution operator @f$f_{\delta_t}:X\to X@f$ advancing one time-step @f$\delta_t@f$
 * @see <a href="https://en.wikipedia.org/wiki/Time_evolution">Time evolution</a>
 *
 * @tparam State      Domain/range of @f$f_{\delta_t}@f$
 * @tparam Duration   Representation of @f$\delta_t@f$
 */
template<class State>
using OneStep = Operator<State, State>;


/**
 * @brief Operator @f$f^{-1}:Y\to X@f$ of a Forward operator @f$f:X\to Y@f$
 *
 * @note We solve the accumulative form of the Forward operator (i.e. @f$\tilde{f}(x) := f(x) + y = 0@f$),
 * thus, the domain of @f$f^{-1}@f$ needs to be passed scaled with -1.
 *
 * @tparam DomainType  Domain of @f$f@f$
 * @tparam RangeType   Range of @f$f@f$
 */
template<class Domain, class Range>
using Inverse = Operator<Range, Domain>;

} // namespace Dune::PDELab::inline Experimental

#endif // DUNE_PDELAB_OPERATOR_OPERATOR_HH
