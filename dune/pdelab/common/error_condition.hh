#ifndef DUNE_PDELAB_COMMON_ERROR_CONDITION_HH
#define DUNE_PDELAB_COMMON_ERROR_CONDITION_HH

#include <dune/common/exceptions.hh>

#include <system_error>
#include <version>
#include <cstdlib>

namespace Dune::PDELab::inline Experimental {

/**
 * @brief Bad error condition exception
 * @details This class is used when an error condition is thrown. Similar to std::bad_expected_access
 */
class BadErrorCondition : public Dune::Exception {
public:
  //! Constructs a new bad error condition
  explicit BadErrorCondition(std::error_condition ec) : _ec{ec} {}

  //! Returns a reference to the stored value.
  [[nodiscard]] const std::error_condition& error() const& noexcept { return _ec; }

  //! Returns a reference to the stored value.
  [[nodiscard]] std::error_condition& error() & noexcept { return _ec; }

  //! Returns a reference to the stored value.
  [[nodiscard]] const std::error_condition&& error() && noexcept {return move(_ec); }
private:
  std::error_condition _ec;
};

/**
 * @brief Throwable extension of std::error_condition
 * @details This is a simple class to communicate dissapointments (i.e. error
 * conditions) as return values, but with an extension to throw errors
 * when they cannot be handled locally. This class is helpful in cases where
 * dissapointments are not exceptional and need to follow a deterministic path
 * in code (e.g. failure of iterative solver of an ill-possed problem).
 * @note When a value of this type is returned from a function, the user of the
 * function needs to handle the condition locally or defer it to another scope
 * with `or_throw`. Otherwise, the user will receive a compiler warning for
 * discarding the value. This ensures that dissapointments are inmediately
 * handled.
 *
 * @code{.cpp}
 * ErrorCondition foo() { ... }; // foo() communicates dissapoitments through its return type
 *
 * if (foo()) { // handle dissapoitment in local scope (deterministic path)
 *   // handle dissapointment
 * }
 *
 * try {
 *   foo().or_thorw(); // defer dissapoitment to parent scopes (non-deterministic path)
 * } catch (BadErrorCondition& ec) {
 *    // handle dissapointment
 * }
 *
 * foo(); // error condition is discarded, compiler warns that this should not be the case
 * @endcode
 */
struct [[nodiscard]] ErrorCondition : public std::error_condition {
  using std::error_condition::error_condition;

  /**
   * @brief Throw when error condition value is non-zero.
   * @details In case of an error condition, throws a BadErrorCondition, or
   * if C++ does not have exceptions enabled, program aborts.
   *
   * @throws BadErrorCondition
   */
  void or_throw() const {
    if (*this) [[unlikely]] {
#ifdef __cpp_exceptions
      BadErrorCondition exception{*this};
      exception.message(this->message());
      throw exception;
#else
      std::abort();
#endif
    }
  }
};

} // namespace Dune::PDELab::inline Experimental

#endif // DUNE_PDELAB_COMMON_ERROR_CONDITION_HH
