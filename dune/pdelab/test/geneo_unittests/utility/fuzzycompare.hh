#ifndef FUZZY_COMPARE_HH
#define FUZZY_COMPARE_HH

#include <cmath>
#include <algorithm>


/*
 * floating point comparisons based on the algorithms proposed in “The Art of
 * Computer Programming, Volume II: Seminumerical Algorithms (Addison-Wesley,
 * 1969)”
 */


namespace Utility
{

  // Comparison with 1.0 for floating point types
  template<typename T>
  bool is_numeric_one(const T& x, double eps=1e-8)
  {
    return (std::abs(x - 1.0) <= std::max(std::abs(x), 1.0) * eps);
  }


  // comparison of two floating point numbers with absolute precision aEps
  // and relative precision rEps
  template<typename T>
  bool is_equal(const T& x1, const T& x2, double aEps=1e-12, double rEps=1e-8)
  {
    double diff{ std::abs(x1 - x2) };
    if (diff <= aEps)
      return true;

    return (diff <= (std::max(std::abs(x1), std::abs(x2)) * rEps));
  }

} // namespace Utility

#endif
