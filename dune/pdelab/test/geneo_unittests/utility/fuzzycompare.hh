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


  // comparison of two istl block matrices of block size 1
  template<class Matrix>
  bool matrices_equal(const Matrix& m1, const Matrix& m2, double aEps=1e-12,
                      double rEps=1e-8)
  {
    const std::size_t rows1{ m1.N() };
    const std::size_t cols1{ m1.M() };

    const std::size_t rows2{ m2.N() };
    const std::size_t cols2{ m2.M() };

    // check if matrices have the same size
    if (rows1 != rows2 || cols1 != cols2)
      return false;

    bool same(false);

    for (int i{ 0 }; i < rows1; ++i)
      for (int j{ 0 }; j < cols1; ++j)
        if (m1.exists(i, j))
        {
          // check if matrices have the same pattern
          if (!m2.exists(i, j))
            return false;

          // check if matrices have the same entries
          same = Utility::is_equal(m1[i][j], m2[i][j], aEps, rEps);

          if (!same)
            return false;
        }

    return true;
  }

} // namespace Utility

#endif
