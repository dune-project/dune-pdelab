// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_COMMON_RANGE_HH
#define DUNE_PDELAB_COMMON_RANGE_HH

#include<vector>

namespace Dune {
  namespace PDELab {

    //! Generate a vector with a range of numbers
    /**
     * This is mostly useful to supply default function arguments.
     *
     * The generated vector will consist of the numbers begin,
     * begin+increment etc. through (but not including) passed_the_end.  For
     * instance, rangeVector(23, 42, 7) will yield the vector (23, 30, 37).
     * The increment must be strictly positive at the moment.
     */
    template<class T>
    std::vector<T> rangeVector(T begin, T passed_the_end, T increment = 1) {
      std::vector<T> tmp;
      tmp.reserve((passed_the_end-begin)/increment);
      for(T i = begin; i < passed_the_end; i+=increment)
        tmp.push_back(i);
      return tmp;
    }

    //! Generate a vector with a range of numbers
    /**
     * Equivalent to rangeVector(0, passed_the_end).
     */
    template<class T>
    std::vector<T> rangeVector(T passed_the_end)
    { return rangeVector(T(0), passed_the_end); }

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_COMMON_RANGE_HH
