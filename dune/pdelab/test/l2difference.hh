#ifndef DUNE_PDELAB_TEST_L2DIFFERENCE_HH
#define DUNE_PDELAB_TEST_L2DIFFERENCE_HH

#include <cmath>

#include <dune/geometry/type.hh>
#include <dune/geometry/quadraturerules.hh>

#include "../common/functionwrappers.hh"

#include "l2norm.hh"

// Calculate the squared L2 differerence of two functions
template<typename U, typename V> 
double l2difference2 (const U& u, const V &v, int qorder=1)
{
  return l2norm2
    (Dune::PDELab::makePointwiseGridFunctionAdapter
     (Dune::PDELab::PointwiseSumAdapterEngine(),
      u,
      Dune::PDELab::makePointwiseGridFunctionAdapter
      (Dune::PDELab::makePointwiseScaleAdapterEngine
      ((typename V::Traits::RangeFieldType)(-1)),
       v)),
     qorder);
}

// Calculate the L2 differerence of two functions
template<typename U, typename V> 
double l2difference (const U& u, const V &v, int qorder=1)
{
  return std::sqrt(l2difference2(u, v, qorder));
}

#endif // DUNE_PDELAB_TEST_L2DIFFERENCE_HH
