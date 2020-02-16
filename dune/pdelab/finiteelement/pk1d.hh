// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// Pk in one dimension with k as runtime variable

#ifndef DUNE_PDELAB_FINITEELEMENT_PK1D_HH
#define DUNE_PDELAB_FINITEELEMENT_PK1D_HH

#warning This file is deprecated. Please use LagrangeLocalFiniteElement from dune-localfunctions directly!

#include <dune/geometry/type.hh>

#include <dune/localfunctions/lagrange.hh>
#include <dune/localfunctions/lagrange/equidistantpoints.hh>

namespace Dune {

  /** \brief Define the Pk Lagrange basis functions in 1d on the reference interval
   *
   * Unlike the corresponding implementation in dune-localfunctions, the order k
   * is a run-time parameter here!
   *
   *  \tparam D Type to represent domain.
   *  \tparam R Type to represent range.
   */
  template<class D, class R>
  class Pk1dLocalFiniteElement
  : public LagrangeLocalFiniteElement<EquidistantPointSet,1,D,R>
  {

    using Base = LagrangeLocalFiniteElement<EquidistantPointSet,1,D,R>;
  public:

    Pk1dLocalFiniteElement (std::size_t k)
      : Base(GeometryTypes::cube(1), k)
    {}

  };
}

#endif // DUNE_PDELAB_FINITEELEMENT_PK1D_HH
