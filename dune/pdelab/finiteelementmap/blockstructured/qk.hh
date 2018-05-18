// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_BLOCKSTRUCTURED_QK_LOCALFINITEELEMENT_HH
#define DUNE_BLOCKSTRUCTURED_QK_LOCALFINITEELEMENT_HH

#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include "qklocalinterpolation.hh"
#include "qklocalbasis.hh"
#include "qklocalcoefficients.hh"

namespace Dune
{
  namespace Blockstructured {
    /** \brief General Lagrange finite element for cubes with arbitrary dimension and polynomial order
     *   \note The general class QkLocalCoefficients is available for k>0 in dimensions 2 and 3 only
     *
     * \tparam D type used for domain coordinates
     * \tparam R type used for function values
     * \tparam d dimension of the reference element
     * \tparam k polynomial order
     */
    template<class D, class R, std::size_t d, std::size_t k, std::size_t blocks>
    class QkLocalFiniteElement {

      typedef QkLocalBasis<D, R, k, d, blocks> LocalBasis;
      typedef QkLocalCoefficients<k, d, blocks> LocalCoefficients;
      typedef QkLocalInterpolation<k, d, blocks, LocalBasis> LocalInterpolation;

    public:

      /** \todo Please doc me !
       */
      typedef Dune::LocalFiniteElementTraits <LocalBasis, LocalCoefficients, LocalInterpolation> Traits;

      /** \todo Please doc me !
       */
      QkLocalFiniteElement() {}

      /** \todo Please doc me !
       */
      const typename Traits::LocalBasisType &localBasis() const {
        return basis;
      }

      /** \todo Please doc me !
       */
      const typename Traits::LocalCoefficientsType &localCoefficients() const {
        return coefficients;
      }

      /** \todo Please doc me !
       */
      const typename Traits::LocalInterpolationType &localInterpolation() const {
        return interpolation;
      }

      /** \brief Number of shape functions in this finite element */
      unsigned int size() const {
        return basis.size();
      }

      /** \todo Please doc me !
       */
      static constexpr GeometryType type() {
        return GeometryTypes::cube(d);
      }

    private:
      LocalBasis basis;
      LocalCoefficients coefficients;
      LocalInterpolation interpolation;
    };
  }

}

#endif
