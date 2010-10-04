// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_FINITEELEMENT_TRAITS_HH
#define DUNE_PDELAB_FINITEELEMENT_TRAITS_HH

#include <cstddef>

#include <dune/common/static_assert.hh>
#include <dune/common/typetraits.hh>

namespace Dune {
  namespace PDELab {

    template<class FiniteElement, class = void>
    struct FiniteElementTraits {
      typedef typename FiniteElement::Traits::Basis Basis;
      typedef typename FiniteElement::Traits::Interpolation Interpolation;
      typedef typename FiniteElement::Traits::Coefficients Coefficients;

      static const Basis &basis(const FiniteElement& fe)
      { return fe.basis(); }
      static const Interpolation &interpolation(const FiniteElement& fe)
      { return fe.interpolation(); }
      static const Coefficients &coefficients(const FiniteElement& fe)
      { return fe.coefficients(); }
    };

    template<class FiniteElement>
    struct FiniteElementTraits
      < FiniteElement,
        typename enable_if<AlwaysTrue<typename FiniteElement::Traits::
                                      LocalBasisType>::value>::type>
    {
      typedef typename FiniteElement::Traits::LocalBasisType Basis;
      typedef typename FiniteElement::Traits::LocalInterpolationType
        Interpolation;
      typedef typename FiniteElement::Traits::LocalCoefficientsType
        Coefficients;

      static const Basis &basis(const FiniteElement& fe)
      { return fe.localBasis(); }
      static const Interpolation &interpolation(const FiniteElement& fe)
      { return fe.localInterpolation(); }
      static const Coefficients &coefficients(const FiniteElement& fe)
      { return fe.localCoefficients(); }
    };

    template<class Basis, class = void>
    struct BasisTraits :
      public Basis::Traits
    { };

    template<class Basis>
    struct BasisTraits<Basis,
                       typename enable_if<Basis::Traits::dimDomain>::type>
    {
      typedef typename Basis::Traits::DomainFieldType DomainField;
      static const std::size_t dimDomainLocal = Basis::Traits::dimDomain;
      static const std::size_t dimDomainGlobal = Basis::Traits::dimDomain;
      typedef typename Basis::Traits::DomainType DomainLocal;
      typedef typename Basis::Traits::DomainType DomainGlobal;

      typedef typename Basis::Traits::RangeFieldType RangeField;
      static const std::size_t dimRange = Basis::Traits::dimRange;
      typedef typename Basis::Traits::RangeType Range;

      typedef typename Basis::Traits::JacobianType Jacobian;

      static const std::size_t diffOrder = Basis::Traits::diffOrder;
    };

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_FINITEELEMENT_TRAITS_HH
