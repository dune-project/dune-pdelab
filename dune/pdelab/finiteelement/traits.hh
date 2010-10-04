// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_FINITEELEMENT_TRAITS_HH
#define DUNE_PDELAB_FINITEELEMENT_TRAITS_HH

#include <cstddef>
#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/shared_ptr.hh>
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

      typedef shared_ptr<const FiniteElement> Store;
      static void setStore(Store& store, const FiniteElement& fe)
      { store.reset(new FiniteElement(fe)); }
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

      typedef const FiniteElement *Store;
      static void setStore(Store& store, const FiniteElement& fe)
      { store = &fe; }
    };

    template<class Basis, class = void>
    struct BasisTraits :
      public Basis::Traits
    {
      typedef typename Basis::Traits::DomainLocal DomainLocal;
      typedef typename Basis::Traits::RangeField RangeField;

      template<typename Geometry>
      static void gradient(const Basis& basis, const Geometry& geometry,
                           const DomainLocal& xl,
                           std::vector<FieldMatrix<RangeField, 1,
                               Geometry::coorddimension> >& grad)
      {
        grad.resize(basis.size());
        basis.evaluateJacobian(xl, grad);
      }
    };

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

      template<typename Geometry>
      static void gradient(const Basis& basis, const Geometry& geometry,
                           const DomainLocal& xl,
                           std::vector<FieldMatrix<RangeField, 1,
                               Geometry::coorddimension> >& grad)
      {
        std::vector<Jacobian> lgrad(basis.size());
        basis.evaluateJacobian(xl, lgrad);

        const typename Geometry::Jacobian& jac =
          geometry.jacobianInverseTransposed(xl);

        grad.resize(basis.size());
        for(std::size_t i = 0; i < basis.size(); ++i)
          jac.mv(lgrad[i][0], grad[i][0]);
      }
    };

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_FINITEELEMENT_TRAITS_HH
