// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_FINITEELEMENT_SIMPLEWRAPPER_HH
#define DUNE_PDELAB_FINITEELEMENT_SIMPLEWRAPPER_HH

#include <cstddef>
#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/geometrytype.hh>
#include <dune/common/static_assert.hh>

namespace Dune {
  namespace PDELab {

    template<class LocalBasis, class Geometry>
    class SimpleLocalBasisWrapper {
      dune_static_assert(LocalBasis::Traits::dimRange == 1,
                         "SimpleLocalBasisWrapper can only wrap a scalar "
                         "local basis.");

    public:
      struct Traits {
        typedef typename LocalBasis::Traits::DomainFieldType DomainField;
        static const std::size_t dimDomainLocal =
          LocalBasis::Traits::dimDomain;
        static const std::size_t dimDomainGlobal = Geometry::coorddimension;
        typedef typename LocalBasis::Traits::DomainType DomainLocal;
        typedef FieldVector<DomainField, dimDomainGlobal> DomainGlobal;

        typedef typename LocalBasis::Traits::RangeFieldType RangeField;
        static const std::size_t dimRange = LocalBasis::Traits::dimRange;
        typedef typename LocalBasis::Traits::RangeType Range;

        typedef FieldMatrix<RangeField, dimRange, dimDomainGlobal> Jacobian;

        static const std::size_t diffOrder = LocalBasis::Traits::diffOrder;
      };

    private:
      const LocalBasis& localBasis;
      const Geometry& geometry;

    public:
      SimpleLocalBasisWrapper(const LocalBasis& localBasis_,
                              const Geometry& geometry_) :
        localBasis(localBasis_), geometry(geometry_)
      { }

      std::size_t size() const { return localBasis.size(); }
      std::size_t order() const {
        if(geometry.affine())
          // affine linear
          return localBasis.order();
        else
          // assume at most order dim
          return localBasis.order() + Traits::dimDomainGlobal - 1;
      }

      void evaluateFunction(const typename Traits::DomainLocal& in,
                            std::vector<typename Traits::Range>& out) const
      {
        localBasis.evaluateFunction(in, out);
      }

      void evaluateJacobian(const typename Traits::DomainLocal& in,
                            std::vector<typename Traits::Jacobian>& out) const
      {
        std::vector<typename LocalBasis::Traits::JacobianType>
          localJacobian(size());
        localBasis.evaluateJacobian(in, localJacobian);

        const typename Geometry::Jacobian &geoJacobian =
          geometry.jacobianInverseTransposed(in);

        out.resize(size());
        for(std::size_t i = 0; i < size(); ++i)
          geoJacobian.mv(localJacobian[i][0], out[i][0]);
      }
    };

    template<class LocalInterpolation, class Traits_>
    class SimpleLocalInterpolationWrapper {
      const LocalInterpolation& localInterpolation;

    public:
      typedef Traits_ Traits;

      SimpleLocalInterpolationWrapper
      ( const LocalInterpolation& localInterpolation_) :
        localInterpolation(localInterpolation_)
      { }

      template<class Function, class Coeff>
      void interpolate(const Function& function, std::vector<Coeff>& out) const
      {
        localInterpolation.interpolate(function, out);
      }
    };

    template<class LocalFiniteElement, class Geometry>
    struct SimpleLocalFiniteElementWrapper {
      struct Traits {
        typedef SimpleLocalBasisWrapper<typename LocalFiniteElement::Traits::
                                        LocalBasisType, Geometry> Basis;
        typedef SimpleLocalInterpolationWrapper<typename LocalFiniteElement::
                                                Traits::LocalInterpolationType,
                                                typename Basis::Traits>
          Interpolation;
        typedef typename LocalFiniteElement::Traits::LocalCoefficientsType
          Coefficients;
      };

    private:
      const LocalFiniteElement &localFE;
      typename Traits::Basis basis_;
      typename Traits::Interpolation interpolation_;

    public:
      SimpleLocalFiniteElementWrapper(const LocalFiniteElement& localFE_,
                                      const Geometry &geometry) :
        localFE(localFE_),
        basis_(localFE.localBasis(), geometry),
        interpolation_(localFE.localInterpolation())
      { }

      const typename Traits::Basis& basis() const { return basis_; }
      const typename Traits::Interpolation& interpolation() const
      { return interpolation_; }
      const typename Traits::Coefficients& coefficients() const
      { return localFE.localCoefficients(); }
      GeometryType type() const { return localFE.type(); }
    };

    template<class LocalFiniteElement, class Geometry>
    class SimpleLocalFiniteElementWrapperFactory {
      const LocalFiniteElement& localFE;

    public:
      typedef SimpleLocalFiniteElementWrapper<LocalFiniteElement, Geometry>
        FiniteElement;

      SimpleLocalFiniteElementWrapperFactory
      (const LocalFiniteElement &localFE_) : localFE(localFE_) {}

      const FiniteElement make(const Geometry& geometry) {
        return FiniteElement(localFE, geometry);
      }
    };

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_FINITEELEMENT_SIMPLEWRAPPER_HH
