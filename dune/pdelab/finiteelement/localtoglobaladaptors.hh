// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_FINITEELEMENT_LOCALTOGLOBALADAPTORS_HH
#define DUNE_PDELAB_FINITEELEMENT_LOCALTOGLOBALADAPTORS_HH

#include <cstddef>
#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/geometrytype.hh>
#include <dune/common/static_assert.hh>
#include <dune/common/typetraits.hh>

namespace Dune {
  namespace PDELab {

    //! Traits class for local-to-global basis adaptors
    /**
     * \tparam LocalBasisTraits Traits class of the LocalBasis to be adapted.
     * \tparam dimDomainGlobal_ Dimension of the global coordinates,
     *                          i.e. Geometry::coorddimension, if the global
     *                          coordinates are determined by a Geometry.
     *
     * \implements BasisInterface::Traits
     */
    template<class LocalBasisTraits, std::size_t dimDomainGlobal_>
    struct LocalToGlobalBasisAdaptorTraits {
      typedef typename LocalBasisTraits::DomainFieldType DomainField;
      static const std::size_t dimDomainLocal = LocalBasisTraits::dimDomain;
      static const std::size_t dimDomainGlobal = dimDomainGlobal_;
      typedef typename LocalBasisTraits::DomainType DomainLocal;
      typedef FieldVector<DomainField, dimDomainGlobal> DomainGlobal;

      typedef typename LocalBasisTraits::RangeFieldType RangeField;
      static const std::size_t dimRange = LocalBasisTraits::dimRange;
      typedef typename LocalBasisTraits::RangeType Range;

      typedef FieldMatrix<RangeField, dimRange, dimDomainGlobal> Jacobian;

      static const std::size_t diffOrder = LocalBasisTraits::diffOrder;
    };

    //! Convert a simple scalar local basis into a global basis
    /**
     * The local basis must be scalar, i.e. LocalBasis::Traits::dimRange must
     * be 1.  It's values are not transformed.
     *
     * For scalar function \f$f\f$, the gradient is equivalent to the
     * transposed Jacobian \f$\nabla f|_x = J_f^T(x)\f$.  The Jacobian is thus
     * transformed using
     * \f[
     *   \nabla f|_{\mu(\hat x)} =
     *       \hat J_\mu^{-T}(\hat x) \cdot \hat\nabla\hat f|_{\hat x}
     * \f]
     * Here the hat \f$\hat{\phantom x}\f$ denotes local quantities and
     * \f$\mu\f$ denotes the local-to-global map of the geometry.
     *
     * \tparam LocalBasis Type of the local basis to adapt.
     * \tparam Geometry   Type of the local-to-global transformation.
     *
     * \implements BasisInterface
     */
    template<class LocalBasis, class Geometry>
    class ScalarLocalToGlobalBasisAdaptor {
      dune_static_assert(LocalBasis::Traits::dimRange == 1,
                         "ScalarLocalToGlobalBasisAdaptor can only wrap a "
                         "scalar local basis.");
      dune_static_assert((is_same<typename LocalBasis::Traits::DomainFieldType,
                                  typename Geometry::ctype>::value),
                         "ScalarLocalToGlobalBasisAdaptor: LocalBasis must "
                         "use the same ctype as Geometry");
      dune_static_assert
        ( static_cast<std::size_t>(LocalBasis::Traits::dimDomain) ==
            static_cast<std::size_t>(Geometry::mydimension),
          "ScalarLocalToGlobalBasisAdaptor: LocalBasis domain dimension must "
          "match local dimension of Geometry");

      const LocalBasis& localBasis;
      const Geometry& geometry;

    public:
      typedef LocalToGlobalBasisAdaptorTraits<typename LocalBasis::Traits,
                                              Geometry::coorddimension> Traits;

      //! construct a ScalarLocalToGlobalBasisAdaptor
      /**
       * \param localBasis_ The local basis object to adapt.
       * \param geometry_   The geometry object to use for adaption.
       *
       * \note This class stores the references passed here.  Any use of this
       *       class after these references have become invalid results in
       *       undefined behaviour.  The exception is that the destructor of
       *       this class may still be called.
       */
      ScalarLocalToGlobalBasisAdaptor(const LocalBasis& localBasis_,
                                      const Geometry& geometry_) :
        localBasis(localBasis_), geometry(geometry_)
      { }

      std::size_t size() const { return localBasis.size(); }
      //! return maximum polynomial order of the base function
      /**
       * This is to determine the required quadrature order.  For an affine
       * geometry this is the same order as for the local basis.  For other
       * geometries this returns the order of the local basis plus the global
       * dimension minus 1.  The assumtion for non-affine geometries is that
       * they are still multi-linear.
       */
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

    //! Convert a local interpolation into a global interpolation
    /**
     * \tparam LocalInterpolation Type of the local interpolation to adapt.
     * \tparam Traits_            Traits of the corresposnding basis class.
     *
     * \implements InterpolationInterface
     */
    template<class LocalInterpolation, class Traits_>
    class LocalToGlobalInterpolationAdaptor {
      const LocalInterpolation& localInterpolation;

    public:
      typedef Traits_ Traits;

      //! construct a LocalToGlobalInterpolationAdaptor
      /**
       * \param localInterpolation_ The local interpolation object to adapt.
       *
       * \note This class stores the reference to the local interpolation
       *       object passed here.  Any use of this class after the reference
       *       have become invalid results in undefined behaviour.  The
       *       exception is that the destructor of this class may still be
       *       called.
       */
      LocalToGlobalInterpolationAdaptor
      ( const LocalInterpolation& localInterpolation_) :
        localInterpolation(localInterpolation_)
      { }

      template<class Function, class Coeff>
      void interpolate(const Function& function, std::vector<Coeff>& out) const
      {
        localInterpolation.interpolate(function, out);
      }
    };

    //! \brief Convert a simple scalar local finite element into a global
    //!        finite element
    /**
     * The local finite elememt must be scalar,
     * i.e. LocalBasis::Traits::dimRange must be 1.  It's values are not
     * transformed, but the Jacobian is (see ScalarLocalToGlobalBasisAdaptor).
     *
     * \tparam LocalFiniteElement Type of the local finite element to adapt.
     * \tparam Geometry           Type of the local-to-global transformation.
     *
     * \implements FiniteElementInterface
     */
    template<class LocalFiniteElement, class Geometry>
    struct ScalarLocalToGlobalFiniteElementAdaptor {
      /**
       * \implements FiniteElementInterface::Traits
       */
      struct Traits {
        typedef ScalarLocalToGlobalBasisAdaptor<typename LocalFiniteElement::
                Traits::LocalBasisType, Geometry> Basis;
        typedef LocalToGlobalInterpolationAdaptor<typename LocalFiniteElement::
                Traits::LocalInterpolationType, typename Basis::Traits>
                Interpolation;
        typedef typename LocalFiniteElement::Traits::LocalCoefficientsType
                Coefficients;
      };

    private:
      const LocalFiniteElement &localFE;
      typename Traits::Basis basis_;
      typename Traits::Interpolation interpolation_;

    public:
      //! construct a ScalarLocalToGlobalFiniteElementAdaptor
      /**
       * \param localFE_  The local finite element object to adapt.
       * \param geometry_ The geometry object to use for adaption.
       *
       * \note This class stores the references passed here.  Any use of this
       *       class after these references have become invalid results in
       *       undefined behaviour.  The exception is that the destructor of
       *       this class may still be called.
       */
      ScalarLocalToGlobalFiniteElementAdaptor
      ( const LocalFiniteElement& localFE_, const Geometry &geometry) :
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

    //! Factory for ScalarLocalToGlobalFiniteElementAdaptor objects
    /**
     * Constructs ScalarLocalToGlobalFiniteElementAdaptor objects given a
     * local finite element object and a geometry.  This class is restricted
     * to one base variant of the local finite element.
     *
     * \tparam LocalFiniteElement Type of the local finite element to adapt.
     * \tparam Geometry           Type of the local-to-global transformation.
     *
     * \implements FiniteElementFactoryInterface
     */
    template<class LocalFiniteElement, class Geometry>
    class ScalarLocalToGlobalFiniteElementAdaptorFactory {
      const LocalFiniteElement& localFE;

    public:
      typedef ScalarLocalToGlobalFiniteElementAdaptor<LocalFiniteElement,
        Geometry> FiniteElement;

      //! construct a ScalarLocalToGlobalFiniteElementAdaptorFactory
      /**
       * \param localFE_ The local finite element object to adapt.
       *
       * \note This class stores the reference to the local finite element
       *       object passed here.  Any use of this class after this reference
       *       has become invalid results in undefined behaviour.  The
       *       exception is that the destructor of this class may still be
       *       called.
       */
      ScalarLocalToGlobalFiniteElementAdaptorFactory
      (const LocalFiniteElement &localFE_) : localFE(localFE_) {}

      //! construct ScalarLocalToGlobalFiniteElementAdaptor
      /**
       * \param geometry The geometry object to use for adaption.
       *
       * \note The returned object stores the reference to the geometry passed
       *       here as well as references to internal data of this factory.
       *       Any use of the returned value after the geometry reference or
       *       the factory object was become invalid results in undefined
       *       behaviour.  The exception is that the destructor of the
       *       returned value may still be called.
       */
      const FiniteElement make(const Geometry& geometry) {
        return FiniteElement(localFE, geometry);
      }
    };

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_FINITEELEMENT_LOCALTOGLOBALADAPTORS_HH
