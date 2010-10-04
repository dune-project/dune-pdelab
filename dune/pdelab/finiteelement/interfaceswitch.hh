// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_FINITEELEMENT_INTERFACESWITCH_HH
#define DUNE_PDELAB_FINITEELEMENT_INTERFACESWITCH_HH

#include <cstddef>
#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/common/static_assert.hh>
#include <dune/common/typetraits.hh>

namespace Dune {
  namespace PDELab {

    //! \brief Switch for uniform treatment of finite element with either the
    //!        local or the global interface
    /**
     * \tparam FiniteElement Type of the finite element to handle.
     * \tparam Dummy         Dummy parameter for SFINAE.  This must be left at
     *                       the default value of \c void.
     *
     * \note The local interface is detected by the presence of the type
     *       FiniteElement::Traits::LocalBasisType.
     */
    template<class FiniteElement, class Dummy = void>
    struct FiniteElementInterfaceSwitch {
      //! export the type of the basis
      typedef typename FiniteElement::Traits::Basis Basis;
      //! export the type of the interpolation
      typedef typename FiniteElement::Traits::Interpolation Interpolation;
      //! export the type of the coefficients
      typedef typename FiniteElement::Traits::Coefficients Coefficients;

      //! access basis
      static const Basis &basis(const FiniteElement& fe)
      { return fe.basis(); }
      //! access interpolation
      static const Interpolation &interpolation(const FiniteElement& fe)
      { return fe.interpolation(); }
      //! access coefficients
      static const Coefficients &coefficients(const FiniteElement& fe)
      { return fe.coefficients(); }

      //! Type for storing finite elements
      /**
       * Some algorithms use one variable to store (a pointer) a finite
       * element and update that pointer while iterating through the grid.
       * This works well for local finite elements, since they exists in a
       * finite number of variants, which can be stored somewhere and don't
       * need to change for the duration of the algorithm, so we can always
       * store a simple pointer.  For global finite elements we have to store
       * the object itself however, and we must make sure that we destroy the
       * object when we are done with it.  Since global finite elements are
       * not assignable in general, we needs to copy-construct them for each
       * grid element we visit.
       *
       * To accommodate both interfaces, we define a store: for local finite
       * elements it is a simple pointer, and if we want to store a finite
       * element in it we simply store its address.  For global finite
       * elements we use a shared_ptr, and if we want to store a finite
       * element in it we allocate a new object and initialise it with the
       * copy-constructor.  For local finite elements we don't need to do
       * anything when we are done with it, global finite elements are
       * automatically destructed by the shared_ptr when we store a new one or
       * when the shared_ptr itself is destroyed.  Access to the finite
       * element is done by simply dereferencing the store in both cases.
       */
      typedef shared_ptr<const FiniteElement> Store;
      //! Store a finite element in the store.
      /**
       * For local finite elements this means storing the address of the
       * passed reference, for global finite element this means creating a new
       * object with allocation and copy-construction and storing that.
       */
      static void setStore(Store& store, const FiniteElement& fe)
      { store.reset(new FiniteElement(fe)); }
    };

#ifndef DOXYGEN
    //! \brief Switch for uniform treatment of finite element with either the
    //!        local or the global interface
    template<class FiniteElement>
    struct FiniteElementInterfaceSwitch
      < FiniteElement,
        typename enable_if<AlwaysTrue<typename FiniteElement::Traits::
                                      LocalBasisType>::value>::type>
    {
      //! export the type of the basis
      typedef typename FiniteElement::Traits::LocalBasisType Basis;
      //! export the type of the interpolation
      typedef typename FiniteElement::Traits::LocalInterpolationType
        Interpolation;
      //! export the type of the coefficients
      typedef typename FiniteElement::Traits::LocalCoefficientsType
        Coefficients;

      //! access basis
      static const Basis &basis(const FiniteElement& fe)
      { return fe.localBasis(); }
      //! access interpolation
      static const Interpolation &interpolation(const FiniteElement& fe)
      { return fe.localInterpolation(); }
      //! access coefficients
      static const Coefficients &coefficients(const FiniteElement& fe)
      { return fe.localCoefficients(); }

      //! Type for storing finite elements
      typedef const FiniteElement *Store;
      //! Store a finite element in the store.
      static void setStore(Store& store, const FiniteElement& fe)
      { store = &fe; }
    };
#endif // !DOXYGEN

    //! Switch for uniform treatment of local and global basis classes
    /**
     * \tparam Basis Type of the basis to handle.
     * \tparam Dummy Dummy parameter for SFINAE.  This must be left at the
     *               default value of \c void.
     *
     * We don't provide any uniform access to the types and constants
     * pertaining to the global domain.  Providing this would require the
     * Geometry as template parameter as well, and the user code can build
     * them itself if it needs them with the help of the geometry.  The
     * omitted types are \c DomainGlobal and \c Jacobian, the omitted constant
     * is \c dimDomainGlobal.
     *
     * \note The local interface is assumed if the constant
     *       Basis::Traits::dimDomain exists and has a value greater than 0.
     */
    template<class Basis, class Dummy = void>
    struct BasisInterfaceSwitch {
      //! export field types of the coordinates
      typedef typename Basis::Traits::DomainField DomainField;
      //! export dimension of local coordinates
      static const std::size_t dimDomainLocal = Basis::Traits::dimDomainLocal;
      //! export vector type of the local coordinates
      typedef typename Basis::Traits::DomainLocal DomainLocal;

      //! export field type of the values
      typedef typename Basis::Traits::RangeField RangeField;
      //! export dimension of the values
      static const std::size_t dimRange = Basis::Traits::dimRange;
      //! export vector type of the values
      typedef typename Basis::Traits::Range Range;

      //! export number of supported differentiations
      static const std::size_t diffOrder = Basis::Traits::diffOrder;

      //! Compute global gradient for scalar valued bases
      /**
       * \param basis    The basis to get the derivatives from.
       * \param geometry The geometry to use to transform the derivatives (for
       *                 a local basis, unused in the case of a global basis).
       * \param xl       The local coordinates where to evaluate the gradient.
       * \param grad     The result (will be resized to the appropriate number
       *                 of entries.
       *
       * \note This make sense only for a scalar valued basis.
       */
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

#ifndef DOXYGEN
    //! Switch for uniform treatment of local and global basis classes
    template<class Basis>
    struct BasisInterfaceSwitch
    < Basis, typename enable_if<Basis::Traits::dimDomain>::type>
    {
      //! export field types of the coordinates
      typedef typename Basis::Traits::DomainFieldType DomainField;
      //! export dimension of local coordinates
      static const std::size_t dimDomainLocal = Basis::Traits::dimDomain;
      //! export vector type of the local coordinates
      typedef typename Basis::Traits::DomainType DomainLocal;

      //! export field type of the values
      typedef typename Basis::Traits::RangeFieldType RangeField;
      //! export dimension of the values
      static const std::size_t dimRange = Basis::Traits::dimRange;
      //! export vector type of the values
      typedef typename Basis::Traits::RangeType Range;

      //! export number of supported differentiations
      static const std::size_t diffOrder = Basis::Traits::diffOrder;

      //! Compute global gradient for scalar valued bases
      template<typename Geometry>
      static void gradient(const Basis& basis, const Geometry& geometry,
                           const DomainLocal& xl,
                           std::vector<FieldMatrix<RangeField, 1,
                               Geometry::coorddimension> >& grad)
      {
        std::vector<typename Basis::Traits::JacobianType> lgrad(basis.size());
        basis.evaluateJacobian(xl, lgrad);

        const typename Geometry::Jacobian& jac =
          geometry.jacobianInverseTransposed(xl);

        grad.resize(basis.size());
        for(std::size_t i = 0; i < basis.size(); ++i)
          jac.mv(lgrad[i][0], grad[i][0]);
      }
    };
#endif // !DOXYGEN

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_FINITEELEMENT_INTERFACESWITCH_HH
