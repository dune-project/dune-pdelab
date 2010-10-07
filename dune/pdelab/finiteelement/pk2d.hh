// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_FINITEELEMENT_PK2D_HH
#define DUNE_PDELAB_FINITEELEMENT_PK2D_HH

#include <cstddef>

#include <dune/common/geometrytype.hh>

#include <dune/localfunctions/lagrange/pk2d/pk2dlocalbasis.hh>
#include <dune/localfunctions/lagrange/pk2d/pk2dlocalcoefficients.hh>
#include <dune/localfunctions/lagrange/pk2d/pk2dlocalinterpolation.hh>

#include <dune/pdelab/finiteelement/localtoglobaladaptors.hh>

namespace Dune {
  namespace PDELab {

    //! Langrange finite element of arbitrary order on triangles
    /**
     * \tparam Geometry Geometry for the local to global transformation.
     * \tparam RF       Field type of the range.
     * \tparam k        Maximum polynomial order of the base functions.
     *
     * \implements FiniteElementInterface
     */
    template<class Geometry, class RF, std::size_t k>
    class Pk2DFiniteElement {
      typedef typename Geometry::ctype DF;
      typedef Pk2DLocalBasis<DF,RF,k> LocalBasis;
      typedef Pk2DLocalInterpolation<LocalBasis> LocalInterpolation;

    public:
      /**
       * \implements FiniteElementInterface::Traits
       */
      struct Traits {
        typedef ScalarLocalToGlobalBasisAdaptor<LocalBasis, Geometry> Basis;
        typedef LocalToGlobalInterpolationAdaptor<
          LocalInterpolation,
          typename Basis::Traits
          > Interpolation;
        typedef Pk2DLocalCoefficients<k> Coefficients;
      };

    private:
      static const GeometryType gt;
      static const LocalBasis localBasis;
      static const LocalInterpolation localInterpolation;

      typename Traits::Basis basis_;
      typename Traits::Interpolation interpolation_;
      typename Traits::Coefficients coefficients_;

    public:
      //! construct a Pk2DFiniteElement
      /**
       * \param geometry    The geometry object to use for adaption.

       * \param vertexOrder The global ordering of the vertices within the
       *                    grid, used to determine orientation of the edges.
       *                    This vertexOrder object must support codim=0.
       *
       * \note This class stores the reference to the geometry passed here.
       *       Any use of this class after this references has become invalid
       *       results in undefined behaviour.  The exception is that the
       *       destructor of this class may still be called.  The information
       *       contained in the vertexOrder object is extracted and the object
       *       is no longer needed after the contructor returns.
       */
      template<class VertexOrder>
      Pk2DFiniteElement(const Geometry &geometry,
                        const VertexOrder& vertexOrder) :
        basis_(localBasis, geometry), interpolation_(localInterpolation),
        coefficients_(vertexOrder.begin(0, 0))
      { }

      const typename Traits::Basis& basis() const { return basis_; }
      const typename Traits::Interpolation& interpolation() const
      { return interpolation_; }
      const typename Traits::Coefficients& coefficients() const
      { return coefficients_; }
      const GeometryType &type() const { return gt; }
    };

    template<class Geometry, class RF, std::size_t k>
    const GeometryType
    Pk2DFiniteElement<Geometry, RF, k>::gt(GeometryType::simplex, 2);

    template<class Geometry, class RF, std::size_t k>
    const typename Pk2DFiniteElement<Geometry, RF, k>::LocalBasis
    Pk2DFiniteElement<Geometry, RF, k>::localBasis = LocalBasis();

    template<class Geometry, class RF, std::size_t k>
    const typename Pk2DFiniteElement<Geometry, RF, k>::LocalInterpolation
    Pk2DFiniteElement<Geometry, RF, k>::localInterpolation =
      LocalInterpolation();

    //! Factory for Pk2DFiniteElement objects
    /**
     * Constructs Pk2DFiniteElement objects given a geometry and a vertex
     * ordering.
     *
     * \tparam Geometry Geometry for the local to global transformation.
     * \tparam RF       Field type of the range.
     * \tparam k        Maximum polynomial order of the base functions.
     *
     * \implements FiniteElementFactoryInterface
     */
    template<class Geometry, class RF, std::size_t k>
    struct Pk2DFiniteElementFactory {
      typedef Pk2DFiniteElement<Geometry, RF, k> FiniteElement;

      //! construct Pk2DFiniteElementFactory
      /**
       * \param geometry    The geometry object to use for adaption.
       * \param vertexOrder The global ordering of the vertices within the
       *                    grid, used to determine orientation of the edges.
       *                    This vertexOrder object must support codim=0.
       *
       * \note The returned object stores the reference to the geometry passed
       *       here.  Any use of the returned value after this references has
       *       become invalid results in undefined behaviour.  The exception
       *       is that the destructor of this class may still be called.  The
       *       information contained in the vertexOrder object is extracted
       *       and the object is no longer needed after the contructor
       *       returns.  No reference to internal data of the factory is
       *       stored.
       */
      template<class VertexOrder>
      const FiniteElement make(const Geometry& geometry,
                               const VertexOrder& vertexOrder)
      { return FiniteElement(geometry, vertexOrder); }
    };

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_FINITEELEMENT_PK2D_HH
