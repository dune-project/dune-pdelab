// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PkDG2DLOCALFINITEELEMENT_HH
#define DUNE_PkDG2DLOCALFINITEELEMENT_HH

#include <cstddef>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localtoglobaladaptors.hh>
// #include "pkdg2dlocalbasis.hh"
#include <dune/localfunctions/lagrange/pk2d/pk2dlocalbasis.hh>
#include "pkdgndlocalcoefficients.hh"
#include <dune/localfunctions/lagrange/pk2d/pk2dlocalinterpolation.hh>
#include <dune/localfunctions/lagrange/pk.hh>
// #include "pkdg2dlocalinterpolation.hh"

namespace Dune
{

  /** \todo Please doc me !
   */
  template<class D, class R, unsigned int k>
  class PkDG2DLocalFiniteElement
  {
  public:
    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<Pk2DLocalBasis<D,R,k>,
        Dune::PkDGNDLocalCoefficients<k, 2>,
        Pk2DLocalInterpolation<Pk2DLocalBasis<D,R,k> > > Traits;

    /** \todo Please doc me !
     */
    PkDG2DLocalFiniteElement ()
    {}

    /** \todo Please doc me !
     */
    PkDG2DLocalFiniteElement (int variant) :
      coefficients(variant)
    {}

    /** Constructor for six variants with permuted vertices.

        \param vertexmap The permutation of the vertices.  This
        can for instance be generated from the global indices of
        the vertices by reducing those to the integers 0...2
     */
    PkDG2DLocalFiniteElement (const unsigned int vertexmap[3]) :
      coefficients(vertexmap)
    {}

    /** \todo Please doc me !
     */
    const typename Traits::LocalBasisType& localBasis () const
    {
      return basis;
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalCoefficientsType& localCoefficients () const
    {
      return coefficients;
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalInterpolationType& localInterpolation () const
    {
      return interpolation;
    }

    /** \brief Number of shape functions in this finite element */
    unsigned int size () const
    {
      return basis.size();
    }

    /** \todo Please doc me !
     */
    static constexpr GeometryType type ()
    {
      return GeometryTypes::triangle;
    }

  private:
    Pk2DLocalBasis<D,R,k> basis;
    PkDGNDLocalCoefficients<k, 2> coefficients;
    Pk2DLocalInterpolation<Pk2DLocalBasis<D,R,k> > interpolation;
  };

  //! Langrange finite element of arbitrary order on triangles
  /**
   * \tparam Geometry Geometry for the local to global transformation.
   * \tparam RF       Field type of the range.
   * \tparam k        Maximum polynomial order of the base functions.
   *
   * \implements FiniteElementInterface
   */
  template<class Geometry, class RF, std::size_t k>
  class PkDG2DFiniteElement {
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
      typedef PkDGNDLocalCoefficients<k, 2> Coefficients;
    };

  private:
    static const GeometryType gt;
    static const LocalBasis localBasis;
    static const LocalInterpolation localInterpolation;

    typename Traits::Basis basis_;
    typename Traits::Interpolation interpolation_;
    typename Traits::Coefficients coefficients_;

  public:
    //! construct a PkDG2DFiniteElement
    /**
     * \param geometry    The geometry object to use for adaption.
     * \param vertexOrder The global ordering of the vertices within the grid,
     *                    used to determine orientation of the edges.  This
     *                    vertexOrder object must support codim=0.
     *
     * \note This class stores the reference to the geometry passed here.  Any
     *       use of this class after this references has become invalid
     *       results in undefined behaviour.  The exception is that the
     *       destructor of this class may still be called.  The information
     *       contained in the vertexOrder object is extracted and the object
     *       is no longer needed after the contructor returns.
     */
    template<class VertexOrder>
    PkDG2DFiniteElement(const Geometry &geometry,
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
  PkDG2DFiniteElement<Geometry, RF, k>::gt(GeometryTypes::simplex(2));

  template<class Geometry, class RF, std::size_t k>
  const typename PkDG2DFiniteElement<Geometry, RF, k>::LocalBasis
  PkDG2DFiniteElement<Geometry, RF, k>::localBasis = LocalBasis();

  template<class Geometry, class RF, std::size_t k>
  const typename PkDG2DFiniteElement<Geometry, RF, k>::LocalInterpolation
  PkDG2DFiniteElement<Geometry, RF, k>::localInterpolation =
    LocalInterpolation();

  //! Factory for PkDG2DFiniteElement objects
  /**
   * Constructs PkDG2DFiniteElement objects given a geometry and a vertex
   * ordering.
   *
   * \tparam Geometry Geometry for the local to global transformation.
   * \tparam RF       Field type of the range.
   * \tparam k        Maximum polynomial order of the base functions.
   *
   * \implements FiniteElementFactoryInterface
   */
  template<class Geometry, class RF, std::size_t k>
  struct PkDG2DFiniteElementFactory {
    typedef PkDG2DFiniteElement<Geometry, RF, k> FiniteElement;

    //! construct PkDG2DFiniteElementFactory
    /**
     * \param geometry    The geometry object to use for adaption.
     * \param vertexOrder The global ordering of the vertices within the grid,
     *                    used to determine orientation of the edges.  This
     *                    vertexOrder object must support codim=0.
     *
     * \note The returned object stores the reference to the geometry passed
     *       here.  Any use of the returned value after this references has
     *       become invalid results in undefined behaviour.  The exception is
     *       that the destructor of this class may still be called.  The
     *       information contained in the vertexOrder object is extracted and
     *       the object is no longer needed after the contructor returns.  No
     *       reference to internal data of the factory is stored.
     */
    template<class VertexOrder>
    const FiniteElement make(const Geometry& geometry,
                             const VertexOrder& vertexOrder)
    { return FiniteElement(geometry, vertexOrder); }
  };
}

#endif
