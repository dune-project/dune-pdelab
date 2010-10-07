// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=8 sw=2 et sts=2:

#ifndef DUNE_PDELAB_FINITEELEMENT_EDGES0_5_HH
#define DUNE_PDELAB_FINITEELEMENT_EDGES0_5_HH

#include <cstddef>
#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/geometrytype.hh>

#include <dune/localfunctions/common/localkey.hh>
#include <dune/localfunctions/lagrange/p1/p1localbasis.hh>

#include <dune/grid/common/genericreferenceelements.hh>

#include <dune/pdelab/finiteelement/localtoglobaladaptors.hh>

namespace Dune {
  namespace PDELab {

    template<std::size_t dim, class DF = double>
    struct EdgeS0_5Common {
      static const GenericReferenceElement<DF, dim>& refelem;
      static const std::size_t s;
    };

    template<std::size_t dim, class DF>
    const GenericReferenceElement<DF, dim>& EdgeS0_5Common<dim,DF>::
    refelem(GenericReferenceElements<DF, dim>::simplex());

    template<std::size_t dim, typename DF>
    const std::size_t EdgeS0_5Common<dim,DF>::s(refelem.size(dim-1));

    //////////////////////////////////////////////////////////////////////
    //
    //  LocalBasis
    //

    //! LocalBasis for order 0.5 (lowest order) edge elements on simplices
    /**
     * @ingroup LocalBasisImplementation
     *
     * \tparam Geometry Type of the local-to-global map.
     * \tparam RF       Type to represent the field in the range.
     *
     * \nosubgrouping
     */
    template<class Geometry, class RF>
    class EdgeS0_5Basis :
      private EdgeS0_5Common<Geometry::mydimension, typename Geometry::ctype>
    {
    public:
      //! \brief export type traits for function signature
      struct Traits {
        typedef typename Geometry::ctype DomainField;
        static const std::size_t dimDomainLocal = Geometry::mydimension;
        static const std::size_t dimDomainGlobal = Geometry::coorddimension;
        typedef FieldVector<DomainField, dimDomainLocal> DomainLocal;
        typedef FieldVector<DomainField, dimDomainGlobal> DomainGlobal;

        typedef RF RangeField;
        static const std::size_t dimRange = dimDomainLocal;
        typedef FieldVector<RangeField, dimRange> Range;

        typedef FieldMatrix<RangeField, dimRange, dimDomainGlobal> Jacobian;

        static const std::size_t diffOrder = 1;
      };

    private:
      typedef Dune::P1LocalBasis<
        typename Traits::DomainField,
        typename Traits::RangeField,
        Traits::dimDomainLocal
        > P1LocalBasis;
      typedef ScalarLocalToGlobalBasisAdaptor<P1LocalBasis, Geometry> P1Basis;

      static const P1LocalBasis& p1LocalBasis;
      static const std::size_t dim = Traits::dimDomainLocal;

      typedef EdgeS0_5Common<dim, typename Geometry::ctype> Base;
      using Base::refelem;
      using Base::s;

      // global values of the Jacobians (gradients) of the p1 basis
      std::vector<typename P1Basis::Traits::Jacobian> p1j;
      // edge sizes and orientations
      std::vector<typename Traits::DomainField> edgel;

    public:
      //! Construct an EdgeS0_5Basis
      /**
       * \param geo         Geometry of the element to contruct a local basis
       *                    for.
       * \param vertexOrder Vertex ordering information.  Only the vertex
       *                    order on the dim=1 sub-entities (edges) is
       *                    required.
       */
      template<typename VertexOrder>
      EdgeS0_5Basis(const Geometry& geo, const VertexOrder& vertexOrder) :
        p1j(s, typename P1Basis::Traits::Jacobian(0)), edgel(s)
      {
        // use some arbitrary position to evaluate jacobians, they are
        // constant
        static const typename Traits::DomainLocal xl(0);

        // precompute Jacobian (gradients) of the p1 element
        P1Basis(p1LocalBasis, geo).evaluateJacobian(xl, p1j);

        // calculate edge sizes and orientations
        for(std::size_t i = 0; i < s; ++i) {
          edgel[i] = (geo.corner(refelem.subEntity(i,dim-1,0,dim))-
                      geo.corner(refelem.subEntity(i,dim-1,1,dim))
                      ).two_norm();
          const typename VertexOrder::iterator& edgeVertexOrder =
            vertexOrder.begin(dim-1, i);
          if(edgeVertexOrder[0] > edgeVertexOrder[1])
            edgel[i] *= -1;
        }
      }

      //! number of shape functions
      std::size_t size () const { return s; }

      //! Evaluate all shape functions
      void evaluateFunction(const typename Traits::DomainLocal& xl,
                            std::vector<typename Traits::Range>& out) const
      {
        out.assign(s, typename Traits::Range(0));

        // compute p1 values -- use the local basis directly for that, local
        // and global values are identical for scalars
        std::vector<typename P1LocalBasis::Traits::RangeType> p1v;
        p1LocalBasis.evaluateFunction(xl, p1v);

        for(std::size_t i = 0; i < s; i++) {
          const std::size_t i0 = refelem.subEntity(i,dim-1,0,dim);
          const std::size_t i1 = refelem.subEntity(i,dim-1,1,dim);
          out[i].axpy( p1v[i0], p1j[i1][0]);
          out[i].axpy(-p1v[i1], p1j[i0][0]);
          out[i] *= edgel[i];
        }
      }

      //! Evaluate all Jacobians
      void evaluateJacobian(const typename Traits::DomainLocal&,
                            std::vector<typename Traits::Jacobian>& out)
        const
      {
        out.resize(s);

        for(std::size_t i = 0; i < s; i++) {
          const std::size_t i0 = refelem.subEntity(i,dim-1,0,dim);
          const std::size_t i1 = refelem.subEntity(i,dim-1,1,dim);
          for(std::size_t j = 0; j < dim; j++)
            for(std::size_t k = 0; k < dim; k++)
              out[i][j][k] = edgel[i] *
                (p1j[i0][0][k]*p1j[i1][0][j]-p1j[i1][0][k]*p1j[i0][0][j]);
        }
      }

      //! Polynomial order of the shape functions
      std::size_t order () const { return 1; }
    };

    template<class Geometry, class RF>
    const typename EdgeS0_5Basis<Geometry, RF>::P1LocalBasis&
    EdgeS0_5Basis<Geometry, RF>::p1LocalBasis = P1LocalBasis();

    //////////////////////////////////////////////////////////////////////
    //
    // Coefficients
    //

    //! Coefficients for lowest order edge elements on simplices
    /**
     * \nosubgrouping
     * \implements CoefficientsInterface
     *
     * \tparam dim Dimension of both domain and range.
     */
    template<std::size_t dim>
    class EdgeS0_5Coefficients :
      private EdgeS0_5Common<dim>
    {
      using EdgeS0_5Common<dim>::s;

      std::vector<LocalKey> li;

    public:
      //! Standard constructor
      EdgeS0_5Coefficients() : li(s)
      {
        for(std::size_t i = 0; i < s; i++)
          li[i] = LocalKey(i, dim-1, 0);
      }

      //! number of coefficients
      std::size_t size () const { return s; }

      //! get i'th index
      const LocalKey& localKey(std::size_t i) const { return li[i]; }
    };

    //////////////////////////////////////////////////////////////////////
    //
    // Interpolation
    //

    //! Interpolation for lowest order edge elements on simplices
    /**
     * \tparam Geometry Type of the local-to-global map.
     * \tparam RF       Type to represent the field in the range.
     *
     * \nosubgrouping
     */
    template<class Geometry, class Traits_>
    class EdgeS0_5Interpolation :
      private EdgeS0_5Common<Traits_::dimDomainLocal,
                            typename Traits_::DomainField>
    {
    public:
      typedef Traits_ Traits;

    private:
      static const std::size_t dim = Traits::dimDomainLocal;
      typedef EdgeS0_5Common<dim, typename Traits::DomainField> Base;
      using Base::refelem;
      using Base::s;

      std::vector<typename Traits::DomainGlobal> edgev;

    public:
      //! constructor
      /**
       * \param geo         Geometry of the element to contruct a local basis
       *                    for.
       * \param vertexOrder Vertex ordering information.  Only the vertex
       *                    order on the dim=1 sub-entities (edges) is
       *                    required.
       */
      template<typename VertexOrder>
      EdgeS0_5Interpolation(const Geometry& geo,
                            const VertexOrder& vertexOrder)
        : edgev(s)
      {
        for(std::size_t i = 0; i < s; ++i) {
          const std::size_t i0 = refelem.subEntity(i,dim-1,0,dim);
          const std::size_t i1 = refelem.subEntity(i,dim-1,1,dim);

          edgev[i] = geo.corner(i1);
          edgev[i] -= geo.corner(i0);
          edgev[i] /= edgev[i].two_norm();

          const typename VertexOrder::iterator& edgeVertexOrder =
            vertexOrder.begin(dim-1, i);
          if(edgeVertexOrder[0] > edgeVertexOrder[1])
            edgev[i] *= -1;
        }
      }

      //! Interpolation of a function
      template<typename F, typename C>
      void interpolate(const F& f, std::vector<C>& out) const {
        typename Traits::Range y;

        out.resize(s);

        for(std::size_t i = 0; i < s; ++i) {
          f.evaluate(refelem.position(i,dim-1), y);

          out[i] = y * edgev[i];
        }
      }
    };

    //////////////////////////////////////////////////////////////////////
    //
    //  FiniteElement
    //

    //! FiniteElement for lowest order edge elements on simplices
    /**
     * Uses the representation
     * \f[
     *    \mathbf N^i=(L^{i_0}\nabla L^{i_1}-
     *                 L^{i_1}\nabla L^{i_0})\ell^i
     * \f]
     * where \f$L^k\f$ is the P1 shape function for vertex \f$k\f$,
     * \f$i_0\f$ and \f$i_1\f$ are the indices of the vertices of edge
     * \f$i\f$ and \f$\ell^i\f$ is the length of edge \f$i\f$.
     *
     * \tparam D   Type to represent the field in the domain.
     * \tparam R   Type to represent the field in the range.
     * \tparam dim Dimension of both domain and range.
     *
     * \nosubgrouping
     */
    template<class Geometry, class RF>
    class EdgeS0_5FiniteElement
    {
    public:
      /**
       * \implements FiniteElementInterface::Traits
       */
      struct Traits {
        typedef EdgeS0_5Basis<Geometry, RF> Basis;
        typedef EdgeS0_5Interpolation<Geometry,
                                      typename Basis::Traits> Interpolation;
        typedef EdgeS0_5Coefficients<Geometry::mydimension> Coefficients;
      };

    private:
      typename Traits::Basis basis_;
      typename Traits::Interpolation interpolation_;
      static const typename Traits::Coefficients& coefficients_;
      static const GeometryType gt;

    public:
      //! Constructor
      /**
       * \param geo
       * \param orient Orientation of this element
       */
      template<class VertexOrder>
      EdgeS0_5FiniteElement(const Geometry& geo,
                            const VertexOrder& vertexOrder) :
        basis_(geo, vertexOrder), interpolation_(geo, vertexOrder)
      { }

      //! return reference to the basis object
      const typename Traits::Basis& basis() const { return basis_; }
      //! return reference to the interpolation object
      const typename Traits::Interpolation& interpolation() const
      { return interpolation_; }
      //! return reference to the coefficients object
      const typename Traits::Coefficients& coefficients() const
      { return coefficients_; }
      //! return geometry type of this element
      const GeometryType& type() const { return gt; }
    };

    template<class Geometry, class RF>
    const typename EdgeS0_5FiniteElement<Geometry, RF>::Traits::Coefficients&
    EdgeS0_5FiniteElement<Geometry, RF>::coefficients_ =
      typename Traits::Coefficients();

    template<class Geometry, class RF>
    const GeometryType
    EdgeS0_5FiniteElement<Geometry, RF>::gt(GeometryType::simplex,
                                            Geometry::mydimension);

    ////////////////////////////////////////////////////////////////////////
    //
    // Factory
    //

    //! Factory for EdgeS0_5FiniteElement objects
    /**
     * Constructs EdgeS0_5FiniteElement objects given a geometry and a vertex
     * ordering.
     *
     * \tparam Geometry Geometry for the local to global transformation.
     * \tparam RF       Field type of the range.
     *
     * \implements FiniteElementFactoryInterface
     */
    template<class Geometry, class RF>
    struct EdgeS0_5FiniteElementFactory {
      typedef EdgeS0_5FiniteElement<Geometry, RF> FiniteElement;

      //! construct the factory
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

#endif // DUNE_PDELAB_FINITEELEMENT_EDGES0_5_HH
