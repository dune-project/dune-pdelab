// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// DG tensor product basis with Legendre polynomials

#ifndef DUNE_PDELAB_FINITEELEMENT_QKDGLEGENDRE_HH
#define DUNE_PDELAB_FINITEELEMENT_QKDGLEGENDRE_HH

#include <numeric>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/deprecated.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localkey.hh>
#include <dune/localfunctions/common/localtoglobaladaptors.hh>
#include <dune/localfunctions/common/localinterpolation.hh>

namespace Dune
{

  namespace LegendreStuff
  {
    // This is the size class
    // k is the polynomial degree,
    // n is the space dimension
    template<int k, int n>
    struct LegendreSize
    {
      enum{
        value=(k+1)*LegendreSize<k,n-1>::value
      };
    };

    template<>
    struct LegendreSize<0,1>
    {
      enum{
        value=1
      };
    };

    template<int k>
    struct LegendreSize<k,1>
    {
      enum{
        value=k+1
      };
    };

    template<int n>
    struct LegendreSize<0,n>
    {
      enum{
        value=1
      };
    };

    template<int k, int d>
    Dune::FieldVector<int,d> multiindex (int i)
    {
      Dune::FieldVector<int,d> alpha;
      for (int j=0; j<d; j++)
        {
          alpha[j] = i % (k+1);
          i = i/(k+1);
        }
      return alpha;
    }

    //! The 1d Legendre Polynomials (k=0,1 are specialized below)
    template<class D, class R, int k>
    class LegendrePolynomials1d
    {
    public:
      //! evaluate all polynomials at point x
      void p (D x, std::vector<R>& value) const
      {
        value.resize(k+1);
        value[0] = 1;
        value[1] = 2*x-1;
        for (int n=2; n<=k; n++)
          value[n] = ((2*n-1)*(2*x-1)*value[n-1]-(n-1)*value[n-2])/n;
      }

      // ith Lagrange polynomial of degree k in one dimension
      R p (int i, D x) const
      {
        std::vector<R> v(k+1);
        p(x,v);
        return v[i];
      }

      // derivative of all polynomials
      void dp (D x, std::vector<R>& derivative) const
      {
        R value[k+1];
        derivative.resize(k+1);
        value[0] = 1;
        derivative[0] = 0.0;
        value[1] = 2*x-1;
        derivative[1] = 2.0;
        for (int n=2; n<=k; n++)
          {
            value[n] = ((2*n-1)*(2*x-1)*value[n-1]-(n-1)*value[n-2])/n;
            derivative[n] = (2*x-1)*derivative[n-1] + 2*n*value[n-1];
          }
      }

      // value and derivative of all polynomials
      void pdp (D x, std::vector<R>& value, std::vector<R>& derivative) const
      {
        value.resize(k+1);
        derivative.resize(k+1);
        value[0] = 1;
        derivative[0] = 0.0;
        value[1] = 2*x-1;
        derivative[1] = 2.0;
        for (int n=2; n<=k; n++)
          {
            value[n] = ((2*n-1)*(2*x-1)*value[n-1]-(n-1)*value[n-2])/n;
            derivative[n] = (2*x-1)*derivative[n-1] + 2*n*value[n-1];
          }
      }

      // derivative of ith Lagrange polynomial of degree k in one dimension
      R dp (int i, D x) const
      {
        std::vector<R> v(k+1);
        dp(x,v);
        return v[i];
      }
    };

    template<class D, class R>
    class LegendrePolynomials1d<D,R,0>
    {
    public:
      //! evaluate all polynomials at point x
      void p (D x, std::vector<R>& value) const
      {
        value.resize(1);
        value[0] = 1.0;
      }

      // ith Lagrange polynomial of degree k in one dimension
      R p (int i, D x) const
      {
        return 1.0;
      }

      // derivative of all polynomials
      void dp (D x, std::vector<R>& derivative) const
      {
        derivative.resize(1);
        derivative[0] = 0.0;
      }

      // derivative of ith Lagrange polynomial of degree k in one dimension
      R dp (int i, D x) const
      {
        return 0.0;
      }

      // value and derivative of all polynomials
      void pdp (D x, std::vector<R>& value, std::vector<R>& derivative) const
      {
        value.resize(1);
        derivative.resize(1);
        value[0] = 1.0;
        derivative[0] = 0.0;
      }
    };

    template<class D, class R>
    class LegendrePolynomials1d<D,R,1>
    {
    public:
      //! evaluate all polynomials at point x
      void p (D x, std::vector<R>& value) const
      {
        value.resize(2);
        value[0] = 1.0;
        value[1] = 2*x-1;
      }

      // ith Lagrange polynomial of degree k in one dimension
      R p (int i, D x) const
      {
        return (1-i) + i*(2*x-1);
      }

      // derivative of all polynomials
      void dp (D x, std::vector<R>& derivative) const
      {
        derivative.resize(2);
        derivative[0] = 0.0;
        derivative[1] = 2.0;
      }

      // derivative of ith Lagrange polynomial of degree k in one dimension
      R dp (int i, D x) const
      {
        return (1-i)*0 + i*(2);
      }

      // value and derivative of all polynomials
      void pdp (D x, std::vector<R>& value, std::vector<R>& derivative) const
      {
        value.resize(2);
        derivative.resize(2);
        value[0] = 1.0;
        value[1] = 2*x-1;
        derivative[0] = 0.0;
        derivative[1] = 2.0;
      }
    };

    /**@ingroup LocalBasisImplementation
       \brief Lagrange shape functions of order k on the reference cube.

       Also known as \f$Q^k\f$.

       \tparam D Type to represent the field in the domain.
       \tparam R Type to represent the field in the range.
       \tparam k Polynomial degree
       \tparam d Dimension of the cube

       \nosubgrouping
    */
    template<class D, class R, int k, int d>
    class DGLegendreLocalBasis
    {
      enum { n = LegendreSize<k,d>::value };
      LegendrePolynomials1d<D,R,k> poly;
      mutable std::vector<std::vector<R> > v;
      mutable std::vector<std::vector<R> > a;

    public:
      typedef LocalBasisTraits<D,d,Dune::FieldVector<D,d>,R,1,Dune::FieldVector<R,1>,Dune::FieldMatrix<R,1,d> > Traits;

      DGLegendreLocalBasis () : v(d,std::vector<R>(k+1,0.0)), a(d,std::vector<R>(k+1,0.0))
      {}

      //! \brief number of shape functions
      unsigned int size () const
      {
        return n;
      }

      //! \brief Evaluate all shape functions
      inline void evaluateFunction (const typename Traits::DomainType& x,
                                    std::vector<typename Traits::RangeType>& value) const
      {
        // resize output vector
        value.resize(n);

        // compute values of 1d basis functions in each direction
        for (size_t j=0; j<d; j++) poly.p(x[j],v[j]);

        // now compute the product
        for (size_t i=0; i<n; i++)
          {
            // convert index i to multiindex
            Dune::FieldVector<int,d> alpha(multiindex<k,d>(i));

            // initialize product
            value[i] = 1.0;

            // dimension by dimension
            for (int j=0; j<d; j++) value[i] *= v[j][alpha[j]];
          }
      }

      //! \brief Evaluate Jacobian of all shape functions
      inline void
      evaluateJacobian (const typename Traits::DomainType& x,         // position
                        std::vector<typename Traits::JacobianType>& value) const      // return value
      {
        // resize output vector
        value.resize(size());

        // compute values of 1d basis functions in each direction
        for (size_t j=0; j<d; j++) poly.pdp(x[j],v[j],a[j]);

        // Loop over all shape functions
        for (size_t i=0; i<n; i++)
          {
            // convert index i to multiindex
            Dune::FieldVector<int,d> alpha(multiindex<k,d>(i));

            // Loop over all coordinate directions
            for (int j=0; j<d; j++)
              {
                // Initialize: the overall expression is a product
                value[i][0][j] = a[j][alpha[j]];

                // rest of the product
                for (int l=0; l<d; l++)
                  if (l!=j)
                    value[i][0][j] *= v[l][alpha[l]];
              }
          }
      }

      //! \brief Evaluate partial derivative of all shape functions
      void partial(const std::array<unsigned int, Traits::dimDomain>& DUNE_UNUSED(order),
                   const typename Traits::DomainType& DUNE_UNUSED(in),
                   std::vector<typename Traits::RangeType>& DUNE_UNUSED(out)) const {
        auto totalOrder = std::accumulate(order.begin(), order.end(), 0);
        if (totalOrder == 0) {
          evaluateFunction(in, out);
        } else {
          DUNE_THROW(NotImplemented, "Desired derivative order is not implemented");
        }
      }

      //! \brief Polynomial order of the shape functions
      unsigned int order () const
      {
        return k;
      }
    };


    //! determine degrees of freedom
    template<class D, class R, int k, int d>
    class DGLegendreLocalInterpolation
    {
      enum { n = LegendreSize<k,d>::value };
      typedef DGLegendreLocalBasis<D,R,k,d> LB;
      const LB lb;
      Dune::GeometryType gt;

    public:

      DGLegendreLocalInterpolation () : gt(GeometryTypes::cube(d))
      {}

      //! \brief Local interpolation of a function
      template<typename F, typename C>
      void interpolate (const F& ff, std::vector<C>& out) const
      {
        // select quadrature rule
        typedef typename LB::Traits::RangeType RangeType;

        auto&& f = Impl::makeFunctionWithCallOperator<typename LB::Traits::DomainType>(ff);
        const Dune::QuadratureRule<R,d>&
          rule = Dune::QuadratureRules<R,d>::rule(gt,2*k);

        // prepare result
        out.resize(n);
        std::vector<R> diagonal(n);
        for (int i=0; i<n; i++) { out[i] = 0.0; diagonal[i] = 0.0; }

        // loop over quadrature points
        for (typename Dune::QuadratureRule<R,d>::const_iterator
               it=rule.begin(); it!=rule.end(); ++it)
          {
            // evaluate function at quadrature point
            typename LB::Traits::DomainType x;
            RangeType y;
            for (int i=0; i<d; i++) x[i] = it->position()[i];
            y = f(x);

            // evaluate the basis
            std::vector<RangeType> phi(n);
            lb.evaluateFunction(it->position(),phi);

            // do integration
            for (int i=0; i<n; i++) {
              out[i] += y*phi[i]*it->weight();
              diagonal[i] += phi[i]*phi[i]*it->weight();
            }
          }
        for (int i=0; i<n; i++) out[i] /= diagonal[i];
      }
    };

    /**@ingroup LocalLayoutImplementation
       \brief Layout map for Q1 elements

       \nosubgrouping
       \implements Dune::LocalCoefficientsVirtualImp
    */
    template<int k, int d>
    class DGLegendreLocalCoefficients
    {
      enum { n = LegendreSize<k,d>::value };

    public:
      //! \brief Standard constructor
      DGLegendreLocalCoefficients () : li(n)
      {
        for (std::size_t i=0; i<n; i++)
          li[i] = LocalKey(0,0,i);
      }

      //! number of coefficients
      std::size_t size () const
      {
        return n;
      }

      //! get i'th index
      const LocalKey& localKey (std::size_t i) const
      {
        return li[i];
      }

    private:
      std::vector<LocalKey> li;
    };

  } // end of LegendreStuff namespace

  /** \todo Please doc me !
   */
  template<class D, class R, int k, int d>
  class QkDGLegendreLocalFiniteElement
  {
    typedef LegendreStuff::DGLegendreLocalBasis<D,R,k,d> LocalBasis;
    typedef LegendreStuff::DGLegendreLocalCoefficients<k,d> LocalCoefficients;
    typedef LegendreStuff::DGLegendreLocalInterpolation<D,R,k,d> LocalInterpolation;

  public:
    // static number of basis functions
    enum { n = LegendreStuff::LegendreSize<k,d>::value };

    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<LocalBasis,LocalCoefficients,LocalInterpolation> Traits;

    /** \todo Please doc me !
     */
    QkDGLegendreLocalFiniteElement ()
      : gt(GeometryTypes::cube(d))
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

    /** \todo Please doc me !
     */
    std::size_t size() const
    {
      return basis.size();
    }

    /** \todo Please doc me !
     */
    GeometryType type () const
    {
      return gt;
    }

    QkDGLegendreLocalFiniteElement* clone () const
    {
      return new QkDGLegendreLocalFiniteElement(*this);
    }

  private:
    LocalBasis basis;
    LocalCoefficients coefficients;
    LocalInterpolation interpolation;
    GeometryType gt;
  };

  //! Factory for global-valued DGLegendre elements
  /**
   * \tparam Geometry Type of the geometry.  Used to extract the domain field
   *                  type and the dimension.
   * \tparam RF       Range field type.
   */
  template<class Geometry, class RF, int k>
  class QkDGLegendreFiniteElementFactory :
    public ScalarLocalToGlobalFiniteElementAdaptorFactory<
    QkDGLegendreLocalFiniteElement<
      typename Geometry::ctype, RF, k, Geometry::mydimension
      >,
    Geometry
    >
  {
    typedef QkDGLegendreLocalFiniteElement<
      typename Geometry::ctype, RF, k, Geometry::mydimension
      > LFE;
    typedef ScalarLocalToGlobalFiniteElementAdaptorFactory<LFE, Geometry> Base;

    static const LFE lfe;

  public:
    //! default constructor
    QkDGLegendreFiniteElementFactory() : Base(lfe) {}
  };

  template<class Geometry, class RF, int k>
  const typename QkDGLegendreFiniteElementFactory<Geometry, RF, k>::LFE
  QkDGLegendreFiniteElementFactory<Geometry, RF, k>::lfe;

} // End Dune namespace

#endif // DUNE_PDELAB_FINITEELEMENT_QKDGLEGENDRE_HH
