// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_FINITEELEMENT_QKDGLAGRANGE_HH
#define DUNE_PDELAB_FINITEELEMENT_QKDGLAGRANGE_HH

#include <numeric>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localkey.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localinterpolation.hh>
#include <dune/localfunctions/common/localtoglobaladaptors.hh>
#include <dune/localfunctions/common/localinterpolation.hh>

namespace Dune
{
  namespace QkStuff
  {
    // This is the main class
    // usage: QkSize<2,3>::value
    // k is the polynomial degree,
    // n is the space dimension
    template<int k, int n>
    struct QkSize
    {
      enum{
        value=(k+1)*QkSize<k,n-1>::value
      };
    };

    template<>
    struct QkSize<0,1>
    {
      enum{
        value=1
      };
    };

    template<int k>
    struct QkSize<k,1>
    {
      enum{
        value=k+1
      };
    };

    template<int n>
    struct QkSize<0,n>
    {
      enum{
        value=1
      };
    };

    // ith Lagrange polynomial of degree k in one dimension
    template<class D, class R, int k>
    R p (int i, D x)
    {
      R result(1.0);
      for (int j=0; j<=k; j++)
        if (j!=i) result *= (k*x-j)/(i-j);
      return result;
    }

    // derivative of ith Lagrange polynomial of degree k in one dimension
    template<class D, class R, int k>
    R dp (int i, D x)
    {
      R result(0.0);

      for (int j=0; j<=k; j++)
        if (j!=i)
          {
            R prod( (k*1.0)/(i-j) );
            for (int l=0; l<=k; l++)
              if (l!=i && l!=j) prod *= (k*x-l)/(i-l);
            result += prod;
          }
      return result;
    }

    //! Lagrange polynomials of degree k at equidistant points as a class
    template<class D, class R, int k>
    class EquidistantLagrangePolynomials
    {
    public:
      // ith Lagrange polynomial of degree k in one dimension
      R p (int i, D x) const
      {
        R result(1.0);
        for (int j=0; j<=k; j++)
          if (j!=i) result *= (k*x-j)/(i-j);
        return result;
      }

      // derivative of ith Lagrange polynomial of degree k in one dimension
      R dp (int i, D x) const
      {
        R result(0.0);

        for (int j=0; j<=k; j++)
          if (j!=i)
            {
              R prod( (k*1.0)/(i-j) );
              for (int l=0; l<=k; l++)
                if (l!=i && l!=j) prod *= (k*x-l)/(i-l);
              result += prod;
            }
        return result;
      }

      // get ith Lagrange point
      R x (int i)
      {
        return i/((1.0)*(k+1));
      }
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
    class QkLocalBasis
    {
      enum{ n = QkSize<k,d>::value };
    public:
      typedef LocalBasisTraits<D,d,Dune::FieldVector<D,d>,R,1,Dune::FieldVector<R,1>,Dune::FieldMatrix<R,1,d> > Traits;

      //! \brief number of shape functions
      unsigned int size () const
      {
        return QkSize<k,d>::value;
      }

      //! \brief Evaluate all shape functions
      inline void evaluateFunction (const typename Traits::DomainType& in,
                                    std::vector<typename Traits::RangeType>& out) const
      {
        out.resize(size());
        for (size_t i=0; i<size(); i++)
          {
            // convert index i to multiindex
            Dune::FieldVector<int,d> alpha(multiindex<k,d>(i));

            // initialize product
            out[i] = 1.0;

            // dimension by dimension
            for (int j=0; j<d; j++)
              out[i] *= p<D,R,k>(alpha[j],in[j]);
          }
      }

      //! \brief Evaluate Jacobian of all shape functions
      inline void
      evaluateJacobian (const typename Traits::DomainType& in,         // position
                        std::vector<typename Traits::JacobianType>& out) const      // return value
      {
        out.resize(size());

        // Loop over all shape functions
        for (size_t i=0; i<size(); i++)
          {
            // convert index i to multiindex
            Dune::FieldVector<int,d> alpha(multiindex<k,d>(i));

            // Loop over all coordinate directions
            for (int j=0; j<d; j++)
              {
                // Initialize: the overall expression is a product
                // if j-th bit of i is set to -1, else 1
                out[i][0][j] = dp<D,R,k>(alpha[j],in[j]);

                // rest of the product
                for (int l=0; l<d; l++)
                  if (l!=j)
                    out[i][0][j] *= p<D,R,k>(alpha[l],in[l]);
              }
          }
      }

      //! \brief Evaluate partial derivative of all shape functions
      void partial(const std::array<unsigned int,Traits::dimDomain>& DUNE_UNUSED(order),
                   const typename Traits::DomainType& DUNE_UNUSED(in),
                   std::vector<typename Traits::RangeType>& DUNE_UNUSED(out)) const
      {
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

    /**@ingroup LocalLayoutImplementation
       \brief Layout map for Q1 elements

       \nosubgrouping
       \implements Dune::LocalCoefficientsVirtualImp
    */
    template<int k, int d>
    class QkDGLocalCoefficients
    {
    public:
      //! \brief Standard constructor
      QkDGLocalCoefficients () : li(QkSize<k,d>::value)
      {
        for (std::size_t i=0; i<QkSize<k,d>::value; i++)
          li[i] = LocalKey(0,0,i);
      }

      //! number of coefficients
      std::size_t size () const
      {
        return QkSize<k,d>::value;
      }

      //! get i'th index
      const LocalKey& localKey (std::size_t i) const
      {
        return li[i];
      }

    private:
      std::vector<LocalKey> li;
    };

    /** \todo Please doc me! */
    template<int k, int d, class LB>
    class QkLocalInterpolation
    {
    public:

      //! \brief Local interpolation of a function
      template<typename F, typename C>
      void interpolate (const F& ff, std::vector<C>& out) const
      {
        typename LB::Traits::DomainType x;
        auto&& f = Impl::makeFunctionWithCallOperator<typename LB::Traits::DomainType>(ff);

        out.resize(QkSize<k,d>::value);

        for (int i=0; i<QkSize<k,d>::value; i++)
          {
            // convert index i to multiindex
            Dune::FieldVector<int,d> alpha(multiindex<k,d>(i));

            // Generate coordinate of the i-th Lagrange point
            for (int j=0; j<d; j++)
              x[j] = (1.0*alpha[j])/k;

            out[i] = f(x);
          }
      }
    };

    /** \todo Please doc me! */
    template<int d, class LB>
    class QkLocalInterpolation<0,d,LB>
    {
    public:
      //! \brief Local interpolation of a function
      template<typename F, typename C>
      void interpolate (const F& ff, std::vector<C>& out) const
      {
        typename LB::Traits::DomainType x(0.5);
        auto&& f = Impl::makeFunctionWithCallOperator<typename LB::Traits::DomainType>(ff);

        out.resize(1);
        out[0] = f(x);
      }
    };

  }

  /** \todo Please doc me !
   */
  template<class D, class R, int k, int d>
  class QkDGLagrangeLocalFiniteElement
  {
    typedef QkStuff::QkLocalBasis<D,R,k,d> LocalBasis;
    typedef QkStuff::QkDGLocalCoefficients<k,d> LocalCoefficients;
    typedef QkStuff::QkLocalInterpolation<k,d,LocalBasis> LocalInterpolation;

  public:
    // static number of basis functions
    enum{ n = QkStuff::QkSize<k,d>::value };

    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<LocalBasis,LocalCoefficients,LocalInterpolation> Traits;

    /** \todo Please doc me !
     */
    QkDGLagrangeLocalFiniteElement ()
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

    QkDGLagrangeLocalFiniteElement* clone () const
    {
      return new QkDGLagrangeLocalFiniteElement(*this);
    }

  private:
    LocalBasis basis;
    LocalCoefficients coefficients;
    LocalInterpolation interpolation;
    GeometryType gt;
  };

  //! Factory for global-valued QkDG elements
  /**
   * \tparam Geometry Type of the geometry.  Used to extract the domain field
   *                  type and the dimension.
   * \tparam RF       Range field type.
   */
  template<class Geometry, class RF, int k>
  class QkDGFiniteElementFactory :
    public ScalarLocalToGlobalFiniteElementAdaptorFactory<
    QkDGLagrangeLocalFiniteElement<
      typename Geometry::ctype, RF, k, Geometry::mydimension
      >,
    Geometry
    >
  {
    typedef QkDGLagrangeLocalFiniteElement<
      typename Geometry::ctype, RF, k, Geometry::mydimension
      > LFE;
    typedef ScalarLocalToGlobalFiniteElementAdaptorFactory<LFE, Geometry> Base;

    static const LFE lfe;

  public:
    //! default constructor
    QkDGFiniteElementFactory() : Base(lfe) {}
  };

  template<class Geometry, class RF, int k>
  const typename QkDGFiniteElementFactory<Geometry, RF, k>::LFE
  QkDGFiniteElementFactory<Geometry, RF, k>::lfe;

}


#endif // DUNE_PDELAB_FINITEELEMENT_QKDGLAGRANGE_HH
