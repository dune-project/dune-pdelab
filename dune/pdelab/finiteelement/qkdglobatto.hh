// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_FINITEELEMENT_QKDGLOBATTO_HH
#define DUNE_PDELAB_FINITEELEMENT_QKDGLOBATTO_HH

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localkey.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localtoglobaladaptors.hh>

namespace Dune
{
  namespace QkStuff
  {

    //! Lagrange polynomials at Gauss-Lobatto points
    template<class D, class R, int k>
    class GaussLobattoLagrangePolynomials
    {
      R xi_gl[k+1];
      R w_gl[k+1];

    public:
      GaussLobattoLagrangePolynomials ()
      {
        int matched_order=-1;
        int matched_size=-1;
        for (int order=1; order<=40; order++)
          {
            const Dune::QuadratureRule<D,1>& rule = Dune::QuadratureRules<D,1>::rule(Dune::GeometryType::cube,order,Dune::QuadratureType::GaussLobatto);
            if (rule.size()==k+1)
              {
                matched_order = order;
                matched_size = rule.size();
                //std::cout << "GL: input order=" << order << " delivered=" << rule.order() << " size=" << rule.size()<< std::endl;
                break;
              }
          }
        if (matched_order<0) DUNE_THROW(Dune::Exception,"could not find Gauss Lobatto rule of appropriate size");
        const Dune::QuadratureRule<D,1>& rule = Dune::QuadratureRules<D,1>::rule(Dune::GeometryType::cube,matched_order,Dune::QuadratureType::GaussLobatto);
        size_t count=0;
        for (typename Dune::QuadratureRule<D,1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            size_t group=count/2;
            size_t member=count%2;
            size_t newj;
            if (member==1) newj=group; else newj=k-group;
            xi_gl[newj] = it->position()[0];
            w_gl[newj] = it->weight();
            count++;
          }
        for (int j=0; j<matched_size/2; j++)
          if (xi_gl[j]>0.5)
            {
              R temp=xi_gl[j];
              xi_gl[j] = xi_gl[k-j];
              xi_gl[k-j] = temp;
              temp=w_gl[j];
              w_gl[j] = w_gl[k-j];
              w_gl[k-j] = temp;
            }
        // for (int i=0; i<k+1; i++)
        //   std::cout << "i=" << i
        //             << " xi=" << xi_gl[i]
        //             << " w=" << w_gl[i]
        //             << std::endl;
      }

      // ith Lagrange polynomial of degree k in one dimension
      R p (int i, D x) const
      {
        R result(1.0);
        for (int j=0; j<=k; j++)
          if (j!=i) result *= (x-xi_gl[j])/(xi_gl[i]-xi_gl[j]);
        return result;
      }

      // derivative of ith Lagrange polynomial of degree k in one dimension
      R dp (int i, D x) const
      {
        R result(0.0);

        for (int j=0; j<=k; j++)
          if (j!=i)
            {
              R prod( 1.0/(xi_gl[i]-xi_gl[j]) );
              for (int l=0; l<=k; l++)
                if (l!=i && l!=j) prod *= (x-xi_gl[l])/(xi_gl[i]-xi_gl[l]);
              result += prod;
            }
        return result;
      }

      // get ith Lagrange point
      R x (int i) const
      {
        return xi_gl[i];
      }

      // get weight of ith Lagrange point
      R w (int i) const
      {
        return w_gl[i];
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
    class QkGLLocalBasis
    {
      enum{ n = QkSize<k,d>::value };
      GaussLobattoLagrangePolynomials<D,R,k> poly;

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
              out[i] *= poly.p(alpha[j],in[j]);
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
                out[i][0][j] = poly.dp(alpha[j],in[j]);

                // rest of the product
                for (int l=0; l<d; l++)
                  if (l!=j)
                    out[i][0][j] *= poly.p(alpha[l],in[l]);
              }
          }
      }

      //! \brief Polynomial order of the shape functions
      unsigned int order () const
      {
        return k;
      }
    };

    /** \todo Please doc me! */
    template<int k, int d, class LB>
    class QkGLLocalInterpolation
    {
      GaussLobattoLagrangePolynomials<double,double,k> poly;

    public:

      //! \brief Local interpolation of a function
      template<typename F, typename C>
      void interpolate (const F& f, std::vector<C>& out) const
      {
        typename LB::Traits::DomainType x;
        typename LB::Traits::RangeType y;

        out.resize(QkSize<k,d>::value);

        for (int i=0; i<QkSize<k,d>::value; i++)
          {
            // convert index i to multiindex
            Dune::FieldVector<int,d> alpha(multiindex<k,d>(i));

            // Generate coordinate of the i-th Lagrange point
            for (int j=0; j<d; j++)
              x[j] = poly.x(alpha[j]);

            f.evaluate(x,y); out[i] = y;
          }
      }
    };

    /** \todo Please doc me! */
    template<int d, class LB>
    class QkGLLocalInterpolation<0,d,LB>
    {
    public:
      //! \brief Local interpolation of a function
      template<typename F, typename C>
      void interpolate (const F& f, std::vector<C>& out) const
      {
        typename LB::Traits::DomainType x(0.5);
        typename LB::Traits::RangeType y;
        f.evaluate(x,y);
        out.resize(1);
        out[0] = y;
      }
    };

  }

  /** \todo Please doc me !
   */
  template<class D, class R, int k, int d>
  class QkDGGLLocalFiniteElement
  {
    typedef QkStuff::QkGLLocalBasis<D,R,k,d> LocalBasis;
    typedef QkStuff::QkDGLocalCoefficients<k,d> LocalCoefficients;
    typedef QkStuff::QkGLLocalInterpolation<k,d,LocalBasis> LocalInterpolation;

  public:
    // static number of basis functions
    enum{ n = QkStuff::QkSize<k,d>::value };

    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<LocalBasis,LocalCoefficients,LocalInterpolation> Traits;

    /** \todo Please doc me !
     */
    QkDGGLLocalFiniteElement ()
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
    GeometryType type () const
    {
      return gt;
    }

    QkDGGLLocalFiniteElement* clone () const
    {
      return new QkDGGLLocalFiniteElement(*this);
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
  class QkDGGLFiniteElementFactory :
    public ScalarLocalToGlobalFiniteElementAdaptorFactory<
    QkDGGLLocalFiniteElement<
      typename Geometry::ctype, RF, k, Geometry::mydimension
      >,
    Geometry
    >
  {
    typedef QkDGGLLocalFiniteElement<
      typename Geometry::ctype, RF, k, Geometry::mydimension
      > LFE;
    typedef ScalarLocalToGlobalFiniteElementAdaptorFactory<LFE, Geometry> Base;

    static const LFE lfe;

  public:
    //! default constructor
    QkDGGLFiniteElementFactory() : Base(lfe) {}
  };

  template<class Geometry, class RF, int k>
  const typename QkDGGLFiniteElementFactory<Geometry, RF, k>::LFE
  QkDGGLFiniteElementFactory<Geometry, RF, k>::lfe;

}

#endif // DUNE_PDELAB_FINITEELEMENT_QKDGLOBATTO_HH
