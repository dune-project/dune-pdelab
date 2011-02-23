// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_QkDG_LOCALFINITEELEMENT_HH
#define DUNE_QkDG_LOCALFINITEELEMENT_HH

#include <dune/common/geometrytype.hh>
#include <dune/common/fvector.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localtoglobaladaptors.hh>

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
      int size () const
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
              x[j] = (1.0*alpha[j])/k;
	  
            f.evaluate(x,y); out[i] = y;
          }
      }
    };

  }

  /** \todo Please doc me !
   */
  template<class D, class R, int k, int d>
  class QkDGLocalFiniteElement
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
    QkDGLocalFiniteElement ()
    {
      gt.makeCube(d);
    }

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

    QkDGLocalFiniteElement* clone () const
    {
      return new QkDGLocalFiniteElement(*this);
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
    QkDGLocalFiniteElement<
      typename Geometry::ctype, RF, k, Geometry::mydimension
      >,
    Geometry
    >
  {
    typedef QkDGLocalFiniteElement<
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


  /** \todo Please doc me !
   * The general class. We only implement some cases below
   */
  template<class D, class R, int k, int d>
  class QkCGLocalFiniteElement {};


  /** \todo Please doc me !
   * The case k=1
   */
  template<class D, class R, int d>
  class QkCGLocalFiniteElement<D,R,1,d>
  {
    enum{ k = 1 };
    typedef QkStuff::QkLocalBasis<D,R,k,d> LocalBasis;

    class QkCGLocalCoefficients 
    {
      enum{ k = 1 };
    public:
      //! \brief Standard constructor
      QkCGLocalCoefficients () : li(QkStuff::QkSize<k,d>::value)
      {
        for (std::size_t i=0; i<QkStuff::QkSize<k,d>::value; i++)
          li[i] = LocalKey(i,d,0);
      }
      
      //! number of coefficients
      std::size_t size () const
      {
        return QkStuff::QkSize<k,d>::value;
      }
      
      //! get i'th index
      const LocalKey& localKey (std::size_t i) const
      {
        return li[i];
      } 
      
    private:
      std::vector<LocalKey> li;
    };

    typedef QkStuff::QkLocalInterpolation<k,d,LocalBasis> LocalInterpolation;

  public:
    // static number of basis functions
    enum{ n = QkStuff::QkSize<k,d>::value };

    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<LocalBasis,QkCGLocalCoefficients,LocalInterpolation> Traits;

    /** \todo Please doc me !
     */
    QkCGLocalFiniteElement ()
    {
      gt.makeCube(d);
    }

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

    QkCGLocalFiniteElement* clone () const
    {
      return new QkCGLocalFiniteElement(*this);
    }

  private:
    LocalBasis basis;
    QkCGLocalCoefficients coefficients;
    LocalInterpolation interpolation;
    GeometryType gt;
  };

  /** \todo Please doc me !
   * The case k=2, d=2
   */
  template<class D, class R>
  class QkCGLocalFiniteElement<D,R,2,2>
  {
    enum{ k = 2 };
    enum{ d = 2 };
    typedef QkStuff::QkLocalBasis<D,R,k,d> LocalBasis;

    class QkCGLocalCoefficients 
    {
      enum{ k = 2 };
      enum{ d = 2 };
    public:
      //! \brief Standard constructor
      QkCGLocalCoefficients () : li(QkStuff::QkSize<k,d>::value)
      {
        li[0] = LocalKey(0,2,0);
        li[1] = LocalKey(2,1,0);
        li[2] = LocalKey(1,2,0);
        li[3] = LocalKey(0,1,0);
        li[4] = LocalKey(0,0,0);
        li[5] = LocalKey(1,1,0);
        li[6] = LocalKey(2,2,0);
        li[7] = LocalKey(3,1,0);
        li[8] = LocalKey(3,2,0);
      }
      
      //! number of coefficients
      std::size_t size () const
      {
        return QkStuff::QkSize<k,d>::value;
      }
      
      //! get i'th index
      const LocalKey& localKey (std::size_t i) const
      {
        return li[i];
      } 
      
    private:
      std::vector<LocalKey> li;
    };

    typedef QkStuff::QkLocalInterpolation<k,d,LocalBasis> LocalInterpolation;

  public:
    // static number of basis functions
    enum{ n = QkStuff::QkSize<k,d>::value };

    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<LocalBasis,QkCGLocalCoefficients,LocalInterpolation> Traits;

    /** \todo Please doc me !
     */
    QkCGLocalFiniteElement ()
    {
      gt.makeCube(d);
    }

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

    QkCGLocalFiniteElement* clone () const
    {
      return new QkCGLocalFiniteElement(*this);
    }

  private:
    LocalBasis basis;
    QkCGLocalCoefficients coefficients;
    LocalInterpolation interpolation;
    GeometryType gt;
  };


  /** \todo Please doc me !
   * The case k=2, d=3
   */
  template<class D, class R>
  class QkCGLocalFiniteElement<D,R,2,3>
  {
    enum{ k = 2 };
    enum{ d = 3 };
    typedef QkStuff::QkLocalBasis<D,R,k,d> LocalBasis;

    class QkCGLocalCoefficients 
    {
      enum{ k = 2 };
      enum{ d = 3 };
    public:
      //! \brief Standard constructor
      QkCGLocalCoefficients () : li(QkStuff::QkSize<k,d>::value)
      {
        li[24] = LocalKey(6,3,0); li[25] = LocalKey(11,2,0); li[26] = LocalKey(7,3,0);
        li[21] = LocalKey(8,2,0); li[22] = LocalKey(5,1,0);  li[23] = LocalKey(9,2,0);
        li[18] = LocalKey(4,3,0); li[19] = LocalKey(10,2,0); li[20] = LocalKey(5,3,0);

        li[15] = LocalKey(2,2,0); li[16] = LocalKey(3,1,0); li[17] = LocalKey(3,2,0);
        li[12] = LocalKey(0,1,0); li[13] = LocalKey(0,0,0); li[14] = LocalKey(1,1,0);
        li[ 9] = LocalKey(0,2,0); li[10] = LocalKey(2,1,0); li[11] = LocalKey(1,2,0);

        li[ 6] = LocalKey(2,3,0); li[ 7] = LocalKey(7,2,0); li[ 8] = LocalKey(3,3,0);
        li[ 3] = LocalKey(4,2,0); li[ 4] = LocalKey(4,1,0); li[ 5] = LocalKey(5,2,0);
        li[ 0] = LocalKey(0,3,0); li[ 1] = LocalKey(6,2,0); li[ 2] = LocalKey(1,3,0);
      }
      
      //! number of coefficients
      std::size_t size () const
      {
        return QkStuff::QkSize<k,d>::value;
      }
      
      //! get i'th index
      const LocalKey& localKey (std::size_t i) const
      {
        return li[i];
      } 
      
    private:
      std::vector<LocalKey> li;
    };

    typedef QkStuff::QkLocalInterpolation<k,d,LocalBasis> LocalInterpolation;

  public:
    // static number of basis functions
    enum{ n = QkStuff::QkSize<k,d>::value };

    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<LocalBasis,QkCGLocalCoefficients,LocalInterpolation> Traits;

    /** \todo Please doc me !
     */
    QkCGLocalFiniteElement ()
    {
      gt.makeCube(d);
    }

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

    QkCGLocalFiniteElement* clone () const
    {
      return new QkCGLocalFiniteElement(*this);
    }

  private:
    LocalBasis basis;
    QkCGLocalCoefficients coefficients;
    LocalInterpolation interpolation;
    GeometryType gt;
  };

}

namespace Dune {
  namespace PDELab {

    //! wrap up element from local functions
    //! \ingroup FiniteElementMap
    template<class D, class R, int k, int d>
    class QkDGLocalFiniteElementMap
      : public Dune::PDELab::SimpleLocalFiniteElementMap< Dune::QkDGLocalFiniteElement<D,R,k,d> >
    {};

    //! wrap up element from local functions
    //! \ingroup FiniteElementMap
    template<class D, class R, int k, int d>
    class QkCGLocalFiniteElementMap
      : public Dune::PDELab::SimpleLocalFiniteElementMap< Dune::QkCGLocalFiniteElement<D,R,k,d> >
    {};

  }
}

#endif
