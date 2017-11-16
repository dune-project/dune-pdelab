// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_FINITEELEMENT_L2ORTHONORMAL_HH
#define DUNE_PDELAB_FINITEELEMENT_L2ORTHONORMAL_HH

#include<iostream>
#include<algorithm>
#include<memory>

#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>
#include<dune/common/gmpfield.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/array.hh>

#include<dune/geometry/referenceelements.hh>
#include<dune/geometry/quadraturerules.hh>
#include<dune/geometry/type.hh>

#include<dune/localfunctions/common/localbasis.hh>
#include<dune/localfunctions/common/localkey.hh>
#include<dune/localfunctions/common/localfiniteelementtraits.hh>

/** @defgroup PB Polynomial Basis
 *  @ingroup FiniteElementMap
 *  @{
 */

/** @file
 *  @brief This file defines polynomial basis functions on the reference
 *         element in a generic way
 */

namespace Dune {
  namespace PB {

#if HAVE_GMP
    typedef GMPField<512> DefaultComputationalFieldType;
#else
    typedef double DefaultComputationalFieldType;
#endif

    //=====================================================
    // TMPs for computing number of basis functions in P_k
    //=====================================================

    template<int k, int d> struct PkSize;

    template<int j, int k, int d>
    struct PkSizeHelper
    {
      enum{
        value = PkSizeHelper<j-1,k,d>::value + PkSize<k-j,d-1>::value
      };
    };

    template<int k, int d>
    struct PkSizeHelper<0,k,d>
    {
      enum{
        value = PkSize<k,d-1>::value
      };
    };

    // This is the main class
    // k is the polynomial degree and d is the space dimension
    // then PkSize<k,d> is the number of polynomials of at most total degree k.
    template<int k, int d>
    struct PkSize
    {
      enum{
        value=PkSizeHelper<k,k,d>::value
      };
    };

    template<>
    struct PkSize<0,1>
    {
      enum{
        value=1
      };
    };

    template<int k>
    struct PkSize<k,1>
    {
      enum{
        value=k+1
      };
    };

    template<int d>
    struct PkSize<0,d>
    {
      enum{
        value=1
      };
    };

    // number of polynomials of exactly degree k
    template<int k, int d>
    struct PkExactSize
    {
      enum{
        value=PkSize<k,d>::value-PkSize<k-1,d>::value
      };
    };

    template<int d>
    struct PkExactSize<0,d>
    {
      enum{
        value=1
      };
    };

    //=====================================================
    // TMPs for computing number of basis functions in Q_k
    //=====================================================

    // This is the main class
    // usage: QkSize<2,3>::value
    // k is the polynomial degree, d is the space dimension
    template<int k, int d>
    struct QkSize
    {
      enum{
        value=(k+1)*QkSize<k,d-1>::value
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

    template<int d>
    struct QkSize<0,d>
    {
      enum{
        value=1
      };
    };

    //=====================================================
    // Type to represent a multiindex in d dimensions
    //=====================================================

    template<int d>
    class MultiIndex : public std::array<int,d>
    {
    public:

      MultiIndex () : std::array<int,d>()
      {
      }

    };

    // the number of polynomials of at most degree k in space dimension d (as run-time function)
    inline int pk_size (int k, int d)
    {
      if (k==0) return 1;
      if (d==1) return k+1;
      int n=0;
      for (int j=0; j<=k; j++)
        n += pk_size(k-j,d-1);
      return n;
    }

    // the number of polynomials of exactly degree k in space dimension d (as run-time function)
    inline int pk_size_exact (int k, int d)
    {
      if (k==0)
        return 1;
      else
        return pk_size(k,d)-pk_size(k-1,d);
    }

    // enumerate all multiindices of degree k and find the i'th
    template<int d>
    void pk_enumerate_multiindex (MultiIndex<d>& alpha, int& count, int norm, int dim, int k, int i)
    {
      // check if dim is valid
      if (dim>=d) return;

      // put all k to current dimension and check
      alpha[dim]=k; count++; // alpha has correct norm
      //      std::cout << "dadi alpha=" << alpha << " count=" << count << " norm=" << norm+k << " dim=" << dim << " k=" << k << " i=" << i << std::endl;
      if (count==i) return; // found the index

      // search recursively
      for (int m=k-1; m>=0; m--)
        {
          alpha[dim]=m;
          //std::cout << "dada alpha=" << alpha << " count=" << count << " norm=" << norm+m << " dim=" << dim << " k=" << k << " i=" << i << std::endl;
          pk_enumerate_multiindex(alpha,count,norm+m,dim+1,k-m,i);
          if (count==i) return;
        }

      // reset if index is not yet found
      alpha[dim]=0;
    }

    // map integer 0<=i<pk_size(k,d) to multiindex
    template<int d>
    void pk_multiindex (int i, MultiIndex<d>& alpha)
    {
      for (int j=0; j<d; j++) alpha[j] = 0; // set alpha to 0
      if (i==0) return;                     // handle case i==0 and k==0
      for (int k=1; k<25; k++)
        if (i>=pk_size(k-1,d) && i<pk_size(k,d))
          {
            int count = -1;
            pk_enumerate_multiindex(alpha,count,0,0,k,i-pk_size(k-1,d));
            return;
          }
      DUNE_THROW(Dune::NotImplemented,"Polynomial degree greater than 24 in pk_multiindex");
    }

    // the number of polynomials of at most degree k in space dimension d (as run-time function)
    inline int qk_size (int k, int d)
    {
      int n = 1;
      for (int i=0; i<d; ++i)
        n *= (k+1);
      return n;
    }

    // map integer 0<=i<qk_size(k,d) to multiindex
    template<int d>
    void qk_multiindex (int i, int k, MultiIndex<d>& alpha)
    {
      for (int j = 0; j<d; ++j) {
        alpha[j] = i % (k+1);
        i /= (k+1);
      }
    }

    //=====================================================
    // Traits classes to group Pk and Qk specifics
    //=====================================================
    enum BasisType {
      Pk, Qk
    };

    template <BasisType basisType>
    struct BasisTraits;

    template <>
    struct BasisTraits<BasisType::Pk>
    {
      template <int k, int d>
      struct Size
      {
        enum{
          value = PkSize<k,d>::value
        };
      };

      template <int k, int d>
      struct Order
      {
        enum{
          value = k
        };
      };

      static int size(int k, int d)
      {
        return pk_size(k,d);
      }

      template <int d>
      static void multiindex(int i, int k, MultiIndex<d>& alpha)
      {
        pk_multiindex(i,alpha);
      }
    };

    template <>
    struct BasisTraits<BasisType::Qk>
    {
      template <int k, int d>
      struct Size
      {
        enum{
          value = QkSize<k,d>::value
        };
      };

      template <int k, int d>
      struct Order
      {
        enum{
          // value = d*k
          value = k
        };
      };

      static int size(int k, int d)
      {
        return qk_size(k,d);
      }

      template <int d>
      static void multiindex(int i, int k, MultiIndex<d>& alpha)
      {
        return qk_multiindex(i,k,alpha);
      }
    };

    //=====================================================
    // Integration kernels for monomials
    //=====================================================

    //! compute binomial coefficient "n over k"
    inline long binomial (long n, long k)
    {
      // pick the shorter version of
      // n*(n-1)*...*(k+1)/((n-k)*(n-k-1)*...*1)
      // and
      // n*(n-1)*...*(n-k+1)/(k*(k-1)*...*1)
      if (2*k>=n)
        {
          long nominator=1;
          for (long i=k+1; i<=n; i++) nominator *= i;
          long denominator=1;
          for (long i=2; i<=n-k; i++) denominator *= i;
          return nominator/denominator;
        }
      else
        {
          long nominator=1;
          for (long i=n-k+1; i<=n; i++) nominator *= i;
          long denominator=1;
          for (long i=2; i<=k; i++) denominator *= i;
          return nominator/denominator;
        }
    }

    /** \brief Integrate monomials over the reference element.
     *
     *  \tparam ComputationFieldType    Type to do computations with. Might be high precission.
     *  \tparam GeometryType::BasicType The reference element
     *  \tparam d                       The space dimension.
     */
    template<typename ComputationFieldType, Dune::GeometryType::BasicType bt, int d>
    class MonomialIntegrator
    {
    public:
      //! integrate one monomial
      ComputationFieldType integrate (const MultiIndex<d>& a) const
      {
        DUNE_THROW(Dune::NotImplemented,"non-specialized version of MonomalIntegrator called. Please implement.");
      }
    };

    /** \brief Integrate monomials over the cube in any d.
     */
    template<typename ComputationFieldType, int d>
    class MonomialIntegrator<ComputationFieldType,Dune::GeometryType::cube,d>
    {
    public:
      //! integrate one monomial
      ComputationFieldType integrate (const MultiIndex<d>& a) const
      {
        ComputationFieldType result(1.0);
        for (int i=0; i<d; i++)
          {
            ComputationFieldType exponent(a[i]+1);
            result = result/exponent;
          }
        return result;
      }
    };

    /** \brief Integrate monomials over the unit interval considered as a simplex.
     */
    template<typename ComputationFieldType>
    class MonomialIntegrator<ComputationFieldType,Dune::GeometryType::simplex,1>
    {
    public:
      //! integrate one monomial
      ComputationFieldType integrate (const MultiIndex<1>& a) const
      {
        ComputationFieldType one(1.0);
        ComputationFieldType exponent0(a[0]+1);
        return one/exponent0;
      }
    };

    /** \brief Integrate monomials over the triangle.
     */
    template<typename ComputationFieldType>
    class MonomialIntegrator<ComputationFieldType,Dune::GeometryType::simplex,2>
    {
    public:
      //! integrate one monomial
      ComputationFieldType integrate (const MultiIndex<2>& a) const
      {
        ComputationFieldType sum(0.0);
        for (int k=0; k<=a[1]+1; k++)
          {
            int sign=1;
            if (k%2==1) sign=-1;
            ComputationFieldType nom(sign*binomial(a[1]+1,k));
            ComputationFieldType denom(a[0]+k+1);
            sum = sum + (nom/denom);
          }
        ComputationFieldType denom(a[1]+1);
        return sum/denom;
      }
    };

    /** \brief Integrate monomials over the tetrahedron.
     */
    template<typename ComputationFieldType>
    class MonomialIntegrator<ComputationFieldType,Dune::GeometryType::simplex,3>
    {
    public:
      //! integrate one monomial
      ComputationFieldType integrate (const MultiIndex<3>& a) const
      {
        ComputationFieldType sumk(0.0);
        for (int k=0; k<=a[2]+1; k++)
          {
            int sign=1;
            if (k%2==1) sign=-1;
            ComputationFieldType nom(sign*binomial(a[2]+1,k));
            ComputationFieldType denom(a[1]+k+1);
            sumk = sumk + (nom/denom);
          }
        ComputationFieldType sumj(0.0);
        for (int j=0; j<=a[1]+a[2]+2; j++)
          {
            int sign=1;
            if (j%2==1) sign=-1;
            ComputationFieldType nom(sign*binomial(a[1]+a[2]+2,j));
            ComputationFieldType denom(a[0]+j+1);
            sumj = sumj + (nom/denom);
          }
        ComputationFieldType denom(a[2]+1);
        return sumk*sumj/denom;
      }
    };

    //=====================================================
    // construct orthonormal basis
    //=====================================================

    //! \brief compute \f$ \prod_{i=0}^{d-1} x_i^{a_i} \f$
    template<typename F, int d>
    struct MonomialEvaluate
    {
      template<typename X, typename A>
      static F compute (const X& x, const A& a)
      {
        F prod(1.0);
        for (int i=0; i<a[d]; i++)
          prod = prod*x[d];
        return prod*MonomialEvaluate<F,d-1>::compute(x,a);
      }
    };

    template<typename F>
    struct MonomialEvaluate<F,0>
    {
      template<typename X, typename A>
      static F compute (const X& x, const A& a)
      {
        F prod(1.0);
        for (int i=0; i<a[0]; i++)
          prod = prod*x[0];
        return prod;
      }
    };

    /** \brief Integrate monomials over the reference element.
     *
     * Computes an L_2 orthonormal basis of P_k on the given reference element.
     * The basis polynomials are stored in a monomial representation. With the
     * matrix coeffs private to this class we have
     * \f[
     *     phi_i(x) = \sum_{j=0}{n_k-1} c[i][j] x^{\alpha_j}  \qquad   (1)
     * \f]
     * with n_k     : the dimension of P_k
     *      alpha_j : the exponents of the j-th monomial
     *
     * The class can be used to evaluate polynomials with any degree l smaller
     * or equal to the compile-time parameter k.
     *
     * Calculating derivatives. From (1) we have
     * \f{align*}{
     *     \partial_s \phi_i(x) &= \sum_{j=0}{n_k-1} c[i][j] \partial_s x^{(\alpha_{j1},...,\alpha_{jd})} \\
     *                          &= \sum_{j=0}{n_k-1} c[i][j] * \alpha_js * x^{\beta_j}
     * \f}
     * where
     *       beta_jr = alpha_jr-1 if r=s and alpha_jr else.
     *
     *  \tparam FieldType               Type to represent coefficients after computation.
     *  \tparam k                       The polynomial degreee.
     *  \tparam d                       The space dimension.
     *  \tparam GeometryType::BasicType The reference element
     *  \tparam ComputationFieldType    Type to do computations with. Might be high precission.
     *  \tparam basisType               Type of the polynomial basis. eiter Pk or Qk
     */
    template<typename FieldType, int k, int d, Dune::GeometryType::BasicType bt, typename ComputationFieldType=FieldType, BasisType basisType = BasisType::Pk>
    class OrthonormalPolynomialBasis
    {
      typedef BasisTraits<basisType> Traits;
    public:
      enum{ n = Traits::template Size<k,d>::value };
      typedef Dune::FieldMatrix<FieldType,n,n> LowprecMat;
      typedef Dune::FieldMatrix<ComputationFieldType,n,n> HighprecMat;

      // construct orthonormal basis
      OrthonormalPolynomialBasis ()
        : coeffs(new LowprecMat)
      {
        for (int i=0; i<d; ++i)
          gradcoeffs[i].reset(new LowprecMat());
        // compute index to multiindex map
        for (int i=0; i<n; i++)
          {
            alpha[i].reset(new MultiIndex<d>());
            Traits::multiindex(i,k,*alpha[i]);
            //std::cout << "i=" << i << " alpha_i=" << alpha[i] << std::endl;
          }

        orthonormalize();
      }

      // construct orthonormal basis from an other basis
      template<class LFE>
      OrthonormalPolynomialBasis (const LFE & lfe)
        : coeffs(new LowprecMat)
      {
        for (int i=0; i<d; ++i)
          gradcoeffs[i].reset(new LowprecMat());
        // compute index to multiindex map
        for (int i=0; i<n; i++)
          {
            alpha[i].reset(new MultiIndex<d>());
            Traits::multiindex(i,k,*alpha[i]);
            //std::cout << "i=" << i << " alpha_i=" << alpha[i] << std::endl;
          }

        orthonormalize();
      }

      // return dimension of P_l
      int size (int l)
      {
        return Traits::size(l,d);
      }

      // evaluate all basis polynomials at given point
      template<typename Point, typename Result>
      void evaluateFunction (const Point& x, Result& r) const
      {
        std::fill(r.begin(),r.end(),0.0);
        for (int j=0; j<n; ++j)
          {
            const FieldType monomial_value = MonomialEvaluate<FieldType,d-1>::compute(x,*alpha[j]);
            for (int i=j; i<n; ++i)
              r[i] += (*coeffs)[i][j] * monomial_value;
          }
      }

      // evaluate all basis polynomials at given point
      template<typename Point, typename Result>
      void evaluateJacobian (const Point& x, Result& r) const
      {
        std::fill(r.begin(),r.end(),0.0);

        for (int j=0; j<n; ++j)
          {
            const FieldType monomial_value = MonomialEvaluate<FieldType,d-1>::compute(x,*alpha[j]);
            for (int i=j; i<n; ++i)
              for (int s=0; s<d; ++s)
                r[i][0][s] += (*gradcoeffs[s])[i][j]*monomial_value;
          }
      }

      // evaluate all basis polynomials at given point up to order l <= k
      template<typename Point, typename Result>
      void evaluateFunction (int l, const Point& x, Result& r) const
      {
        if (l>k)
          DUNE_THROW(Dune::RangeError,"l>k in OrthonormalPolynomialBasis::evaluateFunction");

        for (int i=0; i<Traits::size(l,d); i++)
          {
            FieldType sum(0.0);
            for (int j=0; j<=i; j++)
              sum = sum + (*coeffs)[i][j]*MonomialEvaluate<FieldType,d-1>::compute(x,*alpha[j]);
            r[i] = sum;
          }
      }

      // evaluate all basis polynomials at given point
      template<typename Point, typename Result>
      void evaluateJacobian (int l, const Point& x, Result& r) const
      {
        if (l>k)
          DUNE_THROW(Dune::RangeError,"l>k in OrthonormalPolynomialBasis::evaluateFunction");

        for (int i=0; i<Traits::size(l,d); i++)
          {
            FieldType sum[d];
            for (int s=0; s<d; s++)
              {
                sum[s] = 0.0;
                for (int j=0; j<=i; j++)
                  sum[s] += (*gradcoeffs[s])[i][j]*MonomialEvaluate<FieldType,d-1>::compute(x,*alpha[j]);
              }
            for (int s=0; s<d; s++) r[i][0][s] = sum[s];
          }
      }

    private:
      // store multiindices and coefficients on heap
      std::array<std::shared_ptr<MultiIndex<d> >,n> alpha; // store index to multiindex map
      std::shared_ptr<LowprecMat> coeffs; // coefficients with respect to monomials
      std::array<std::shared_ptr<LowprecMat>,d > gradcoeffs; // coefficients of gradient

      // compute orthonormalized shapefunctions from a given set of coefficients
      void orthonormalize()
      {
        // run Gram-Schmidt orthogonalization procedure in high precission
        gram_schmidt();

        // std::cout << "orthogonal basis monomial representation" << std::endl;
        // for (int i=0; i<n; i++)
        //   {
        //     std::cout << "phi_" << i << ":" ;
        //     for (int j=0; j<=i; j++)
        //       std::cout << " (" << alpha[j] << "," << coeffs[i][j] << ")";
        //     std::cout << std::endl;
        //   }

        // compute coefficients of gradient
        for (int s=0; s<d; s++)
          for (int i=0; i<n; i++)
            for (int j=0; j<=i; j++)
              (*gradcoeffs[s])[i][j] = 0;
        for (int i=0; i<n; i++)
          for (int j=0; j<=i; j++)
            for (int s=0; s<d; s++)
              if ((*alpha[j])[s]>0)
                {
                  MultiIndex<d> beta = *alpha[j]; // get exponents
                  FieldType factor = beta[s];
                  beta[s] -= 1;
                  int l = invert_index(beta);
                  (*gradcoeffs[s])[i][l] += (*coeffs)[i][j]*factor;
                }

        // for (int s=0; s<d; s++)
        //   {
        //     std::cout << "derivative in direction " << s << std::endl;
        //     for (int i=0; i<n; i++)
        //       {
        //         std::cout << "phi_" << i << ":" ;
        //         for (int j=0; j<=i; j++)
        //           std::cout << " (" << alpha[j] << "," << gradcoeffs[s][i][j] << ")";
        //         std::cout << std::endl;
        //       }
        //   }
      }

      // get index from a given multiindex
      int invert_index (MultiIndex<d>& a)
      {
        for (int i=0; i<n; i++)
          {
            bool found(true);
            for (int j=0; j<d; j++)
              if (a[j]!=(*alpha[i])[j]) found=false;
            if (found) return i;
          }
        DUNE_THROW(Dune::RangeError,"index not found in invertindex");
      }

      void gram_schmidt ()
      {
        // allocate a high precission matrix on the heap
        HighprecMat *p = new HighprecMat();
        HighprecMat& c = *p;

        // fill identity matrix
        for (int i=0; i<n; i++)
          for (int j=0; j<n; j++)
            if (i==j)
              c[i][j] = ComputationFieldType(1.0);
            else
              c[i][j] = ComputationFieldType(0.0);

        // the Gram-Schmidt loop
        MonomialIntegrator<ComputationFieldType,bt,d> integrator;
        for (int i=0; i<n; i++)
          {
            // store orthogonalization coefficients for scaling
            ComputationFieldType bi[n];

            // make p_i orthogonal to previous polynomials p_j
            for (int j=0; j<i; j++)
              {
                // p_j is represented with monomials
                bi[j] = ComputationFieldType(0.0);
                for (int l=0; l<=j; l++)
                  {
                    MultiIndex<d> a;
                    for (int m=0; m<d; m++) a[m] = (*alpha[i])[m] + (*alpha[l])[m];
                    bi[j] = bi[j] + c[j][l]*integrator.integrate(a);
                  }
                for (int l=0; l<=j; l++)
                  c[i][l] = c[i][l] - bi[j]*c[j][l];
              }

            // scale ith polynomial
            ComputationFieldType s2(0.0);
            MultiIndex<d> a;
            for (int m=0; m<d; m++) a[m] = (*alpha[i])[m] + (*alpha[i])[m];
            s2 = s2 + integrator.integrate(a);
            for (int j=0; j<i; j++)
              s2 = s2 - bi[j]*bi[j];
            ComputationFieldType s(1.0);
            using std::sqrt;
            s = s/sqrt(s2);
            for (int l=0; l<=i; l++)
              c[i][l] = s*c[i][l];
          }

        // store coefficients in low precission type
        for (int i=0; i<n; i++)
          for (int j=0; j<n; j++)
            (*coeffs)[i][j] = c[i][j];

        delete p;

        //std::cout << coeffs << std::endl;
      }
    };

  } // PB namespace

  // define the local finite element here

  template<class D, class R, int k, int d, Dune::GeometryType::BasicType bt, typename ComputationFieldType=Dune::PB::DefaultComputationalFieldType, PB::BasisType basisType = PB::BasisType::Pk>
  class OPBLocalBasis
  {
    typedef PB::BasisTraits<basisType> BasisTraits;
    typedef Dune::PB::OrthonormalPolynomialBasis<R,k,d,bt,ComputationFieldType,basisType> PolynomialBasis;
    PolynomialBasis opb;
    Dune::GeometryType gt;

  public:
    typedef Dune::LocalBasisTraits<D,d,Dune::FieldVector<D,d>,R,1,Dune::FieldVector<R,1>,Dune::FieldMatrix<R,1,d> > Traits;
    enum{ n = BasisTraits::template Size<k,d>::value };

    OPBLocalBasis (int order_) : opb(), gt(bt,d) {}

    template<class LFE>
    OPBLocalBasis (int order_, const LFE & lfe) : opb(lfe), gt(bt,d) {}

    unsigned int size () const { return n; }

    //! \brief Evaluate all shape functions
    inline void evaluateFunction (const typename Traits::DomainType& in,
                                  std::vector<typename Traits::RangeType>& out) const {
      out.resize(n);
      opb.evaluateFunction(in,out);
    }

    //! \brief Evaluate Jacobian of all shape functions
    inline void
    evaluateJacobian (const typename Traits::DomainType& in,
                      std::vector<typename Traits::JacobianType>& out) const {
      out.resize(n);
      opb.evaluateJacobian(in,out);
    }

    //! \brief Polynomial order of the shape functions
    unsigned int order () const  {
      return BasisTraits::template Order<k,d>::value;
    }

    Dune::GeometryType type () const { return gt; }
  };

  template<int k, int d, PB::BasisType basisType = PB::BasisType::Pk>
  class OPBLocalCoefficients
  {
    enum{ n = PB::BasisTraits<basisType>::template Size<k,d>::value };
  public:
    OPBLocalCoefficients (int order_) : li(n)  {
      for (int i=0; i<n; i++) li[i] = Dune::LocalKey(0,0,i);
    }

    //! number of coefficients
    std::size_t size () const { return n; }

    //! map index i to local key
    const Dune::LocalKey& localKey (int i) const  {
      return li[i];
    }

  private:
    std::vector<Dune::LocalKey> li;
  };

  template<class LB>
  class OPBLocalInterpolation
  {
    const LB& lb;

  public:
    OPBLocalInterpolation (const LB& lb_, int order_) : lb(lb_) {}

    //! \brief Local interpolation of a function
    template<typename F, typename C>
    void interpolate (const F& f, std::vector<C>& out) const
    {
      // select quadrature rule
      typedef typename FieldTraits<typename LB::Traits::RangeFieldType>::real_type RealFieldType;

      typedef typename LB::Traits::RangeType RangeType;
      const int d = LB::Traits::dimDomain;
      const Dune::QuadratureRule<RealFieldType,d>&
        rule = Dune::QuadratureRules<RealFieldType,d>::rule(lb.type(),2*lb.order());

      // prepare result
      out.resize(LB::n);
      for (int i=0; i<LB::n; i++) out[i] = 0.0;

      // loop over quadrature points
      for (typename Dune::QuadratureRule<RealFieldType,d>::const_iterator
             it=rule.begin(); it!=rule.end(); ++it)
        {
          // evaluate function at quadrature point
          typename LB::Traits::DomainType x;
          RangeType y;
          for (int i=0; i<d; i++) x[i] = it->position()[i];
          f.evaluate(x,y);

          // evaluate the basis
          std::vector<RangeType> phi(LB::n);
          lb.evaluateFunction(it->position(),phi);

          // do integration
          for (int i=0; i<LB::n; i++)
            out[i] += y*phi[i]*it->weight();
        }
    }
  };

  template<class D, class R, int k, int d, Dune::GeometryType::BasicType bt, typename ComputationFieldType=Dune::PB::DefaultComputationalFieldType, PB::BasisType basisType = PB::BasisType::Pk>
  class OPBLocalFiniteElement
  {
    Dune::GeometryType gt;
    OPBLocalBasis<D,R,k,d,bt,ComputationFieldType,basisType> basis;
    OPBLocalCoefficients<k,d,basisType> coefficients;
    OPBLocalInterpolation<OPBLocalBasis<D,R,k,d,bt,ComputationFieldType,basisType> > interpolation;
  public:
    typedef Dune::LocalFiniteElementTraits<OPBLocalBasis<D,R,k,d,bt,ComputationFieldType,basisType>,
                                           OPBLocalCoefficients<k,d,basisType>,
                                           OPBLocalInterpolation<OPBLocalBasis<D,R,k,d,bt,ComputationFieldType,basisType> > > Traits;

    OPBLocalFiniteElement ()
      : gt(bt,d), basis(k), coefficients(k), interpolation(basis,k)
    {}

    template<class LFE>
    explicit OPBLocalFiniteElement (const LFE & lfe)
      : gt(bt,d), basis(k, lfe), coefficients(k), interpolation(basis,k)
    {}

    OPBLocalFiniteElement (const OPBLocalFiniteElement & other)
      : gt(other.gt), basis(other.basis), coefficients(other.coefficients), interpolation(basis,k)
    {}

    const typename Traits::LocalBasisType& localBasis () const
    {
      return basis;
    }

    const typename Traits::LocalCoefficientsType& localCoefficients () const
    {
      return coefficients;
    }

    const typename Traits::LocalInterpolationType& localInterpolation () const
    {
      return interpolation;
    }

    Dune::GeometryType type () const { return gt; }

    OPBLocalFiniteElement* clone () const {
      return new OPBLocalFiniteElement(*this);
    }
  };

}

/** @} end documentation */

#endif // DUNE_PDELAB_FINITEELEMENT_L2ORTHONORMAL_HH
