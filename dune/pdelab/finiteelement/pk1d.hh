// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// Pk in one dimension with k as runtime variable

#ifndef DUNE_PDELAB_FINITEELEMENT_PK1D_HH
#define DUNE_PDELAB_FINITEELEMENT_PK1D_HH

#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/geometry/type.hh>

#include<dune/localfunctions/common/localbasis.hh>
#include<dune/localfunctions/common/localkey.hh>
#include<dune/localfunctions/common/localfiniteelementtraits.hh>

namespace Dune {

  /** \brief Define the Pk Lagrange basis functions in 1d on the reference interval
   *
   * Unlike the corresponding implementation in dune-localfunctions, the order k
   * is a run-time parameter here!
   *
   *  \tparam D Type to represent domain.
   *  \tparam R Type to represent range.
   */
  template<class D, class R>
  class Pk1dLocalFiniteElement
  {
    //! \brief Class for the basis functions
    class Pk1dLocalBasis
    {
      Dune::GeometryType gt; // store geometry type for the basis
      std::size_t k;         // polynomial degree
      std::size_t n;         // the number of basis functions
      std::vector<R> s;      // Lagrange points on the reference interval
    public:
      typedef Dune::LocalBasisTraits<D,1,Dune::FieldVector<D,1>,R,1,Dune::FieldVector<R,1>,Dune::FieldMatrix<R,1,1>, 1> Traits;

      //! \brief make a basis object for given polynomial degree
      Pk1dLocalBasis (std::size_t k_) : gt(Dune::GeometryType::cube,1), k(k_), n(k_+1), s(n)
      {
        for (std::size_t i=0; i<=k; i++) s[i] = (1.0*i)/k;
      }

      //! \brief return number of basis functions
      std::size_t size () const { return n; }

      //! \brief Evaluate all shape functions at a given point in local coordinates
      inline void evaluateFunction (const typename Traits::DomainType& in,
                                    std::vector<typename Traits::RangeType>& out) const {
        out.resize(n);
        for (std::size_t i=0; i<=k; i++)
          {
            out[i] = 1.0;
            for (std::size_t j=0; j<=k; j++)
              if (i!=j) out[i] *= (in[0]-s[j])/(s[i]-s[j]);
          }
      }

      //! \brief Evaluate Jacobian of all shape functions
      inline void
      evaluateJacobian (const typename Traits::DomainType& in,
                        std::vector<typename Traits::JacobianType>& out) const {
        out.resize(n);
        for (std::size_t i=0; i<=k; i++) // derivative of basis function i
          {
            out[i][0][0] = 0.0;
            R factor = 1.0;
            R denominator = 1.0;
            for (std::size_t j=0; j<=k; j++)
              {
                if (j==i) continue; // treat factor (x-s_j)
                denominator *= s[i]-s[j];
                R a=1.0;            // product of remaining factors (might be empty)
                for (std::size_t l=j+1; l<=k; l++)
                  {
                    if (l==i) continue;
                    a *= in[0]-s[l];
                  }
                out[i][0][0] += factor*a;
                factor *= in[0]-s[j];
              }
            out[i][0][0] /= denominator;
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

      //! \brief Polynomial order of the basis functions
      unsigned int order () const  {
        return k;
      }

      //! \brief return geometry type
      Dune::GeometryType type () const { return gt; }
    };

    //! \brief Class for the basis functions
    class Pk1dLocalCoefficients
    {
    public:
      Pk1dLocalCoefficients (std::size_t k_) : k(k_), n(k_+1), li(k_+1)  {
        li[0] = Dune::LocalKey(0,1,0);
        for (int i=1; i<int(k); i++) li[i] = Dune::LocalKey(0,0,i-1);
        li[k] = Dune::LocalKey(1,1,0);
      }

      //! number of coefficients
      std::size_t size () const { return n; }

      //! map index i to local key
      const Dune::LocalKey& localKey (int i) const  {
        return li[i];
      }

    private:
      std::size_t k;                  // polynomial degree
      std::size_t n;                  // the number of basis functions
      std::vector<Dune::LocalKey> li; // assignment of basis function to subentities
    };

    //! \brief Class for interpolating a given function by the basis
    template<typename LB>
    class Pk1dLocalInterpolation
    {
    public:
      Pk1dLocalInterpolation (std::size_t k_) : k(k_), n(k_+1) {}

      //! \brief Local interpolation of a function
      template<typename F, typename C>
      void interpolate (const F& f, std::vector<C>& out) const
      {
        out.resize(n);
        typename LB::Traits::DomainType x;
        typename LB::Traits::RangeType y;
        for (int i=0; i<=int(k); i++)
          {
            x[0] = (1.0*i)/k; // the point to evaluate
            f.evaluate(x,y);
            out[i] = y[0];
          }
      }
    private:
      std::size_t k;                  // polynomial degree
      std::size_t n;                  // the number of basis functions
    };

    Dune::GeometryType gt;
    Pk1dLocalBasis basis;
    Pk1dLocalCoefficients coefficients;
    Pk1dLocalInterpolation<Pk1dLocalBasis> interpolation;

  public:
    typedef Dune::LocalFiniteElementTraits<Pk1dLocalBasis,
                                           Pk1dLocalCoefficients,
                                           Pk1dLocalInterpolation<Pk1dLocalBasis> > Traits;

    Pk1dLocalFiniteElement (std::size_t k)
      : gt(Dune::GeometryType::cube,1), basis(k), coefficients(k), interpolation(k)
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

    /** \todo Please doc me !
     */
    std::size_t size() const
    {
      return basis.size();
    }

    Dune::GeometryType type () const { return gt; }

    Pk1dLocalFiniteElement* clone () const {
      return new Pk1dLocalFiniteElement(*this);
    }
  };
}

#endif // DUNE_PDELAB_FINITEELEMENT_PK1D_HH
