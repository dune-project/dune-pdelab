// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PkDG2DLOCALBASIS_HH
#define DUNE_PkDG2DLOCALBASIS_HH

#include <numeric>

#include <dune/common/fmatrix.hh>

#include <dune/localfunctions/common/localbasis.hh>

namespace Dune
{
  /**@ingroup LocalBasisImplementation
         \brief Lagrange shape functions of arbitrary order on the reference triangle.

         Lagrange shape functions of arbitrary order have the property that
         \f$\hat\phi^i(x_j) = \delta_{i,j}\f$ for certain points \f$x_j\f$.

         \tparam D Type to represent the field in the domain.
         \tparam R Type to represent the field in the range.
         \tparam k Polynomial order.

         \nosubgrouping
   */
  template<class D, class R, unsigned int k>
  class PkDG2DLocalBasis
  {
  public:

    /** \brief Export the number of degrees of freedom */
    enum {N = (k+1)*(k+2)/2};

    /** \brief Export the element order
       OS: Surprising that we need to export this both statically and dynamically!
     */
    enum {O = k};

    typedef LocalBasisTraits<D,2,Dune::FieldVector<D,2>,R,1,Dune::FieldVector<R,1>,
        Dune::FieldMatrix<R,1,2> > Traits;

    //! \brief Standard constructor
    PkDG2DLocalBasis ()
    {
      for (unsigned int i=0; i<=k; i++)
        pos_[i] = (1.0*i)/std::max(k,(unsigned int)1);
    }

    //! \brief number of shape functions
    unsigned int size () const
    {
      return N;
    }

    //! \brief Evaluate all shape functions
    inline void evaluateFunction (const typename Traits::DomainType& x,
                                  std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(N);
      // specialization for k==0, not clear whether that is needed
      if (k==0) {
        out[0] = 1;
        return;
      }

      int n=0;
      for (unsigned int j=0; j<=k; j++)
        for (unsigned int i=0; i<=k-j; i++)
        {
          out[n] = 1.0;
          for (unsigned int alpha=0; alpha<i; alpha++)
            out[n] *= (x[0]-pos_[alpha])/(pos_[i]-pos_[alpha]);
          for (unsigned int beta=0; beta<j; beta++)
            out[n] *= (x[1]-pos_[beta])/(pos_[j]-pos_[beta]);
          for (unsigned int gamma=i+j+1; gamma<=k; gamma++)
            out[n] *= (pos_[gamma]-x[0]-x[1])/(pos_[gamma]-pos_[i]-pos_[j]);
          n++;
        }
    }

    //! \brief Evaluate Jacobian of all shape functions
    inline void
    evaluateJacobian (const typename Traits::DomainType& x,       // position
                      std::vector<typename Traits::JacobianType>& out) const                        // return value
    {
      out.resize(N);

      // specialization for k==0, not clear whether that is needed
      if (k==0) {
        out[0][0][0] = 0; out[0][0][1] = 0;
        return;
      }

      int n=0;
      for (unsigned int j=0; j<=k; j++)
        for (unsigned int i=0; i<=k-j; i++)
        {
          // x_0 derivative
          out[n][0][0] = 0.0;
          R factor=1.0;
          for (unsigned int beta=0; beta<j; beta++)
            factor *= (x[1]-pos_[beta])/(pos_[j]-pos_[beta]);
          for (unsigned int a=0; a<i; a++)
          {
            R product=factor;
            for (unsigned int alpha=0; alpha<i; alpha++)
              if (alpha==a)
                product *= D(1)/(pos_[i]-pos_[alpha]);
              else
                product *= (x[0]-pos_[alpha])/(pos_[i]-pos_[alpha]);
            for (unsigned int gamma=i+j+1; gamma<=k; gamma++)
              product *= (pos_[gamma]-x[0]-x[1])/(pos_[gamma]-pos_[i]-pos_[j]);
            out[n][0][0] += product;
          }
          for (unsigned int c=i+j+1; c<=k; c++)
          {
            R product=factor;
            for (unsigned int alpha=0; alpha<i; alpha++)
              product *= (x[0]-pos_[alpha])/(pos_[i]-pos_[alpha]);
            for (unsigned int gamma=i+j+1; gamma<=k; gamma++)
              if (gamma==c)
                product *= -D(1)/(pos_[gamma]-pos_[i]-pos_[j]);
              else
                product *= (pos_[gamma]-x[0]-x[1])/(pos_[gamma]-pos_[i]-pos_[j]);
            out[n][0][0] += product;
          }

          // x_1 derivative
          out[n][0][1] = 0.0;
          factor = 1.0;
          for (unsigned int alpha=0; alpha<i; alpha++)
            factor *= (x[0]-pos_[alpha])/(pos_[i]-pos_[alpha]);
          for (unsigned int b=0; b<j; b++)
          {
            R product=factor;
            for (unsigned int beta=0; beta<j; beta++)
              if (beta==b)
                product *= D(1)/(pos_[j]-pos_[beta]);
              else
                product *= (x[1]-pos_[beta])/(pos_[j]-pos_[beta]);
            for (unsigned int gamma=i+j+1; gamma<=k; gamma++)
              product *= (pos_[gamma]-x[0]-x[1])/(pos_[gamma]-pos_[i]-pos_[j]);
            out[n][0][1] += product;
          }
          for (unsigned int c=i+j+1; c<=k; c++)
          {
            R product=factor;
            for (unsigned int beta=0; beta<j; beta++)
              product *= (x[1]-pos_[beta])/(pos_[j]-pos_[beta]);
            for (unsigned int gamma=i+j+1; gamma<=k; gamma++)
              if (gamma==c)
                product *= -D(1)/(pos_[gamma]-pos_[i]-pos_[j]);
              else
                product *= (pos_[gamma]-x[0]-x[1])/(pos_[gamma]-pos_[i]-pos_[j]);
            out[n][0][1] += product;
          }

          n++;
        }

    }

    /** \brief Evaluate partial derivatives of any order of all shape functions
     * \param order Order of the partial derivatives, in the classic multi-index notation
     * \param in Position where to evaluate the derivatives
     * \param[out] out Return value: the desired partial derivatives
     */
    void partial(const std::array<unsigned int,2>& order,
                 const typename Traits::DomainType& in,
                 std::vector<typename Traits::RangeType>& out) const
    {
      auto totalOrder = std::accumulate(order.begin(), order.end(), 0);

      switch (totalOrder)
      {
        case 0:
          evaluateFunction(in,out);
          break;
        case 1:
        {
          int direction = std::find(order.begin(), order.end(), 1)-order.begin();

          out.resize(N);

          int n=0;
          for (unsigned int j=0; j<=k; j++)
          {
            for (unsigned int i=0; i<=k-j; i++, n++)
            {
              out[n] = 0.0;
              for (unsigned int no1=0; no1 < k; no1++)
              {
                R factor = lagrangianFactorDerivative(direction, no1, i, j, in);
                for (unsigned int no2=0; no2 < k; no2++)
                  if (no1 != no2)
                    factor *= lagrangianFactor(no2, i, j, in);

                out[n] += factor;
              }
            }
          }

          break;
        }
        case 2:
        {
          out.resize(N);

          // specialization for k<2, not clear whether that is needed
          if (k<2)
          {
            std::fill(out.begin(), out.end(), 0.0);
            return;
          }

          std::array<int,2> directions;
          unsigned int counter = 0;
          auto nonconstOrder = order;  // need a copy that I can modify
          for (int i=0; i<2; i++)
          {
            while (nonconstOrder[i])
            {
              directions[counter++] = i;
              nonconstOrder[i]--;
            }
          }

          //f = prod_{i} f_i -> dxa dxb f = sum_{i} {dxa f_i sum_{k \neq i} dxb f_k prod_{l \neq k,i} f_l
          int n=0;
          for (unsigned int j=0; j<=k; j++)
          {
            for (unsigned int i=0; i<=k-j; i++, n++)
            {
              R res = 0.0;

              for (unsigned int no1=0; no1 < k; no1++)
              {
                R factor1 = lagrangianFactorDerivative(directions[0], no1, i, j, in);
                for (unsigned int no2=0; no2 < k; no2++)
                {
                  if (no1 == no2)
                    continue;
                  R factor2 = factor1*lagrangianFactorDerivative(directions[1], no2, i, j, in);
                  for (unsigned int no3=0; no3 < k; no3++)
                  {
                    if (no3 == no1 || no3 == no2)
                      continue;
                    factor2 *= lagrangianFactor(no3, i, j, in);
                  }
                  res += factor2;
                }
              }
              out[n] = res;
            }
          }

          break;
        }
        default:
          DUNE_THROW(NotImplemented, "Desired derivative order is not implemented");
      }
    }

    //! \brief Polynomial order of the shape functions
    unsigned int order () const
    {
      return k;
    }

  private:
  /** \brief Returns a single Lagrangian factor of l_ij evaluated at x */
  typename Traits::RangeType lagrangianFactor(const int no, const int i, const int j, const typename Traits::DomainType& x) const
  {
    if ( no < i)
      return (x[0]-pos_[no])/(pos_[i]-pos_[no]);
    if (no < i+j)
      return (x[1]-pos_[no-i])/(pos_[j]-pos_[no-i]);
    return (pos_[no+1]-x[0]-x[1])/(pos_[no+1]-pos_[i]-pos_[j]);
  }

  /** \brief Returns the derivative of a single Lagrangian factor of l_ij evaluated at x
   * \param direction Derive in x-direction if this is 0, otherwise derive in y direction
   */
  typename Traits::RangeType lagrangianFactorDerivative(const int direction, const int no, const int i, const int j, const typename Traits::DomainType& x) const
  {
    using T = typename Traits::RangeType;
    if ( no < i)
      return (direction == 0) ? T(1.0/(pos_[i]-pos_[no])) : T(0);

    if (no < i+j)
      return (direction == 0) ? T(0) : T(1.0/(pos_[j]-pos_[no-i]));

    return -1.0/(pos_[no+1]-pos_[i]-pos_[j]);
  }

    D pos_[k+1]; // positions on the interval
  };

}
#endif
