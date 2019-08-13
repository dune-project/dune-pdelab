// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PkDG1DLOCALBASIS_HH
#define DUNE_PkDG1DLOCALBASIS_HH

#include <dune/common/fmatrix.hh>

#include <dune/localfunctions/common/localbasis.hh>

namespace Dune
{
  /**@ingroup LocalBasisImplementation
     \brief Lagrange shape functions of arbitrary order on the 1D reference triangle.

     Lagrange shape functions of arbitrary order have the property that
     \f$\hat\phi^i(x_j) = \delta_{i,j}\f$ for certain points \f$x_j\f$.

     \tparam D Type to represent the field in the domain.
     \tparam R Type to represent the field in the range.
     \tparam k Polynomial order.

     \nosubgrouping
   */
  template<class D, class R, unsigned int k>
  class PkDG1DLocalBasis
  {
  public:

    /** \brief Export the number of degrees of freedom */
    enum {N = k+1};

    /** \brief Export the element order */
    enum {O = k};

    typedef LocalBasisTraits<D,
        1,
        Dune::FieldVector<D,1>,
        R,
        1,
        Dune::FieldVector<R,1>,
        Dune::FieldMatrix<R,1,1>
        > Traits;

    //! \brief Standard constructor
    PkDG1DLocalBasis ()
    {
      for (unsigned int i=0; i<=k; i++)
        pos[i] = (1.0*i)/std::max(k,(unsigned int)1);
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

      for (unsigned int i=0; i<N; i++)
      {
        out[i] = 1.0;
        for (unsigned int alpha=0; alpha<i; alpha++)
          out[i] *= (x[0]-pos[alpha])/(pos[i]-pos[alpha]);
        for (unsigned int gamma=i+1; gamma<=k; gamma++)
          out[i] *= (x[0]-pos[gamma])/(pos[i]-pos[gamma]);
      }
    }

    //! \brief Evaluate Jacobian of all shape functions
    inline void
    evaluateJacobian (const typename Traits::DomainType& x,             // position
                      std::vector<typename Traits::JacobianType>& out) const          // return value
    {
      out.resize(N);

      for (unsigned int i=0; i<=k; i++) {

        // x_0 derivative
        out[i][0][0] = 0.0;
        R factor=1.0;
        for (unsigned int a=0; a<i; a++)
        {
          R product=factor;
          for (unsigned int alpha=0; alpha<i; alpha++)
            product *= (alpha==a) ? 1.0/(pos[i]-pos[alpha])
                       : (x[0]-pos[alpha])/(pos[i]-pos[alpha]);
          for (unsigned int gamma=i+1; gamma<=k; gamma++)
            product *= (pos[gamma]-x[0])/(pos[gamma]-pos[i]);
          out[i][0][0] += product;
        }
        for (unsigned int c=i+1; c<=k; c++)
        {
          R product=factor;
          for (unsigned int alpha=0; alpha<i; alpha++)
            product *= (x[0]-pos[alpha])/(pos[i]-pos[alpha]);
          for (unsigned int gamma=i+1; gamma<=k; gamma++)
            product *= (gamma==c) ? -1.0/(pos[gamma]-pos[i])
                       : (pos[gamma]-x[0])/(pos[gamma]-pos[i]);
          out[i][0][0] += product;
        }
      }

    }

    /** \brief Evaluate partial derivatives of any order of all shape functions
     * \param order Order of the partial derivatives, in the classic multi-index notation
     * \param in Position where to evaluate the derivatives
     * \param[out] out Return value: the desired partial derivatives
     */
    void partial(const std::array<unsigned int,1>& order,
                 const typename Traits::DomainType& in,
                 std::vector<typename Traits::RangeType>& out) const
    {
      switch (order[0])
      {
        case 0:
          evaluateFunction(in,out);
          break;
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
    R pos[k+1];     // positions on the interval
  };

}
#endif
