// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_BLOCKSTRUCTURED_QKLOCALBASIS_HH
#define DUNE_BLOCKSTRUCTURED_QKLOCALBASIS_HH

#include <numeric>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/power.hh>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>


namespace Dune
{
  namespace Blockstructured {
    /**@ingroup LocalBasisImplementation
       \brief Lagrange shape functions of order k on the reference cube.

       Also known as \f$Q^k\f$.

       \tparam D Type to represent the field in the domain.
       \tparam R Type to represent the field in the range.
       \tparam k Polynomial degree
       \tparam d Dimension of the cube

       \nosubgrouping
     */
    template<class D, class R, std::size_t k, std::size_t d, std::size_t blocks>
    class QkLocalBasis {

      static constexpr std::size_t DOFs1d = k * blocks + 1;

      // ith Lagrange polynomial of degree k in one dimension
      static R p(int i, D x) {
        R result(1.0);
        for (int j = 0; j <= DOFs1d - 1; j++)
          if (j != i) result *= ((DOFs1d - 1) * x - j) / (i - j);
        return result;
      }

      // derivative of ith Lagrange polynomial of degree k in one dimension
      static R dp(int i, D x) {
        R result(0.0);

        for (int j = 0; j <= DOFs1d - 1; j++)
          if (j != i) {
            R prod(((DOFs1d - 1) * 1.0) / (i - j));
            for (int l = 0; l <= DOFs1d - 1; l++)
              if (l != i && l != j)
                prod *= ((DOFs1d - 1) * x - l) / (i - l);
            result += prod;
          }
        return result;
      }

      // Second derivative of j-th Lagrange polynomial of degree k in one dimension
      // Formula and notation taken from https://en.wikipedia.org/wiki/Lagrange_polynomial#Derivatives
      static R ddp(int j, D x) {
        R result(0.0);

        for (int i = 0; i <= DOFs1d - 1; i++) {
          if (i == j)
            continue;

          R sum(0);

          for (int m = 0; m <= DOFs1d - 1; m++) {
            if (m == i || m == j)
              continue;

            R prod(((DOFs1d - 1) * 1.0) / (j - m));
            for (int l = 0; l <= DOFs1d - 1; l++)
              if (l != i && l != j && l != m)
                prod *= ((DOFs1d - 1) * x - l) / (j - l);
            sum += prod;
          }

          result += sum * (((DOFs1d - 1) * 1.0) / (j - i));
        }

        return result;
      }

      // Return i as a d-digit number in the DOFs1d -nary system
      static Dune::FieldVector<int, d> multiindex(int i) {
        Dune::FieldVector<int, d> alpha;
        for (int j = 0; j < d; j++) {
          alpha[j] = i % DOFs1d;
          i = i / DOFs1d;
        }
        return alpha;
      }

    public:
      typedef LocalBasisTraits<D, d, Dune::FieldVector<D, d>, R, 1, Dune::FieldVector<R, 1>, Dune::FieldMatrix<R, 1, d> > Traits;

      //! \brief number of shape functions
      unsigned int size() const {
        return StaticPower<DOFs1d, d>::power;
      }

      //! \brief Evaluate all shape functions
      inline void evaluateFunction(const typename Traits::DomainType &in,
                                   std::vector<typename Traits::RangeType> &out) const {
        out.resize(size());
        for (size_t i = 0; i < size(); i++) {
          // convert index i to multiindex
          Dune::FieldVector<int, d> alpha(multiindex(i));

          // initialize product
          out[i] = 1.0;

          // dimension by dimension
          for (int j = 0; j < d; j++)
            out[i] *= p(alpha[j], in[j]);
        }
      }

      /** \brief Evaluate Jacobian of all shape functions
       * \param in position where to evaluate
       * \param out The return value
       */
      inline void
      evaluateJacobian(const typename Traits::DomainType &in,
                       std::vector<typename Traits::JacobianType> &out) const {
        out.resize(size());

        // Loop over all shape functions
        for (size_t i = 0; i < size(); i++) {
          // convert index i to multiindex
          Dune::FieldVector<int, d> alpha(multiindex(i));

          // Loop over all coordinate directions
          for (int j = 0; j < d; j++) {
            // Initialize: the overall expression is a product
            // if j-th bit of i is set to -1, else 1
            out[i][0][j] = dp(alpha[j], in[j]);

            // rest of the product
            for (int l = 0; l < d; l++)
              if (l != j)
                out[i][0][j] *= p(alpha[l], in[l]);
          }
        }
      }

      /** \brief Evaluate partial derivatives of any order of all shape functions
       * \param order Order of the partial derivatives, in the classic multi-index notation
       * \param in Position where to evaluate the derivatives
       * \param[out] out Return value: the desired partial derivatives
       */
      inline void partial(const std::array<unsigned int, d> &order,
                          const typename Traits::DomainType &in,
                          std::vector<typename Traits::RangeType> &out) const {
        out.resize(size());

        // Loop over all shape functions
        for (size_t i = 0; i < size(); i++) {
          // convert index i to multiindex
          Dune::FieldVector<int, d> alpha(multiindex(i));

          // Initialize: the overall expression is a product
          out[i][0] = 1.0;

          // rest of the product
          for (std::size_t l = 0; l < d; l++) {
            switch (order[l]) {
              case 0:
                out[i][0] *= p(alpha[l], in[l]);
                break;
              case 1:
                out[i][0] *= dp(alpha[l], in[l]);
                break;
              case 2:
                out[i][0] *= ddp(alpha[l], in[l]);
                break;
              default:
                DUNE_THROW(NotImplemented, "Desired derivative order is not implemented");
            }
          }
        }
      }

      //! \brief Polynomial order of the shape functions
      unsigned int order() const {
        return k;
      }
    };
  }
}

#endif
