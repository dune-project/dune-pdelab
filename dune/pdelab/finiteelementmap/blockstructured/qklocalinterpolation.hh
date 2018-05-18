// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_BLOCKSTRUCTURED_QKLOCALINTERPOLATION_HH
#define DUNE_BLOCKSTRUCTURED_QKLOCALINTERPOLATION_HH

#include <dune/common/fvector.hh>
#include <dune/common/power.hh>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localinterpolation.hh>


namespace Dune
{
  namespace Blockstructured {
    /** \todo Please doc me! */
    template<std::size_t k, std::size_t d, std::size_t blocks, class LB>
    class QkLocalInterpolation {

      static constexpr std::size_t DOFs1d = k * blocks + 1;

      // Return i as a d-digit number in the (k+1)-nary system
      static Dune::FieldVector<int, d> multiindex(int i) {
        Dune::FieldVector<int, d> alpha;
        for (int j = 0; j < d; j++) {
          alpha[j] = i % DOFs1d;
          i = i / DOFs1d;
        }
        return alpha;
      }

    public:

      //! \brief Local interpolation of a function
      template<typename F, typename C>
      void interpolate(const F &ff, std::vector<C> &out) const {
        typename LB::Traits::DomainType x;

        auto &&f = Impl::makeFunctionWithCallOperator<typename LB::Traits::DomainType>(ff);

        out.resize(StaticPower<DOFs1d, d>::power);

        for (int i = 0; i < StaticPower<DOFs1d, d>::power; i++) {
          // convert index i to multiindex
          Dune::FieldVector<int, d> alpha(multiindex(i));

          // Generate coordinate of the i-th Lagrange point
          for (int j = 0; j < d; j++)
            x[j] = (1.0 * alpha[j]) / (DOFs1d - 1);

          out[i] = f(x);
        }
      }
    };

    /** \todo Please doc me! */
    template<std::size_t d, class LB>
    class QkLocalInterpolation<0, d, 1, LB> {
    public:
      //! \brief Local interpolation of a function
      template<typename F, typename C>
      void interpolate(const F &ff, std::vector<C> &out) const {
        typename LB::Traits::DomainType x(0);

        auto &&f = Impl::makeFunctionWithCallOperator<typename LB::Traits::DomainType>(ff);

        out.resize(1);
        out[0] = f(x);
      }
    };
  }
}


#endif
