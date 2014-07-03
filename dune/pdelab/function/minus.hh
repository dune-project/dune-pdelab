// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_FUNCTION_MINUS_HH
#define DUNE_PDELAB_FUNCTION_MINUS_HH

#include <cstddef>

#include <dune/pdelab/common/function.hh>

namespace Dune {
  namespace PDELab {

    //! Substract two GridFunctions
    /**
     * \tparam GF1 Type for the GridFunction on the left side of the minus
     * \tparam GF2 Type for the GridFunction on the right side of the minus
     */
    template<typename GF1, typename GF2>
    class MinusGridFunctionAdapter
      : public GridFunctionBase<typename GF1::Traits,
                                MinusGridFunctionAdapter<GF1,GF2> >
    {
      static_assert(std::size_t(GF1::Traits::dimRange) ==
                    std::size_t(GF2::Traits::dimRange),
                    "Range dimension must match for both operands of a "
                    "MinusGridFunctionAdapter");
      typedef typename GF1::Traits T;
      typedef GridFunctionBase<T, MinusGridFunctionAdapter<GF1,GF2> >
              Base;

      GF1& gf1;
      GF2& gf2;

    public:
      typedef typename Base::Traits Traits;

      MinusGridFunctionAdapter(GF1& gf1_, GF2& gf2_)
        : gf1(gf1_), gf2(gf2_)
      { }

      void evaluate(const typename Traits::ElementType &e,
                    const typename Traits::DomainType &x,
                    typename Traits::RangeType &y) const {
        gf1.evaluate(e,x,y);
        typename GF2::Traits::RangeType y2;
        gf2.evaluate(e,x,y2);
        y -= y2;
      }

      const typename Traits::GridViewType& getGridView() const {
        return gf1.getGridView();
      }

      template<typename Time>
      void setTime(Time time) {
        gf1.setTime(time);
        gf2.setTime(time);
      }
    };

  } // namspace PDELab
} // namspace Dune

#endif // DUNE_PDELAB_FUNCTION_MINUS_HH
