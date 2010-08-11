// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_FUNCTION_SCALARSCALED_HH
#define DUNE_PDELAB_FUNCTION_SCALARSCALED_HH

#include <dune/pdelab/common/function.hh>

namespace Dune {
  namespace PDELab {

    //! Scale a GridFunction by a constant
    /**
     * \tparam GF The type of the GridFunction to scale
     */
    template<typename GF>
    class ScalarScaledGridFunctionAdapter
      : public GridFunctionBase<typename GF::Traits,
                                ScalarScaledGridFunctionAdapter<GF> >
    {
      typedef typename GF::Traits T;
      typedef GridFunctionBase<T, ScalarScaledGridFunctionAdapter<GF> >
              Base;
      typedef typename T::RangeFieldType RF;

      RF factor;
      GF& gf;

    public:
      typedef typename Base::Traits Traits;

      ScalarScaledGridFunctionAdapter(RF factor_, GF& gf_)
        : factor(factor_), gf(gf_)
      { }

      void evaluate(const typename Traits::ElementType &e,
                    const typename Traits::DomainType &x,
                    typename Traits::RangeType &y) const {
        gf.evaluate(e,x,y);
        y *= factor;
      }

      const typename Traits::GridViewType& getGridView() const {
        return gf.getGridView();
      }

      template<typename Time>
      void setTime(Time time) { gf.setTime(time); }
    };

  } // namspace PDELab
} // namspace Dune

#endif // DUNE_PDELAB_FUNCTION_SCALARSCALED_HH
