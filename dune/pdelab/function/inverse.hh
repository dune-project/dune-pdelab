// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_FUNCTION_INVERSE_HH
#define DUNE_PDELAB_FUNCTION_INVERSE_HH

#include <dune/pdelab/common/function.hh>

namespace Dune {
  namespace PDELab {

    //! Take inverse of a GridFunction
    /**
     * \tparam GF The type of the GridFunction to invert
     */
    template<typename GF>
    class InverseGridFunctionAdapter
      : public GridFunctionBase<typename GF::Traits,
                                InverseGridFunctionAdapter<GF> >
    {
      static_assert(GF::Traits::dimRange == 1, "Dimension of range must "
                    "be 1 to take the inverse");

      typedef typename GF::Traits T;
      typedef GridFunctionBase<T, InverseGridFunctionAdapter<GF> >
              Base;
      typedef typename T::RangeFieldType RF;

      GF& gf;

    public:
      typedef typename Base::Traits Traits;

      InverseGridFunctionAdapter(GF& gf_)
        : gf(gf_)
      { }

      void evaluate(const typename Traits::ElementType &e,
                    const typename Traits::DomainType &x,
                    typename Traits::RangeType &y) const {
        gf.evaluate(e,x,y);
        y = 1/y;
      }

      const typename Traits::GridViewType& getGridView() const {
        return gf.getGridView();
      }

      template<typename Time>
      void setTime(Time time) { gf.setTime(time); }
    };

  } // namspace PDELab
} // namspace Dune

#endif // DUNE_PDELAB_FUNCTION_INVERSE_HH
