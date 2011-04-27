// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_FUNCTION_SQR_HH
#define DUNE_PDELAB_FUNCTION_SQR_HH

#include <dune/common/fvector.hh>

#include <dune/pdelab/common/function.hh>

namespace Dune {
  namespace PDELab {

    //! Take square of a GridFunction
    /**
     * \tparam GF The type of the GridFunction to square
     */
    template<typename GF>
    class SqrGridFunctionAdapter :
      public GridFunctionBase<
        GridFunctionTraits<
          typename GF::Traits::GridViewType,
          typename GF::Traits::RangeFieldType, 1,
          FieldVector<typename GF::Traits::RangeFieldType, 1>
          >,
        SqrGridFunctionAdapter<GF>
      >
    {
      typedef GridFunctionTraits<
        typename GF::Traits::GridViewType,
        typename GF::Traits::RangeFieldType, 1,
        FieldVector<typename GF::Traits::RangeFieldType, 1>
        > T;
      typedef GridFunctionBase<T, SqrGridFunctionAdapter<GF> >
              Base;
      typedef typename T::RangeFieldType RF;

      GF& gf;

    public:
      typedef typename Base::Traits Traits;

      SqrGridFunctionAdapter(GF& gf_)
        : gf(gf_)
      { }

      void evaluate(const typename Traits::ElementType &e,
                    const typename Traits::DomainType &x,
                    typename Traits::RangeType &y) const
      {
        typename GF::Traits::RangeType y_;
        gf.evaluate(e,x,y_);
        y[0] = y_*y_;
      }

      const typename Traits::GridViewType& getGridView() const {
        return gf.getGridView();
      }

      template<typename Time>
      void setTime(Time time) { gf.setTime(time); }
    };

  } // namspace PDELab
} // namspace Dune

#endif // DUNE_PDELAB_FUNCTION_SQR_HH
