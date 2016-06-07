// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_FUNCTION_PRODUCT_HH
#define DUNE_PDELAB_FUNCTION_PRODUCT_HH

// #include <cstddef>

#include <dune/common/fvector.hh>
#include <dune/common/typetraits.hh>

#include <dune/pdelab/common/function.hh>

namespace Dune {
  namespace PDELab {

    //! Product of two GridFunctions
    template<typename GF1, typename GF2, class = void>
    class ProductGridFunctionAdapter :
      public GridFunctionBase<
        GridFunctionTraits<
          typename GF1::Traits::GridViewType,
          typename GF1::Traits::RangeFieldType, 1,
          FieldVector<typename GF1::Traits::RangeFieldType, 1> >,
        ProductGridFunctionAdapter<GF1,GF2> >
    {
      static_assert(unsigned(GF1::Traits::dimRange) ==
                    unsigned(GF2::Traits::dimRange),
                    "ProductGridFunctionAdapter: Operands must have "
                    "matching range dimensions, or one operand must be "
                    "scalar-valued.");

      typedef GridFunctionTraits<
        typename GF1::Traits::GridViewType,
        typename GF1::Traits::RangeFieldType, 1,
        FieldVector<typename GF1::Traits::RangeFieldType, 1> > T;
      typedef GridFunctionBase<T, ProductGridFunctionAdapter<GF1,GF2> > Base;

      GF1& gf1;
      GF2& gf2;

    public:
      typedef typename Base::Traits Traits;

      ProductGridFunctionAdapter(GF1& gf1_, GF2& gf2_)
        : gf1(gf1_), gf2(gf2_)
      { }

      void evaluate(const typename Traits::ElementType &e,
                    const typename Traits::DomainType &x,
                    typename Traits::RangeType &y) const {
        typename GF1::Traits::RangeType y1;
        gf1.evaluate(e,x,y1);
        typename GF2::Traits::RangeType y2;
        gf2.evaluate(e,x,y2);
        y = y1 * y2;
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

    //! Product of two GridFunctions
    template<typename GF1, typename GF2>
    class ProductGridFunctionAdapter<
      GF1, GF2,
      typename std::enable_if<
        GF1::Traits::dimRange == 1 && GF2::Traits::dimRange != 1
        >::type> :
      public GridFunctionBase<typename GF2::Traits,
                              ProductGridFunctionAdapter<GF1,GF2> >
    {
      typedef typename GF2::Traits T;
      typedef GridFunctionBase<T, ProductGridFunctionAdapter<GF1,GF2> > Base;

      GF1& gf1;
      GF2& gf2;

    public:
      typedef typename Base::Traits Traits;

      ProductGridFunctionAdapter(GF1& gf1_, GF2& gf2_)
        : gf1(gf1_), gf2(gf2_)
      { }

      void evaluate(const typename Traits::ElementType &e,
                    const typename Traits::DomainType &x,
                    typename Traits::RangeType &y) const {
        typename GF1::Traits::RangeType y1;
        gf1.evaluate(e,x,y1);
        gf2.evaluate(e,x,y);
        y *= y1;
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

    //! Product of two GridFunctions
    template<typename GF1, typename GF2>
    class ProductGridFunctionAdapter<
      GF1, GF2,
      typename std::enable_if<
        GF1::Traits::dimRange != 1 && GF2::Traits::dimRange == 1
        >::type> :
      public ProductGridFunctionAdapter<GF2, GF1>
    {
    public:
      ProductGridFunctionAdapter(GF1& gf1, GF2& gf2)
        : ProductGridFunctionAdapter<GF2, GF1>(gf2, gf1)
      { }
    };

  } // namspace PDELab
} // namspace Dune

#endif // DUNE_PDELAB_FUNCTION_PRODUCT_HH
