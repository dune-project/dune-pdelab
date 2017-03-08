// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_GRIDFUNCTIONADAPTER_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_GRIDFUNCTIONADAPTER_HH

#include <dune/pdelab/common/function.hh>

namespace Dune {
  namespace PDELab {

    /*! \brief Adapter returning f1(x)-f2(x) for two given grid functions

      \tparam T1  a grid function type
      \tparam T2  a grid function type
    */
    template<typename T1, typename T2>
    class DifferenceAdapter
      : public Dune::PDELab::GridFunctionBase<
      Dune::PDELab::GridFunctionTraits<typename T1::Traits::GridViewType,
                                       typename T1::Traits::RangeFieldType,
                                       1,
                                       Dune::FieldVector<typename T1::Traits::RangeFieldType,1> >,
      DifferenceAdapter<T1,T2> >
    {
    public:
      typedef Dune::PDELab::GridFunctionTraits<typename T1::Traits::GridViewType,
                                               typename T1::Traits::RangeFieldType,
                                               1,Dune::FieldVector<typename T1::Traits::RangeFieldType,1> > Traits;

      //! constructor
      DifferenceAdapter (const T1& t1_, const T2& t2_) : t1(t1_), t2(t2_) {}

      //! \copydoc GridFunctionBase::evaluate()
      inline void evaluate (const typename Traits::ElementType& e,
                            const typename Traits::DomainType& x,
                            typename Traits::RangeType& y) const
      {
        typename Traits::RangeType y1;
        t1.evaluate(e,x,y1);
        typename Traits::RangeType y2;
        t2.evaluate(e,x,y2);
        y1 -= y2;
        y = y1;
      }

      inline const typename Traits::GridViewType& getGridView () const
      {
        return t1.getGridView();
      }

    private:
      const T1& t1;
      const T2& t2;
    };


    /*! \brief Adapter returning ||f1(x)-f2(x)||^2 for two given grid functions

      \tparam T1  a grid function type
      \tparam T2  a grid function type
    */
    template<typename T1, typename T2>
    class DifferenceSquaredAdapter
      : public Dune::PDELab::GridFunctionBase<
      Dune::PDELab::GridFunctionTraits<typename T1::Traits::GridViewType,
                                       typename T1::Traits::RangeFieldType,
                                       1,
                                       Dune::FieldVector<typename T1::Traits::RangeFieldType,1> >,
      DifferenceSquaredAdapter<T1,T2> >
    {
    public:
      typedef Dune::PDELab::GridFunctionTraits<typename T1::Traits::GridViewType,
                                               typename T1::Traits::RangeFieldType,
                                               1,Dune::FieldVector<typename T1::Traits::RangeFieldType,1> > Traits;

      //! constructor
      DifferenceSquaredAdapter (const T1& t1_, const T2& t2_) : t1(t1_), t2(t2_) {}

      //! \copydoc GridFunctionBase::evaluate()
      inline void evaluate (const typename Traits::ElementType& e,
                            const typename Traits::DomainType& x,
                            typename Traits::RangeType& y) const
      {
        typename T1::Traits::RangeType y1;
        t1.evaluate(e,x,y1);
        typename T2::Traits::RangeType y2;
        t2.evaluate(e,x,y2);
        y1 -= y2;
        y = y1.two_norm2();
      }

      inline const typename Traits::GridViewType& getGridView () const
      {
        return t1.getGridView();
      }

    private:
      const T1& t1;
      const T2& t2;
    };

  }
}
#endif
