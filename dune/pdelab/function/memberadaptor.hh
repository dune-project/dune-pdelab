// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_FUNCTION_MEMBERADAPTOR_HH
#define DUNE_PDELAB_FUNCTION_MEMBERADAPTOR_HH

#include <cstddef>

#include <dune/common/fvector.hh>

#include <dune/pdelab/common/function.hh>

namespace Dune {
  namespace PDELab {

    //! GridFunction returning a constant value everywhere
    /**
     * \tparam GV       The type of the GridView
     * \tparam RF       The type of the range field
     * \tparam dimRange The dimension of the Range
     */
    template<class Member, class Class,
             class GV, class RF, std::size_t dimR = 1>
    class MemberFunctionToGridFunctionAdaptor :
      public GridFunctionBase<
        GridFunctionTraits< GV, RF, dimR, FieldVector<RF, dimR> >,
        MemberFunctionToGridFunctionAdaptor<Member, Class, GV, RF, dimR>
      >
    {
    public:
      typedef GridFunctionTraits< GV, RF, dimR, FieldVector<RF, dimR> > Traits;

    private:
      typedef GridFunctionBase<
        Traits,
        MemberFunctionToGridFunctionAdaptor<Member, Class, GV, RF, dimR>
        > Base;

      const Class &obj;
      Member member;
      const GV &gv;

    public:
      //! Construct a Const GridFunction
      /**
       * \param gv   The GridView to use.  It is passed as a reference to
       *             AnalyticGridFunctionBase (look there for the requirements
       *             of this argument).
       * \param val_ The value tu return on evaluation.  This class stores a
       *             copy of that value.
       */
      MemberFunctionToGridFunctionAdaptor(const Class &obj_, Member member_,
                                          const GV& gv_) :
        obj(obj_), member(member_), gv(gv_)
      { }

      void evaluate(const typename Traits::ElementType &e,
                    const typename Traits::DomainType &x,
                    typename Traits::RangeType &y) const
      {
        (obj.*member)(e, x, y);
      }

      const GV& getGridView() const { return gv; }
    };

    template<class RF, std::size_t dimRange,
             class GV, class Class, class Member>
    MemberFunctionToGridFunctionAdaptor<Member, Class, GV, RF, dimRange>
    makeMemberFunctionToGridFunctionAdaptor(const Class &obj, Member member,
                                            const GV &gv)
    {
      return MemberFunctionToGridFunctionAdaptor
        <Member, Class, GV, RF, dimRange>(obj, member, gv);
    }

  } // namspace PDELab
} // namspace Dune

#endif // DUNE_PDELAB_FUNCTION_MEMBERADAPTOR_HH
