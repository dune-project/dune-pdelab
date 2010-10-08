// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_FUNCTION_MEMBERADAPTOR_HH
#define DUNE_PDELAB_FUNCTION_MEMBERADAPTOR_HH

#include <cstddef>

#include <dune/common/fvector.hh>

#include <dune/pdelab/common/function.hh>

namespace Dune {
  namespace PDELab {

    //! GridFunction implemented by a member function of some class
    /**
     * \tparam Member   Member function pointer type.
     * \tparam Class    Type of the class containing the member.
     * \tparam GV       The type of the GridView
     * \tparam RF       The type of the range field
     * \tparam dimRange The dimension of the Range.
     *
     * The member function must support the signature
     * \code
(obj.*member)(const FieldVector<typename GV::ctype, GV::dimension> &xl,
              FieldVector<RF, dimRange> &y) const;
     * \endcode
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
      //! export traits class
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
      //! Construct an adaptor object
      /**
       * \param obj_    Class object to call the member function on.
       * \param member_ Member function pointer to the member to use.
       * \param gv_     The GridView to use.
       *
       * This class store the \c obj_ and \c gv_ references internally and
       * becomes invalid when they become invalid.
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

      //! get reference to the internal gridview.
      const GV& getGridView() const { return gv; }
    };

    //! easy construction of a MemberFunctionToGridFunctionAdaptor
    /**
     * \relates MemberFunctionToGridFunctionAdaptor
     *
     * \code
#include <dune/pdelab/function/memberadaptor.hh>
     * \endcode
     *
     * The only required template parameters are \c RF and \c dimRange, all
     * other template parameters can be deduced automatically.
     */
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
