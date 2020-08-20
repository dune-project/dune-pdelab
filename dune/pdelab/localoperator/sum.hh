// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_LOCALOPERATOR_SUM_HH
#define DUNE_PDELAB_LOCALOPERATOR_SUM_HH

#include <dune/pdelab/localoperator/combinedoperator.hh>

#include <utility>


namespace Dune {
  namespace PDELab {

    //! \addtogroup LocalOperator
    //! \ingroup PDELab
    //! \{

    //! A local operator to take the sum of other local operators
    /**
     * \nosubgrouping
     *
     * \tparam Args variadic list of local operators
     */
    template<typename... Args>
    class InstationarySumLocalOperator :
      public CombinedOperator<InstationarySumLocalOperator<Args...>, Args...>
    {
      template<typename F, typename... FArgs>
      void applyLops(F && f, FArgs &... args) const
      {
        Hybrid::forEach(std::make_index_sequence<sizeof...(Args)>{},
          [&](auto i){f(*Hybrid::elementAt(this->lops, i), args...);});
      }

      using Base = CombinedOperator<InstationarySumLocalOperator<Args...>, Args...>;
      using ArgPtrs = typename Base::ArgPtrs;

      friend Base;
    public:

      //! \brief Default-construct an InstationarySumLocalOperator. Expects the operators
      //!        to be added later through the setSummand method.
      InstationarySumLocalOperator()
      {}

      //! \brief construct a InstationarySumLocalOperator from a set of
      //!        local operators
      InstationarySumLocalOperator(Args&... lops)
        : Base(lops...)
      {}

      //! \brief construct a InstationarySumLocalOperator from a set of
      //!        local operators (rvalue reference)
      InstationarySumLocalOperator(Args&&... lops)
        : Base(std::forward<Args>(lops)...)
      {}

    protected:
      InstationarySumLocalOperator(ArgPtrs&& lops)
        : Base(std::forward<ArgPtrs>(lops))
      { }

    };

    //! \brief deprecated specialization of \see InstationarySumLocalOperator
    /**
     * \nosubgrouping
     *
     * \tparam Args Tuple of local operators.  Must fulfill \c
     *              tuple_size<Args>::value>=1.
     */
    template<typename... Args>
    class InstationarySumLocalOperator<std::tuple<Args...>> :
      public InstationarySumLocalOperator<Args...>
    {
      using Base = InstationarySumLocalOperator<Args...>;
      using ArgRefs = typename Base::ArgRefs;

    public:

      //! \brief Default-construct an InstationarySumLocalOperator. Expects the operators
      //!        to be added later through the setSummand method.
      [[deprecated("The specialization InstationarySumLocalOperator<Tuple<...>> is"
            "deprecated and will be removed after PDELab 2.7.")]]
      InstationarySumLocalOperator()
      {}

      //! \brief construct a InstationarySumLocalOperator from a tuple of
      //!        local operators
      [[deprecated("The specialization InstationarySumLocalOperator<Tuple<...>> is"
            "deprecated and will be removed after PDELab 2.7.")]]
      InstationarySumLocalOperator(const ArgRefs& lops)
        : Base(genericTransformTuple(lops,
            [](auto & l){return stackobject_to_shared_ptr(l);}))
      { }
    };

  }
}

#endif // DUNE_PDELAB_LOCALOPERATOR_SUM_HH
