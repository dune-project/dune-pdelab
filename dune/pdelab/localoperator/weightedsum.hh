// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_LOCALOPERATOR_WEIGHTEDSUM_HH
#define DUNE_PDELAB_LOCALOPERATOR_WEIGHTEDSUM_HH

#include <dune/common/concept.hh>
#include <dune/pdelab/localoperator/sum.hh>

namespace Dune {
  namespace PDELab {
    //! \addtogroup LocalOperator
    //! \ingroup PDELab
    //! \{

    //! A local operator to take the weighted sum of other local operators
    /**
     * \nosubgrouping
     *
     * If the weight for one summand is zero, calls to that local operators
     * evaluation and pattern methods are eliminated at run-time.
     *
     * \tparam K    Type of the scaling factors.
     * \tparam Args variadic list of local operators
     */
    template<typename K, typename... Args>
    class WeightedSumLocalOperator :
      public CombinedOperator<WeightedSumLocalOperator<K, Args...>, Args...>
    {
    protected:
      using Base = CombinedOperator<WeightedSumLocalOperator<K, Args...>, Args...>;
      using ArgPtrs = typename Base::ArgPtrs;
      using ArgRefs = typename Base::ArgRefs;
      friend Base;

      using Weights = FieldVector<K, sizeof...(Args)>;
      Weights weights;

      // concept check for a weighted container, like e.g. LocalVector or LocalMatrix
      struct WeightedContainer {
        template<class C>
        auto require(C& c) -> decltype(
          Concept::requireType<typename C::weight_type>(),
          const_cast<C&>(c).weight()
          // c.setWeight(std::declval<typename C::weight_type>())
          );
      };

      template<typename... FArgs>
      void getWeights(FieldVector<K, sizeof...(FArgs)> & aweights,
        std::tuple<FArgs...> fargs) const
      {
        Hybrid::forEach(std::make_index_sequence<sizeof...(FArgs)>{},
          [&](auto j){
            const auto & a = get<j>(fargs);
            Hybrid::ifElse(models<WeightedContainer,decltype(a)>(),
              [&](auto id){
                aweights[j] = id(a).weight();});
          });
      }

      template<typename... FArgs>
      void setWeights(const FieldVector<K, sizeof...(FArgs)> & aweights,
        std::tuple<FArgs...> fargs) const
      {
        Hybrid::forEach(std::make_index_sequence<sizeof...(FArgs)>{},
          [&](auto j){
            auto & a = get<j>(fargs);
            Hybrid::ifElse(models<WeightedContainer,decltype(a)>(),
              [&](auto id){
                id(a).setWeight(aweights[j]);});
          });
      }

      template<typename F, typename... FArgs>
      void applyLops(F && f, FArgs &... fargs) const
      {
        // remember weights
        FieldVector<K, sizeof...(FArgs)> aweights(K(0));
        FieldVector<K, sizeof...(FArgs)> current_weights;
        getWeights(aweights, std::forward_as_tuple(fargs...));
        Hybrid::forEach(std::make_index_sequence<sizeof...(Args)>{},
          [&](auto i){
            if(weights[i] != K(0)) {
              // set weights
              current_weights = aweights;
              current_weights *= weights[i];
              setWeights(current_weights, std::forward_as_tuple(fargs...));
              f(*Hybrid::elementAt(this->lops, i), fargs...);}}
          );
        // reset weights
        setWeights(aweights, std::forward_as_tuple(fargs...));
      }

    public:
      //////////////////////////////////////////////////////////////////////
      //
      //! \name Construction and modification
      //! \{
      //

      //! construct a WeightedSumLocalOperator
      /**
       * No summand local operators are set.  They must be initialized with
       * setSummand() before the constructed object is used.
       */
      WeightedSumLocalOperator (const Weights& weights_ = Weights(1))
        : weights(weights_)
      { }

      //! construct a WeightedSumLocalOperator from a set
      //! of local operators
      WeightedSumLocalOperator (Args&... lops_,
        const Weights& weights_ = Weights(1))
        : Base(lops_...), weights(weights_)
      { }

      //! construct a WeightedSumLocalOperator from a set
      //! of local operators (rvalue reference)
      WeightedSumLocalOperator (Args&&... lops_,
        const Weights& weights_ = Weights(1))
        : Base(std::forward<Args>(lops_)...), weights(weights_)
      { }

    protected:
      WeightedSumLocalOperator(ArgPtrs&& lops, const Weights& weights_)
        : Base(std::forward<ArgPtrs>(lops)), weights(weights_)
      { }

    public:
      //! set the weight for the i'th component of the sum
      void setWeight(K w, std::size_t i)
      { weights[i] = w; }

      //! get the weight for the i'th component of the sum
      K getWeight(std::size_t i)
      { return weights[i]; }

      //! \} Construction and modification
    };

    //! \brief deprecated specialization of \see WeightedSumLocalOperator
    /**
     * \nosubgrouping
     *
     * If the weight for one summand is zero, calls to that local operators
     * evaluation and pattern methods are eliminated at run-time.
     *
     * \tparam K    Type of the weighting factors.
     * \tparam Args Tuple of local operators.  Must fulfill \c
     *              tuple_size<Args>::value>=1.
     */
    template<typename K, typename... Args>
    class WeightedSumLocalOperator<K, std::tuple<Args...>> :
      public WeightedSumLocalOperator<K, Args...>
    {
      using Base = WeightedSumLocalOperator<K, Args...>;
      using ArgRefs = typename Base::ArgRefs;
      using Weights = typename Base::Weights;

    public:
      //! construct a WeightedSumLocalOperator
      /**
       * No summand local operators are set.  They must be initialized with
       * setSummand() before the constructed object is used.
       */
      [[deprecated("The specialization WeightedSumLocalOperator<K,Tuple<...>> is"
            "deprecated and will be removed after PDELab 2.7.")]]
      WeightedSumLocalOperator (const Weights& weights_ = Weights(1))
        : Base(weights_)
      { }

      //! construct a WeightedSumLocalOperator from a tuple of local operators
      [[deprecated("The specialization WeightedSumLocalOperator<K,Tuple<...>> is"
            "deprecated and will be removed after PDELab 2.7.")]]
      WeightedSumLocalOperator (const ArgRefs& lops_, const Weights& weights_ = Weights(1))
        : Base(genericTransformTuple(lops_,
            [](auto & l){return stackobject_to_shared_ptr(l);}), weights_)
      { }
    };

  }
}

#endif // DUNE_PDELAB_LOCALOPERATOR_WEIGHTEDSUM_HH
