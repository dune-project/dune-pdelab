#ifndef DUNE_PDELAB_OPERATOR_FORWARD_INSTATIONARY_TRAITS_HH
#define DUNE_PDELAB_OPERATOR_FORWARD_INSTATIONARY_TRAITS_HH


namespace Dune::PDELab::inline Experimental {


enum class DurationPosition {StiffnessNumerator, MassDenominator};

template<DurationPosition dt_position = DurationPosition::StiffnessNumerator>
class InstationaryTraits {

public:

  template<class Duration>
  using MassFactor = std::conditional_t<(dt_position == DurationPosition::MassDenominator), decltype(double{}/Duration{}), double>;

  template<class Duration>
  using StiffnessFactor = std::conditional_t<(dt_position == DurationPosition::StiffnessNumerator), decltype(double{}*Duration{}), double>;

  template<class Duration>
  static constexpr MassFactor<Duration> massFactor(Duration duration)
  {
    if constexpr (dt_position == DurationPosition::MassDenominator)
      return 1. / duration;
    else
      return 1.;
  }

  template<class Duration>
  static constexpr StiffnessFactor<Duration> stiffnessFactor(Duration duration)
  {
    if constexpr (dt_position == DurationPosition::StiffnessNumerator)
      return 1. * duration;
    else
      return 1.;
  }

};


} // namespace Dune::PDELab::inline Experimental

#endif // DUNE_PDELAB_OPERATOR_FORWARD_INSTATIONARY_TRAITS_HH
