#ifndef DUNE_PDELAB_OPERATOR_FORWARD_INSTATIONARY_TRAITS_HH
#define DUNE_PDELAB_OPERATOR_FORWARD_INSTATIONARY_TRAITS_HH


namespace Dune::PDELab::inline Experimental {


enum class DurationPosition {StiffnessNumerator, MassDenominator};

template<DurationPosition dt_position = DurationPosition::StiffnessNumerator>
class InstationaryTraits {

public:

  template<class DurationQuantity>
  using MassFactor = std::conditional_t<(dt_position == DurationPosition::MassDenominator), decltype(double{}/DurationQuantity{}), double>;

  template<class DurationQuantity>
  using StiffnessFactor = std::conditional_t<(dt_position == DurationPosition::StiffnessNumerator), decltype(double{}*DurationQuantity{}), double>;

  template<class DurationQuantity>
  static constexpr MassFactor<DurationQuantity> massFactor(DurationQuantity duration)
  {
    if constexpr (dt_position == DurationPosition::MassDenominator)
      return 1. / duration;
    else
      return 1.;
  }

  template<class DurationQuantity>
  static constexpr StiffnessFactor<DurationQuantity> stiffnessFactor(DurationQuantity duration)
  {
    if constexpr (dt_position == DurationPosition::StiffnessNumerator)
      return 1. * duration;
    else
      return 1.;
  }

};


} // namespace Dune::PDELab::inline Experimental

#endif // DUNE_PDELAB_OPERATOR_FORWARD_INSTATIONARY_TRAITS_HH
