#ifndef DUNE_PDELAB_SOLVER_UTILITY_HH
#define DUNE_PDELAB_SOLVER_UTILITY_HH

#include <chrono>
#include <cstdlib>
#include <string_view>

#include <dune/common/exceptions.hh>

#include <dune/logging/fmt.hh>

namespace Dune::PDELab {

  enum class LineSearchStrategy
  {
    none,
    hackbuschReusken,
    hackbuschReuskenAcceptBest
  };

  constexpr std::string_view name(LineSearchStrategy strategy) noexcept
  {
    using namespace std::literals;
    switch (strategy)
    {
    case LineSearchStrategy::none:
      return "none"sv;
    case LineSearchStrategy::hackbuschReusken:
      return "hackbuschReusken"sv;
    case LineSearchStrategy::hackbuschReuskenAcceptBest:
      return "hackbuschReuskenAcceptBest"sv;
    default:
      std::abort();
    }
  }

  LineSearchStrategy lineSearchStrategyFromString(const std::string& name)
  {
    if (name == "none")
      return LineSearchStrategy::none;
    if (name == "hackbusch_reusken")
      return LineSearchStrategy::hackbuschReusken;
    if (name == "hackbusch_reusken_accept_best")
      return LineSearchStrategy::hackbuschReuskenAcceptBest;
    DUNE_THROW(Exception,"Unkown line search strategy: " << name);
  }

  template<typename Clock>
  class Timing
  {

  public:

    using Duration = typename Clock::duration;

    template<typename F>
    Duration operator()(F f)
    {
      auto start = Clock::now();
      std::invoke(f);
      auto end = Clock::now();
      auto elapsed = end - start;
      _total += elapsed;
      return elapsed;
    }

    Duration total() const
    {
      return _total;
    }

    void reset()
    {
      _total = Duration::zero();
    }

  private:

    Duration _total = Duration::zero();

  };

} // end namespace Dune::PDELab

namespace fmt {

  template<>
  struct formatter<Dune::PDELab::LineSearchStrategy>
    : public formatter<std::string_view>
  {

    using Base = formatter<std::string_view>;

    // we just inherit the parsing from the original formatter

    template <typename FormatContext>
    auto format(Dune::PDELab::LineSearchStrategy strategy, FormatContext& ctx) {
      return Base::format(name(strategy),ctx);
    }

  };

} // end namespace fmt

#endif // DUNE_PDELAB_SOLVER_UTILITY_HH
