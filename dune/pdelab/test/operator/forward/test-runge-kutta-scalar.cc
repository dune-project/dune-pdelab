#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/pdelab/operator/forward/runge_kutta.hh>
#include <dune/pdelab/operator/adapter.hh>

#include <gtest/gtest.h>

#include <cmath>

using namespace Dune::PDELab::Experimental;

using Domain = double;
using Range = double;
using TimeQuantity = double;
using DurationQuantity = double;

using RungeKuttaDomain = std::vector<Domain>;
using RungeKuttaRange = std::vector<Range>;

TEST(TestRungeKutta, ExplicitScalar) {

  auto stiffness = [](TimeQuantity t, Domain y) {
    return y - t*t + 1;
  };

  auto sol = [](TimeQuantity t){
    return t*t + 2*t + 1 - 0.5 * std::exp(t);
  };

  std::shared_ptr<Operator<RungeKuttaDomain,RungeKuttaRange>> forward = std::make_shared<OperatorAdapter<RungeKuttaDomain,RungeKuttaRange>>(
    [&](Operator<RungeKuttaDomain,RungeKuttaRange>& f, const RungeKuttaDomain& y, RungeKuttaRange& residual) -> ErrorCondition {
      auto time = f.get<TimeQuantity>("time");
      auto dt = f.get<DurationQuantity>("duration");
      auto tableau = f.get<InstationaryCoefficients>("instationary_coefficients");

      assert(residual.size() == tableau.extent(0));
      assert(y.size() == tableau.extent(1));
      for (std::size_t s = 0; s != residual.size(); ++s) {
        for (std::size_t k = 0; k != tableau.extent(1); ++k) {
          auto stime = time + dt*tableau.timeWeight(k);
          auto a_ij = tableau.massWeight(s, k);
          auto b_ij = tableau.stiffnessWeight(s, k);
          residual[s] += b_ij * dt * stiffness(stime, y[k]) -  a_ij * y[k];
        }
      }
      return {};
    });

  std::shared_ptr<Operator<RungeKuttaRange,RungeKuttaDomain>> inverse = std::make_shared<OperatorAdapter<RungeKuttaRange,RungeKuttaDomain>>(
    [&](Operator<RungeKuttaRange,RungeKuttaDomain>& inv, const RungeKuttaRange& residual, RungeKuttaDomain& y) -> ErrorCondition {
      // the following only works if problem is explicit and the mass term is strictly diagonal!
      auto tableau = inv.get<InstationaryCoefficients>("forward.instationary_coefficients");
      if (y.size() != 1) DUNE_THROW(Dune::NotImplemented,"");
      if (tableau.extent(0) != 1) DUNE_THROW(Dune::NotImplemented,"");
      if (tableau.extent(1) != 1) DUNE_THROW(Dune::NotImplemented,"");
      if (Dune::FloatCmp::ne(tableau.stiffnessWeight(0,0), 0.)) DUNE_THROW(Dune::NotImplemented,"");
      y[0] = residual[0]/tableau.massWeight(0,0);
      return {};
    });

  RungeKutta<RungeKuttaDomain,RungeKuttaRange> runge_kutta;
  TimeQuantity time{0.};
  runge_kutta["initial_residual"] = Range{0.};
  runge_kutta["inverse"] = inverse;
  runge_kutta["inverse.forward"] = std::weak_ptr(forward);
  runge_kutta["inverse.forward.time"] = std::ref(time);
  runge_kutta["inverse.forward.duration"] = DurationQuantity{0.05};
  runge_kutta["instationary_coefficients"] = InstationaryCoefficients{Dune::PDELab::RK4Parameter<double>{}};

  Domain x0 = 0.5, x1;
  std::cout << std::format("Time      Kunge-Kutta    Analytical") << std::endl;
  std::cout << std::format("{:1.4f}    {:1.4f}         {:1.4f}", time, x0, sol(time)) << std::endl;
  while(Dune::FloatCmp::lt(time, TimeQuantity{2.})) {
    runge_kutta.apply(x1 = x0, x0).or_throw();
    EXPECT_NEAR(x0, sol(time), 1e-6);
    std::cout << std::format("{:1.4f}    {:1.4f}         {:1.4f}", time, x0, sol(time)) << std::endl;
  }
  std::cout << runge_kutta << std::endl;
}
