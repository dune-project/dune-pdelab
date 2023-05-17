#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/pdelab/operator/operator.hh>

#include <dune/pdelab/operator/inverse/newton.hh>
#include <dune/pdelab/operator/adapter.hh>

#include <gtest/gtest.h>

#include <cmath>

using namespace Dune::PDELab::Experimental;

struct X2 : public Operator<double, double> {

  X2() = default;
  X2(const X2&) = default;
  X2& operator=(const X2&) = default;

  ErrorCondition apply(const double& x, double& y) override {
    y += x*x;
    return {};
  }

  std::shared_ptr<Operator<double, double>> derivative(const double& x, std::shared_ptr<Operator<double, double>> dx = nullptr) const override {
    if (not dx)
      dx = std::make_shared<OperatorAdapter<double,double>>(
        [](Operator<double,double>& dxdy_inv, const double& u, double& z) -> ErrorCondition {
          z += 2*u;
          return {};
        });
    dx->get("linearization_point") = x;
    return dx;
  };
};


TEST(TestNewton, X2) {
  Newton<double,double,double> newton;
  // find root of x^2=a. See https://en.wikipedia.org/wiki/Newton%27s_method#Square_root
  double a = 612.;
  double abs_tolerance = 1e-5;
  double x;

  try {
    newton["forward"] = std::shared_ptr<Operator<double, double>>(std::make_shared<X2>());
    newton.setNorm([](double v){return std::abs(v);});
    std::shared_ptr<Operator<double,double>> jac_inv_op = std::make_shared<OperatorAdapter<double,double>>(
      [](Operator<double,double>& dxdy_inv, const double& y, double& z) -> ErrorCondition {
        auto x = dxdy_inv.get<double>("forward.linearization_point");
        z += y/(2*x);
        return {};
      });
    newton["dx_inverse"] = jac_inv_op;

    EXPECT_EQ(newton.apply(-a, x = 0), make_error_condition(Convergence::Reason::DivergedByNanOrInf));

    newton["convergence_condition.absolute_tolerance"] = abs_tolerance;

    newton.apply(-a, x = 1).or_throw();
    EXPECT_NEAR(x, std::sqrt(a), abs_tolerance);

    newton.apply(-a, x = -1).or_throw();
    EXPECT_NEAR(x, -std::sqrt(a), abs_tolerance);
  } catch (...) {
    std::cout << newton << std::endl;
    throw;
  }
  std::cout << newton << std::endl;
}
