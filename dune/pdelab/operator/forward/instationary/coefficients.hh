#ifndef DUNE_PDELAB_OPERATOR_FORWARD_INSTATIONARY_COEFFICIENTS_HH
#define DUNE_PDELAB_OPERATOR_FORWARD_INSTATIONARY_COEFFICIENTS_HH

#include <dune/pdelab/instationary/onestepparameter.hh>

#include <dune/common/float_cmp.hh>

#if __has_include(<format>)
#include <format>
#endif
#include <array>
#include <string>

namespace Dune::PDELab::inline Experimental {

//! Object to store coefficients for instationary problems in residual form.
/**
 * The parameters \f$ \alpha,\beta \in \mathbb{R}^{s\times r} \f$ and
 * \f$ d\in \mathbb{R}^r \f$ represent the coefficients of a system of ODE in the following form:
 * \f[
 * \begin{aligned}
 *   \sum_{j=0}^s \left[ \hat{\alpha}_{ij} m_h\left(u_h^{(j)}, v; t^k + d_j\Delta t^k\right)
 *     + \beta_{ij}\Delta t^k r_h \left( u_h^{(j)},v,t^k+d_j\Delta t^k \right)\right] &= 0 & \forall i=1,\ldots,s \quad \forall v\in V_h(t^{k+1})\\
 * \end{aligned}
 * \f]
 * where \f$ m_h\f$ is the temporal residual form (mass operator), \f$
 * r_h \f$ the spatial residual form (stiffness operator), and \f$ \hat{\alpha} := (I-\alpha)\f$.
 * Most notably, these inestationary coefficients may be used to implement a
 * Runge-Kutta method in the modified Shu-Osher form [1], embedded methods, IMEX methods,
 * and multi-stage methods.
 *
 * [1] Chi W. Shu and Stanley Osher. Efficient implementation of essentially
 * non- oscillatory shock-capturing schemes. J. Comput. Phys., 77:439–471
 *
 */
class InstationaryCoefficients {
public:
  enum class Type {Explicit, SemiImplicit, FullyImplicit};

  using MassWeight = double;
  using StiffnessWeight = double;
  using TimeWeight = double;

  InstationaryCoefficients() {}

  // The stage id might be used to distinguish different kinds of methods when they are mereged into a single object
  // (e.g., assign 0/1 ids to implicit and explicit parts in an IMEX method).
  InstationaryCoefficients(
    const std::vector<std::vector<double>>& mass_weight,
    const std::vector<std::vector<double>>& stiffness_weight,
    const std::vector<double>& time_weight,
    const std::vector<int>& stage_id
  ) : _mass_weight(mass_weight)
    , _stiffness_weight(stiffness_weight)
    , _time_weight(time_weight)
    , _stage_id(stage_id)
  {
    assert(_mass_weight.size() == _stiffness_weight.size());
    assert(_mass_weight.size() == _stage_id.size());
    for (std::size_t i = 0; i != _mass_weight.size(); ++i) {
      assert(_mass_weight[i].size() == _stiffness_weight[i].size());
      assert(_mass_weight[i].size() == _time_weight.size());
    }
  }

  InstationaryCoefficients(
    const std::vector<std::vector<double>>& mass_weight,
    const std::vector<std::vector<double>>& stiffness_weight,
    const std::vector<double>& time_weight,
    int stage_id = 0
  ) : InstationaryCoefficients(mass_weight, stiffness_weight, time_weight, std::vector<int>(mass_weight.size(), stage_id))
  {}

  template<std::convertible_to<double> T>
  InstationaryCoefficients(const TimeSteppingParameterInterface<T>& method)
    : _mass_weight(method.s())
    , _stiffness_weight(method.s())
    , _time_weight(method.s()+1)
    , _stage_id(method.s(), 0)
  {
    for (std::size_t i = 0; i != method.s(); ++i) {
      _mass_weight[i].resize(method.s()+1);
      _stiffness_weight[i].resize(method.s()+1);
      for (std::size_t j = 0; j != method.s()+1; ++j) {
        _time_weight[j] = static_cast<double>(method.d(j));
        _mass_weight[i][j] = static_cast<double>(method.a(i+1,j));
        _stiffness_weight[i][j] = static_cast<double>(method.b(i+1,j));
      }
    }
    check();
  }

  std::size_t extent(std::size_t i) const {
    if (i == 0)
      return _mass_weight.size();
    else if (i == 1)
      return _time_weight.size();
    else
      DUNE_THROW(RangeError, "InstationaryCoefficients has rank 2");
  };

  bool doMass(std::size_t i, std::size_t j) const {
    return FloatCmp::ne(massWeight(i,j), double{0});
  }

  bool doStiffness(std::size_t i, std::size_t j) const {
    return FloatCmp::ne(stiffnessWeight(i,j), double{0});
  }

  void check() const {
    assert(extent(0) > 0);
    assert(extent(1) > 0);
    assert(_mass_weight.size() == extent(0));
    assert(_stiffness_weight.size() == extent(0));
    for (std::size_t i = 0; i < extent(0); ++i) {
      assert(_mass_weight[i].size() == extent(1));
      assert(_stiffness_weight[i].size() == extent(1));
      double acc = 0.;
      for (auto val : _mass_weight[i]) acc += val;
      assert(bool(FloatCmp::eq<double, FloatCmp::CmpStyle::absolute>(acc, MassWeight{0.}, 1e-15)));
      assert(FloatCmp::ge(timeWeight(i), TimeWeight{0.}));
      assert(FloatCmp::le(timeWeight(i), TimeWeight{1.}));
    }

    assert(FloatCmp::eq(timeWeight(0), TimeWeight{0}));
    assert(FloatCmp::eq(timeWeight(extent(1)-1), TimeWeight{1.}));
  }


  MassWeight massWeight(std::size_t i, std::size_t j) const {
    assert(i < extent(0));
    assert(j < extent(1));
    return _mass_weight[i][j];
  }

  StiffnessWeight stiffnessWeight(std::size_t i, std::size_t j) const {
    assert(i < extent(0));
    assert(j < extent(1));
    return _stiffness_weight[i][j];
  }

  TimeWeight timeWeight(std::size_t i) const {
    assert(i < extent(1));
    return _time_weight[i];
  }

  int id(std::size_t i) const {
    assert(i < extent(0));
    return _stage_id[i];
  }

  auto type() const {
    check();

    if (extent(0)+1 != extent(1))
      DUNE_THROW(RangeError, "");

    for (std::size_t i = 0; i != extent(0); ++i) {
      for (std::size_t j = i+2; j < extent(1); ++j)
        if (doMass(i,j) or doStiffness(i,j))
          return Type::FullyImplicit;
    }

    for (std::size_t i = 1; i != extent(0); ++i)
      if (doMass(i,i+1) or doStiffness(i,i+1))
        return Type::SemiImplicit;

    return Type::Explicit;
  }

  InstationaryCoefficients slice(std::size_t row_begin, std::size_t row_size, std::size_t col_begin, std::size_t col_size) const {

    assert(row_begin + row_size <= extent(0));
    assert(col_begin + col_size <= extent(1));

    std::vector<std::vector<double>> mass_weight(row_size);
    std::vector<std::vector<double>> stiffness_weight(row_size);

    for(std::size_t row = 0; row != row_size; ++row) {
      auto mass_begin = std::begin(_mass_weight[row_begin+row]) + col_begin;
      auto stiffness_begin = std::begin(_stiffness_weight[row_begin+row]) + col_begin;
      mass_weight[row] = std::vector<double>(mass_begin, mass_begin + col_size);
      stiffness_weight[row] = std::vector<double>(stiffness_begin, stiffness_begin + col_size);
    }
    auto time_begin = std::begin(_time_weight) + col_begin;
    std::vector<double> time_weight(time_begin, time_begin + col_size);
    auto stage_id_begin = std::begin(_stage_id) + row_begin;
    std::vector<int> stage_id(stage_id_begin, stage_id_begin + row_size);
    return InstationaryCoefficients{
      mass_weight,
      stiffness_weight,
      time_weight,
      stage_id
    };
  }

  void print() const {
#if __has_include(<format>)
    auto hline = [&]{std::cout << std::format("{0:>8}+{1:->{2}}+\n", ' ','-', 8*extent(1));};
    hline();
    for (std::size_t row = 0; row != extent(0); ++row) {
      std::cout << ((row == extent(0)/2) ? "(I-α):= |" : "        |");
      for (std::size_t col = 0; col != extent(1); ++col)
        std::cout << std::format(" {: 1.3f}{}", massWeight(row,col), (col+1 == extent(1)) ? " |\n" : ",");
    }
    hline();
    for (std::size_t row = 0; row != extent(0); ++row) {
      std::cout << ((row == extent(0)/2) ? "    β:= |" : "        |");
      for (std::size_t col = 0; col != extent(1); ++col)
        std::cout << std::format(" {: 1.3f}{}", stiffnessWeight(row,col), (col+1 == extent(1)) ? " |\n" : ",");
    }
    hline();
    std::cout << "  d^T:= |";
    for (std::size_t col = 0; col != extent(1); ++col)
      std::cout << std::format(" {: 1.3f}{}", timeWeight(col), (col+1 == extent(1)) ? " |\n" : ",");
    hline();
#else
    DUNE_THROW(
      NotImplemented,
      "To print this class you need an standard library with std::format");
#endif
  }

private:
  std::vector<std::vector<double>> _mass_weight;
  std::vector<std::vector<double>> _stiffness_weight;
  std::vector<double> _time_weight;
  std::vector<int> _stage_id;
  // std::string name;
};

inline auto SSPRK33() {
  return InstationaryCoefficients{
    // mass
    {{ {{   -1,     1,     0,  0}},
       {{-3./4, -1./4,     1,  0}},
       {{-1./3,     0, -2./3,  1}} }},
    // stiffness
    {{ {{1,    0,    0,  0}},
       {{0, 1./4,    0,  0}},
       {{0,    0, 2./3,  0}} }},
    // time
       {{0,    1,   .5,  1}}
  };
}

// https://doi.org/10.1137/S00361429013890
inline auto SSPRK54() {
  return InstationaryCoefficients{
        { {                -1,                 1,                 0,                 0,                 0, 0},
          {-0.444370493651235,-0.555629506348765,                 1,                 0,                 0, 0},
/*(I-α)*/ {-0.620101851488403,                 0,-0.379898148511597,                 1,                 0, 0},
          {-0.178079954393132,                 0,                 0,-0.821920045606868,                 1, 0},
          {                 0,                 0,-0.517231671970585,-0.096059710526147,-0.386708617503269, 1} },

        { { 0.391752226571890,                 0,                 0,                 0,                 0, 0},
          {                 0, 0.368410593050371,                 0,                 0,                 0, 0},
/*  β:=*/ {                 0,                 0, 0.251891774271694,                 0,                 0, 0},
          {                 0,                 0,                 0, 0.544974750228521,                 0, 0},
          {                 0,                 0,                 0, 0.063692468666290, 0.226007483236906, 0} },

/*d^T:=*/ {                 0,  0.39175222700392,  0.58607968896779, 0.47454236302687,   0.93501063100924, 1}
  };
}

} // namespace Dune::PDELab::inline Experimental

#endif // DUNE_PDELAB_OPERATOR_FORWARD_INSTATIONARY_COEFFICIENTS_HH
