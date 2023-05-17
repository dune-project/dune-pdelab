#ifndef DUNE_PDELAB_OPERATOR_FORWARD_INSTATIONARY_ASSEMBLER_HH
#define DUNE_PDELAB_OPERATOR_FORWARD_INSTATIONARY_ASSEMBLER_HH

#include <dune/pdelab/operator/forward/instationary/assembler/forward.hh>
#include <dune/pdelab/operator/forward/instationary/assembler/jacobian_apply.hh>


namespace Dune::PDELab::inline Experimental {

template<class          Coefficients,
         class          Residual,
         Concept::Basis TrialBasis,
         Concept::Basis TestBasis,
         class          MassLocalOperator,
         class          StiffnessLocalOperator,
         class TimePoint = double,
         class Duration  = double,
         DurationPosition dt_position = DurationPosition::StiffnessNumerator>
class InstationaryMatrixFreeAssembler : public InstationaryForwardAssembler<Coefficients, Residual, TrialBasis,TestBasis, MassLocalOperator, StiffnessLocalOperator, TimePoint, Duration, dt_position>
{
  using Base = InstationaryForwardAssembler<Coefficients, Residual, TrialBasis,TestBasis, MassLocalOperator, StiffnessLocalOperator, TimePoint, Duration, dt_position>;
public:
  using Base::Base;

  virtual std::shared_ptr<Operator<Coefficients,Residual>> derivative(const Coefficients& x, std::shared_ptr<Operator<Coefficients,Residual>> reuse_dx = nullptr) const override {
    using Type = InstationaryJacobianApplyAssembler<Coefficients, Residual, TrialBasis,TestBasis, MassLocalOperator, StiffnessLocalOperator, TimePoint, Duration, dt_position>;
    std::shared_ptr<Type> dx;
    if (reuse_dx)
      dx = std::dynamic_pointer_cast<Type>(reuse_dx);
    else
      dx = std::make_unique<Type>(this->_trial, this->_test, this->_mass_lop, this->_stiff_lop);

    dx->get("time_point")                 = this->get("time_point");
    dx->get("duration")                   = this->get("duration");
    dx->get("instationary_coefficients")  = this->get("instationary_coefficients");
    dx->get("linearization_point") = x;
    return dx;
  }
};

template<class          Coefficients,
         class          Residual,
         Concept::Basis TrialBasis,
         Concept::Basis TestBasis,
         class          MassLocalOperator,
         class          StiffnessLocalOperator,
         class TimePoint = double,
         class Duration  = double,
         DurationPosition dt_position = DurationPosition::StiffnessNumerator>
auto makeInstationaryAssembler(const TrialBasis& trial,
                               const TestBasis& test,
                               const MassLocalOperator& mass_lop,
                               const StiffnessLocalOperator& stiff_lop)
{
  using Type = InstationaryMatrixFreeAssembler<Coefficients, Residual, TrialBasis,TestBasis, MassLocalOperator, StiffnessLocalOperator, TimePoint, Duration, dt_position>;
  return std::make_unique<Type>(trial, test, mass_lop, stiff_lop);
}

} // namespace Dune::PDELab::inline Experimental

#endif // DUNE_PDELAB_OPERATOR_FORWARD_INSTATIONARY_ASSEMBLER_HH
