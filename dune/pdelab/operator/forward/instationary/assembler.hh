#ifndef DUNE_PDELAB_OPERATOR_FORWARD_INSTATIONARY_ASSEMBLER_HH
#define DUNE_PDELAB_OPERATOR_FORWARD_INSTATIONARY_ASSEMBLER_HH

#include <dune/pdelab/operator/forward/instationary/assembler/forward.hh>
#include <dune/pdelab/operator/forward/instationary/assembler/jacobian.hh>
#include <dune/pdelab/operator/forward/instationary/assembler/jacobian_apply.hh>


namespace Dune::PDELab::inline Experimental {

template<class          Coefficients,
         class          Residual,
         Concept::Basis TrialBasis,
         Concept::Basis TestBasis,
         class          MassLocalOperator,
         class          StiffnessLocalOperator,
         class TimeQuantity = double,
         class DurationQuantity  = double,
         DurationPosition dt_position = DurationPosition::StiffnessNumerator>
class InstationaryMatrixFreeAssembler : public InstationaryForwardAssembler<Coefficients, Residual, TrialBasis, TestBasis, MassLocalOperator, StiffnessLocalOperator, TimeQuantity, DurationQuantity, dt_position>
{
  using Base = InstationaryForwardAssembler<Coefficients, Residual, TrialBasis, TestBasis, MassLocalOperator, StiffnessLocalOperator, TimeQuantity, DurationQuantity, dt_position>;
public:
  using Base::Base;

  virtual std::shared_ptr<Operator<Coefficients,Residual>> derivative(const Coefficients& x, std::shared_ptr<Operator<Coefficients,Residual>> reuse_dx = nullptr) const override {
    using Type = InstationaryJacobianApplyAssembler<Coefficients, Residual, TrialBasis, TestBasis, MassLocalOperator, StiffnessLocalOperator, TimeQuantity, DurationQuantity, dt_position>;
    std::shared_ptr<Type> dx;
    if (reuse_dx)
      dx = std::dynamic_pointer_cast<Type>(reuse_dx);
    else
      dx = std::make_unique<Type>(this->_trial, this->_test, this->_mass_lop, this->_stiff_lop);

    // TODO: assert that old basis and local operators are same as for this

    dx->get("time")                       = this->get("time");
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
         class          StiffnessLocalOperator>
auto makeInstationaryMatrixFreeAssembler(const TrialBasis& trial,
                                         const TestBasis& test,
                                         const MassLocalOperator& mass_lop,
                                         const StiffnessLocalOperator& stiff_lop)
{
  using Type = InstationaryMatrixFreeAssembler<Coefficients, Residual, TrialBasis, TestBasis, MassLocalOperator, StiffnessLocalOperator, double, double, DurationPosition::StiffnessNumerator>;
  return std::make_unique<Type>(trial, test, mass_lop, stiff_lop);
}



template<class          Coefficients,
         class          Residual,
         Concept::Basis TrialBasis,
         Concept::Basis TestBasis,
         class          MassLocalOperator,
         class          StiffnessLocalOperator,
         class          MatrixContainer,
         class TimeQuantity = double,
         class DurationQuantity  = double,
         DurationPosition dt_position = DurationPosition::StiffnessNumerator>
class InstationaryMatrixBasedAssembler : public InstationaryForwardAssembler<Coefficients, Residual, TrialBasis, TestBasis, MassLocalOperator, StiffnessLocalOperator, TimeQuantity, DurationQuantity, dt_position>
{
  using Base = InstationaryForwardAssembler<Coefficients, Residual, TrialBasis, TestBasis, MassLocalOperator, StiffnessLocalOperator, TimeQuantity, DurationQuantity, dt_position>;
public:

  InstationaryMatrixBasedAssembler(const TrialBasis& trial,
                        const TestBasis& test,
                        const MassLocalOperator& mass_lop,
                        const StiffnessLocalOperator& stiff_lop,
                        std::function<void(const Operator<Coefficients,Residual>&, MatrixContainer&)> container_resize = {})
    : Base{ trial, test, mass_lop, stiff_lop }
    , _container_resize{std::move(container_resize)}
  {}

  virtual std::shared_ptr<Operator<Coefficients,Residual>> derivative(const Coefficients& x, std::shared_ptr<Operator<Coefficients,Residual>> reuse_dx = nullptr) const override {
    using Jacobian = InstationaryJacobianAssembler<Coefficients, Residual, TrialBasis, TestBasis, MassLocalOperator, StiffnessLocalOperator, MatrixContainer, TimeQuantity, DurationQuantity, dt_position>;
    std::shared_ptr<Jacobian> dx;
    if (reuse_dx)
      dx = std::dynamic_pointer_cast<Jacobian>(reuse_dx);
    else
      dx = std::make_unique<Jacobian>(this->_trial, this->_test, this->_mass_lop, this->_stiff_lop);

    dx->get("time")                       = this->get("time");
    dx->get("duration")                   = this->get("duration");
    dx->get("instationary_coefficients")  = this->get("instationary_coefficients");

    if (not reuse_dx)
      _container_resize(*dx, dx->template get<MatrixContainer>("container"));

    dx->linearize(x);
    return dx;
  }

private:
  std::function<void(const Operator<Coefficients,Residual>&, MatrixContainer&)> _container_resize;
};


template<class          Coefficients,
         class          Residual,
         class          MatrixContainer,
         Concept::Basis TrialBasis,
         Concept::Basis TestBasis,
         class          MassLocalOperator,
         class          StiffnessLocalOperator>
auto makeInstationaryMatrixBasedAssembler(const TrialBasis& trial,
                                          const TestBasis& test,
                                          const MassLocalOperator& mass_lop,
                                          const StiffnessLocalOperator& stiff_lop,
                                          std::function<void(const Operator<Coefficients,Residual>&, MatrixContainer&)> container_resize = {})
{
  using Type = InstationaryMatrixBasedAssembler<Coefficients, Residual, TrialBasis, TestBasis, MassLocalOperator, StiffnessLocalOperator, MatrixContainer, double, double, DurationPosition::StiffnessNumerator>;
  return std::make_unique<Type>(trial, test, mass_lop, stiff_lop, container_resize);
}

} // namespace Dune::PDELab::inline Experimental

#endif // DUNE_PDELAB_OPERATOR_FORWARD_INSTATIONARY_ASSEMBLER_HH
