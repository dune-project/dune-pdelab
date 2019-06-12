// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_ISTL_SEQUENTIALCOMPONENTS_HH
#define DUNE_PDELAB_BACKEND_ISTL_SEQUENTIALCOMPONENTS_HH

#include <dune/istl/solvercategory.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solver.hh>
#include <dune/istl/preconditioner.hh>
#include <dune/istl/scalarproducts.hh>

#include <dune/logging/logger.hh>

#include <dune/pdelab/common/exceptions.hh>
#include <dune/pdelab/backend/common/interface.hh>
#include <dune/pdelab/backend/solver.hh>
#include <dune/pdelab/backend/istl/vector.hh>
#include <dune/pdelab/backend/istl/bcrsmatrix.hh>

namespace Dune::PDELab::ISTL::Experimental {

  template<typename Domain, typename Range>
  class LinearOperator
    : public Dune::LinearOperator<Domain,Range>
    , public virtual Backend::LinearSolverComponent<Domain,Range>
  {};

  template<typename Domain, typename Range>
  class Preconditioner
    : public Dune::Preconditioner<Domain,Range>
    , public virtual Backend::LinearSolverComponent<Domain,Range>
  {};

  template<typename Matrix>
  class MatrixBasedPreconditioner
    : public Preconditioner<typename Matrix::Domain,typename Matrix::Range>
    , public Backend::MatrixUsingLinearSolverComponent<Matrix>
  {

    using Base = Backend::MatrixUsingLinearSolverComponent<Matrix>;

  public:

    MatrixBasedPreconditioner(std::shared_ptr<typename Base::MatrixProvider> matrix_provider)
      : Base(std::move(matrix_provider))
    {}

  };

  template<typename GO>
  class OnTheFlyLinearOperator
    : public LinearOperator<typename GO::Traits::Domain,typename GO::Traits::Range>
  {

    using Base = LinearOperator<typename GO::Traits::Domain,typename GO::Traits::Range>;

  public:

    using Domain        = typename GO::Traits::Domain;
    using Range         = typename GO::Traits::Range;
    using Matrix        = typename GO::Traits::Jacobian;
    using Field         = typename Base::Field;
    using Real          = typename Base::Real;
    using OneStepMethod = typename Base::OneStepMethod;
    using domain_type   = Domain;
    using range_type    = Range;
    using field_type    = Field;

    OnTheFlyLinearOperator(GO& go)
      : _go(go)
    {}

    void apply(const Domain& x, Range& y) const override
    {
      y = 0.0;
      _go.jacobian_apply(x,y);
    }

    void applyscaleadd(Field alpha, const Domain& x, Range& y) const override
    {
      if (not _temp)
        _temp = std::make_shared<Range>(y);
      *_temp = 0.0;
      _go.jacobian_apply(x,*_temp);
      y.axpy(alpha,*_temp);
    }

    void setLinearizationPoint(const Domain& linearization_point, bool keep_matrix) override
    {
      if (isNonlinear(_go.localOperator()))
        _go.applyJacobianEngine()->setLinearizationPoint(linearization_point);
    }

    SolverCategory::Category category() const override
    {
      return SolverCategory::sequential;
    }

    void setOneStepMethod(std::shared_ptr<OneStepMethod> method) override
    {
      _go.setOneStepMethod(method);
    }

    int startStep(Real t0, Real dt) override
    {
      return _go.startStep(t0,dt);
    }

    bool acceptStage(int stage, const Domain& solution) override
    {
      return _go.acceptStage(stage,solution);
    }

  private:
    GO& _go;
    mutable std::shared_ptr<Range> _temp;
  };


  template<typename GO>
  class AssembledLinearOperator
    : public LinearOperator<typename GO::Traits::Domain,typename GO::Traits::Range>
    , public Backend::MatrixBasedLinearSolverComponent<typename GO::Traits::Jacobian>
  {

    using Base = LinearOperator<typename GO::Traits::Domain,typename GO::Traits::Range>;

  public:

    using Domain        = typename GO::Traits::Domain;
    using Range         = typename GO::Traits::Range;
    using Matrix        = typename GO::Traits::Jacobian;
    using Field         = typename Base::Field;
    using Real          = typename Base::Real;
    using OneStepMethod = typename Base::OneStepMethod;
    using domain_type   = Domain;
    using range_type    = Range;
    using field_type    = Field;

    AssembledLinearOperator(std::shared_ptr<GO> go)
      : _go(std::move(go))
    {
      if (not isNonlinear(_go->localOperator()))
        updateMatrix();
    }

    void setOneStepMethod(std::shared_ptr<OneStepMethod> method) override
    {
      _go->setOneStepMethod(method);
    }

    int startStep(Real t0, Real dt) override
    {
      _stage = 0;
      return _go->startStep(t0,dt);
    }

    bool acceptStage(int stage, const Domain& solution) override
    {
      _go->acceptStage(stage,solution);
      assert(_stage == stage - 1);
      if (_stage < stage)
      {
        if (not isNonlinear(_go->localOperator()))
          updateMatrix();
        _stage = stage;
        return true;
      }
      else
      {
        return false;
      }
    }

    void apply(const Domain& x, Range& y) const override
    {
      _jacobian->mv(x,y);
    }

    void applyscaleadd(Field alpha, const Domain& x, Range& y) const override
    {
      _jacobian->usmv(alpha,x,y);
    }

    void setLinearizationPoint(const Domain& linearization_point, bool keep_matrix) override
    {
      if (isNonlinear(_go->localOperator()) and not keep_matrix)
      {
        _go->jacobianEngine()->setLinearizationPoint(linearization_point);
        updateMatrix();
      }
    }

    SolverCategory::Category category() const override
    {
      return SolverCategory::sequential;
    }

    void updateMatrix(bool rebuild_pattern = false) override
    {
      if (rebuild_pattern)
        _jacobian.reset();
      ensureMatrix();
      _go->jacobian(*_jacobian);
      for (auto dependent : this->dependents())
        dependent->matrixUpdated(rebuild_pattern);
    }

    bool hasMatrix() const override
    {
      return bool(_jacobian);
    }

    const Matrix& matrix() const override
    {
      if (not _jacobian)
        DUNE_THROW(LinearAlgebraError, "Matrix has not been assembled yet");
      return *_jacobian;
    }

    const Matrix* matrixPointer() const override
    {
      return _jacobian.get();
    }

  private:

    void ensureMatrix()
    {
      if (not _jacobian)
        _jacobian = std::make_shared<Matrix>(*_go);
    }

    std::shared_ptr<GO> _go;
    mutable std::shared_ptr<Matrix> _jacobian;
    int _stage = 0;

  };

  enum class LinearizedOperatorMode
  {
    assembled,
    matrixFree,
    opportunistic
  };

  LinearizedOperatorMode parseLinearizedOperatorMode(std::string_view name)
  {
    if (name == "assembled")
      return LinearizedOperatorMode::assembled;
    if (name == "matrix-free")
      return LinearizedOperatorMode::matrixFree;
    if (name == "opportunistic")
      return LinearizedOperatorMode::opportunistic;
    DUNE_THROW(InvalidArgument,"unknown linearized operator mode: " << name);
  }

  template<typename GO>
  class LinearizedOperator
    : public LinearOperator<typename GO::Traits::Domain,typename GO::Traits::Range>
    , public Backend::MatrixBasedLinearSolverComponent<typename GO::Traits::Jacobian>
  {

    using Base = LinearOperator<typename GO::Traits::Domain,typename GO::Traits::Range>;
    using MatrixProviderBase = Backend::MatrixBasedLinearSolverComponent<typename GO::Traits::Jacobian>;

  public:

    using Domain        = typename GO::Traits::Domain;
    using Range         = typename GO::Traits::Range;
    using Matrix        = typename GO::Traits::Jacobian;
    using Field         = typename Base::Field;
    using Real          = typename Base::Real;
    using OneStepMethod = typename Base::OneStepMethod;
    using domain_type   = Domain;
    using range_type    = Range;
    using field_type    = Field;

    LinearizedOperator(std::shared_ptr<GO> go, LinearizedOperatorMode mode, bool allow_assembly, bool time_dependent = true)
      : _go(std::move(go))
      , _mode(mode)
      , _allow_assembly(allow_assembly)
      , _time_dependent(time_dependent)
    {
      if (_mode == LinearizedOperatorMode::assembled and not isNonlinear(_go->localOperator()))
        updateMatrix();
    }

    void setOneStepMethod(std::shared_ptr<OneStepMethod> method) override
    {
      _go->setOneStepMethod(method);
    }

    int startStep(Real t0, Real dt) override
    {
      _stage = 0;
      return _go->startStep(t0,dt);
    }

    bool acceptStage(int stage, const Domain& solution) override
    {
      _go->acceptStage(stage,solution);
      assert(_stage == stage - 1);
      if (_stage < stage)
      {
        if (not isNonlinear(_go->localOperator()) and _time_dependent)
          updateMatrix();
        _stage = stage;
        return true;
      }
      else
      {
        return false;
      }
    }

    void apply(const Domain& x, Range& y) const override
    {
      if (_mode == LinearizedOperatorMode::assembled or (_mode == LinearizedOperatorMode::opportunistic and _jacobian))
        _jacobian->mv(x,y);
      else
        _go.jacobian_apply(x,y);
    }

    void applyscaleadd(Field alpha, const Domain& x, Range& y) const override
    {
      if (_mode == LinearizedOperatorMode::assembled or (_mode == LinearizedOperatorMode::opportunistic and _jacobian))
        _jacobian->usmv(alpha,x,y);
      else
      {
        if (not _temp)
          _temp = std::make_shared<Range>(y);
        *_temp = 0.0;
        _go.jacobian_apply(x,*_temp);
        y.axpy(alpha,*_temp);
      }
    }

    void setLinearizationPoint(const Domain& linearization_point, bool keep_matrix) override
    {
      if (isNonlinear(_go->localOperator()))
      {
        _go->jacobianEngine()->setLinearizationPoint(linearization_point);
        if (_jacobian and not keep_matrix) {
          updateMatrix();
        }
      }
    }

    SolverCategory::Category category() const override
    {
      return SolverCategory::sequential;
    }

    void updateMatrix(bool rebuild_pattern = false) override
    {
      if (rebuild_pattern)
        _jacobian.reset();
      ensureMatrix();
      _log.trace("Assembling matrix"_fmt);
      _go->jacobian(*_jacobian);
      for (auto dependent : this->dependents())
        dependent->matrixUpdated(rebuild_pattern);
    }

    bool hasMatrix() const override
    {
      return bool(_jacobian);
    }

    void subscribe(typename MatrixProviderBase::Dependent& dependent) override
    {
      if (not _allow_assembly)
        DUNE_THROW(LinearAlgebraError, "Cannot subscribe to matrix of operator that does not provide a matrix");
      if (not _jacobian)
      {
        if (_mode != LinearizedOperatorMode::assembled)
          _log.warning("Requested matrix of linearized operator in matrix-free mode"_fmt);
        if (isNonlinear(_go->localOperator()))
          ensureMatrix();
        else
          updateMatrix();
      }
      MatrixProviderBase::subscribe(dependent);
    }

    const Matrix& matrix() const override
    {
      if (not _jacobian)
        DUNE_THROW(LinearAlgebraError, "Matrix has not been assembled yet");
      return *_jacobian;
    }

    const Matrix* matrixPointer() const override
    {
      return _jacobian.get();
    }

  private:

    void ensureMatrix()
    {
      if (not _jacobian)
      {
        _log.trace("Allocating matrix"_fmt);
        _jacobian = std::make_shared<Matrix>(*_go);
      }
    }

    std::shared_ptr<GO> _go;
    mutable std::shared_ptr<Range> _temp;
    mutable std::shared_ptr<Matrix> _jacobian;
    LinearizedOperatorMode _mode;
    bool _allow_assembly;
    bool _time_dependent;
    int _stage = 0;
    Logging::Logger _log;

  };


  template<typename GO>
  LinearizedOperator<GO> makeLinearizedOperator(std::shared_ptr<GO> go, const ParameterTree& params)
  {
    auto mode = LinearizedOperatorMode::assembled;
    if (params.hasKey("mode"))
      mode = parseLinearizedOperatorMode(params["mode"]);
    bool allow_assembly = params.get("allow_assembly",true);
    bool time_dependent = params.get("time_depdendent",true);
    return {go, mode, allow_assembly, time_dependent};
  }


  template<typename GO>
  class Operator
    : public ResidualEvaluator<typename GO::Traits::Domain,typename GO::Traits::Range>
  {

    using Base = ResidualEvaluator<typename GO::Traits::Domain,typename GO::Traits::Range>;

  public:

    using GridOperator = GO;
    using Domain = typename GridOperator::Traits::Domain;
    using Range  = typename GridOperator::Traits::Range;
    using Real   = typename Base::Real;
    using LinearizedOperator = Dune::PDELab::ISTL::Experimental::LinearizedOperator<GO>;

    void operator()(const Domain& domain, Range& range) override
    {
      _grid_operator->residual(domain,range);
    }

    int startStep(Real time, Real dt) override
    {
      return _grid_operator->startStep(time,dt);
    }

    bool acceptStage(int stage, const Domain& solution) override
    {
      return _grid_operator->acceptStage(stage,solution);
    }

    Range makeRange() override
    {
      return Range(_grid_operator->testSpace());
    }

    Domain makeDomain() override
    {
      return Domain(_grid_operator->trialSpace());
    }

    Operator(std::shared_ptr<GridOperator> grid_operator, const ParameterTree& params)
      : _grid_operator(std::move(grid_operator))
      , _params(params)
    {}

    friend LinearizedOperator linearize(const Operator& op)
    {
      return makeLinearizedOperator(op._grid_operator,op._params);
    }

  private:

    std::shared_ptr<GridOperator> _grid_operator;
    const ParameterTree& _params;

  };


  template<
    typename GridOperator
    >
  Operator(
    std::shared_ptr<GridOperator>
    )
    -> Operator<GridOperator>;

  /*
  template<typename GridOperator>
  auto makeOperator(std::shared_ptr<LO> lo, const Factory& factory, const ParameterTree& params)
  {
    return std::make_shared<
      MatrixFreeSequentialPreconditioner<
        typename LO::Domain,
        typename LO::Range,
        Factory
        >
      >(lo,factory,params);
  }
*/


  template<typename Domain_, typename Range_, typename Factory>
  class MatrixFreeSequentialPreconditioner
    : public Preconditioner<Domain_,Range_>
  {

  public:

    using Domain             = Domain_;
    using Range              = Range_;

    using ISTLPreconditioner = typename Dune::Preconditioner<Backend::Native<Domain>,Backend::Native<Range>>;
    using LinearOperator     = Dune::PDELab::ISTL::Experimental::LinearOperator<Domain_,Range_>;

    void pre(Domain& domain, Range& range) override
    {
      using Backend::native;
      _prec->pre(native(domain),native(range));
    }

    void apply(Domain& domain, const Range& range) override
    {
      using Backend::native;
      _prec->apply(native(domain),native(range));
    }

    void post(Domain& domain) override
    {
      using Backend::native;
      _prec->post(native(domain));
    }

    SolverCategory::Category category() const override
    {
      return SolverCategory::sequential;
    }

    MatrixFreeSequentialPreconditioner(
      std::shared_ptr<LinearOperator> lin_op,
      Factory factory,
      const ParameterTree& params
      )
      : _lin_op(std::move(lin_op))
      , _factory(std::move(factory))
      , _prec(_factory(_lin_op,params,std::shared_ptr<ISTLPreconditioner>(nullptr)))
    {}

  private:

    std::shared_ptr<LinearOperator> _lin_op;
    Factory _factory;
    std::shared_ptr<ISTLPreconditioner> _prec;

  };

  template<
    typename LinearOperator,
    typename Factory
    >
  MatrixFreeSequentialPreconditioner(
    std::shared_ptr<LinearOperator>,
    Factory,
    const ParameterTree&)
    -> MatrixFreeSequentialPreconditioner<
      typename LinearOperator::Domain,
      typename LinearOperator::Range,
      Factory
      >;


  template<typename LO, typename Factory>
  auto makeMatrixFreeSequentialPreconditioner(std::shared_ptr<LO> lo, const Factory& factory, const ParameterTree& params)
  {
    return std::make_shared<
      MatrixFreeSequentialPreconditioner<
        typename LO::Domain,
        typename LO::Range,
        Factory
        >
      >(lo,factory,params);
  }


  template<typename Matrix_, typename Factory>
  class MatrixBasedSequentialPreconditioner
    : public MatrixBasedPreconditioner<Matrix_>
  {

    Factory _factory;

  public:

    using Matrix         = Matrix_;
    using MatrixProvider = Backend::MatrixBasedLinearSolverComponent<Matrix>;
    using Domain         = typename Matrix::Domain;
    using Range          = typename Matrix::Range;

    using ISTLPreconditioner = typename Dune::Preconditioner<Backend::Native<Domain>,Backend::Native<Range>>;

    using ISTLPreconditionerPointer = decltype(
      _factory(
        std::declval<std::shared_ptr<MatrixProvider>>(),
        std::declval<std::shared_ptr<ISTLPreconditioner>>()
        )
      );

    void matrixUpdated(bool pattern_updated) override
    {
      _prec = _factory(this->matrixProviderPointer(),_prec);
    }

    void pre(Domain& domain, Range& range) override
    {
      using Backend::native;
      _prec->pre(native(domain),native(range));
    }

    void apply(Domain& domain, const Range& range) override
    {
      using Backend::native;
      _prec->apply(native(domain),native(range));
    }

    void post(Domain& domain) override
    {
      using Backend::native;
      _prec->post(native(domain));
    }

    SolverCategory::Category category() const override
    {
      return SolverCategory::sequential;
    }

    template<typename MatrixProviderPointer>
    MatrixBasedSequentialPreconditioner(
      MatrixProviderPointer matrix_provider,
      const Factory& factory,
      const ParameterTree& params
      )
      : _factory(factory)
      , _params(params)
    {
      this->setMatrixProvider(matrix_provider);
    }

  private:

    std::shared_ptr<ISTLPreconditioner> _prec;
    const ParameterTree& _params;

  };

  template<
    typename MatrixProviderPointer,
    typename Factory
    >
  MatrixBasedSequentialPreconditioner(
    MatrixProviderPointer,
    const Factory&,
    const ParameterTree&)
    -> MatrixBasedSequentialPreconditioner<
      typename MatrixProviderPointer::element_type::Matrix,
      Factory
      >;


  template<typename LO, typename Factory>
  auto makeMatrixBasedSequentialPreconditioner(std::shared_ptr<LO> lo, const Factory& factory, const ParameterTree& params)
  {
    return std::make_shared<
      MatrixBasedSequentialPreconditioner<
        typename LO::Matrix,
        Factory
        >
      >(lo,factory,params);
  }

  template<typename Domain_, typename Range_, typename Factory>
  class IterativeLinearSolver
    : public Backend::LinearSolver<Domain_,Range_>
  {

    using Base = Backend::LinearSolver<Domain_,Range_>;

  public:

    using Domain               = Domain_;
    using Range                = Range_;
    using ISTLSolver           = Dune::IterativeSolver<Domain,Range>;
    using ISTLPreconditioner   = Dune::Preconditioner<Domain,Range>;
    using ISTLScalarProduct    = Dune::ScalarProduct<Domain>;
    using Preconditioner       = Dune::PDELab::ISTL::Experimental::Preconditioner<Domain,Range>;
    using LinearOperator       = Dune::PDELab::ISTL::Experimental::LinearOperator<Domain,Range>;
    using Real                 = typename Base::Real;
    using OneStepMethod        = typename Base::OneStepMethod;

    int startStep(Real t0, Real dt) override
    {
      int stages = _linear_operator->startStep(t0,dt);
      _preconditioner->startStep(t0,dt);
      return stages;
    }

    bool acceptStage(int stage, const Domain& solution) override
    {
      bool changed = _linear_operator->acceptStage(stage,solution);
      changed |= _preconditioner->acceptStage(stage,solution);
      return changed;
    }

    void setOneStepMethod(std::shared_ptr<OneStepMethod> method) override
    {
      _linear_operator->setOneStepMethod(method);
      _preconditioner->setOneStepMethod(method);
    }

    void setLinearizationPoint(const Domain& linearization_point, bool keep_matrix) override
    {
      _linear_operator->setLinearizationPoint(linearization_point,keep_matrix);
      _preconditioner->setLinearizationPoint(linearization_point,keep_matrix);
    }

    void solve(Domain& solution, Range& rhs) override
    {
      Dune::InverseOperatorResult stat;
      _solver->apply(solution,rhs,stat);
    }

    void solve(Domain& solution, Range& rhs, Real reduction) override
    {
      Dune::InverseOperatorResult stat;
      _solver->apply(solution,rhs,static_cast<double>(reduction),stat);
    }

    Real norm(const Range& range) const override
    {
      return range.two_norm();
    }

    IterativeLinearSolver(
      std::shared_ptr<LinearOperator> linear_operator,
      std::shared_ptr<Preconditioner> preconditioner,
      const Factory& factory,
      const ParameterTree& params,
      Logging::Logger log
      )
      : _factory(factory)
      , _linear_operator(std::move(linear_operator))
      , _preconditioner(std::move(preconditioner))
      , _params(params)
      , _log(log)
      , _solver(_factory(_linear_operator,_preconditioner,_params,_log))
    {}

  private:

    Factory _factory;
    std::shared_ptr<LinearOperator> _linear_operator;
    std::shared_ptr<Preconditioner> _preconditioner;
    const ParameterTree& _params;
    Logging::Logger _log;
    std::shared_ptr<ISTLSolver> _solver;

  };


  template<
    typename LinearOperatorPointer,
    typename PreconditionerPointer,
    typename Factory
    >
  IterativeLinearSolver(
    LinearOperatorPointer,
    PreconditionerPointer,
    const Factory&,
    const ParameterTree&)
    -> IterativeLinearSolver<
      typename LinearOperatorPointer::element_type::Domain,
      typename LinearOperatorPointer::element_type::Range,
      Factory
      >;


  template<typename LO, typename Preconditioner, typename Factory>
  auto makeIterativeLinearSolver(
    std::shared_ptr<LO> lo,
    std::shared_ptr<Preconditioner> preconditioner,
    const Factory& factory,
    const ParameterTree& params)
  {
    return std::make_shared<
      IterativeLinearSolver<
        typename LO::Domain,
        typename LO::Range,
        Factory
        >
      >(lo,preconditioner,factory,params,Logging::Logger());
  }

} // namespace Dune::PDELab::ISTL::Experimental

#endif // DUNE_PDELAB_BACKEND_ISTL_SEQUENTIALCOMPONENTS_HH
