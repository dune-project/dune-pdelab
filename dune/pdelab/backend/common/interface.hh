#ifndef DUNE_PDELAB_BACKEND_COMMON_INTERFACE_HH
#define DUNE_PDELAB_BACKEND_COMMON_INTERFACE_HH

#include <memory>
#include <typeinfo>
#include <unordered_set>

#include <dune/common/ftraits.hh>
#include <dune/pdelab/instationary/onestepmethods.hh>

namespace Dune::PDELab::Backend {


  template<typename Domain_>
  struct TimestepAwareComponent
  {

    using Domain        = Domain_;
    using Real          = typename Dune::FieldTraits<typename Domain::ElementType>::real_type;
    using OneStepMethod = OneStep::Method<Real>;

    virtual ~TimestepAwareComponent()
    {}

    virtual int startStep(Real t0, Real dt)
    {
      return -1;
    }

    virtual bool acceptStage(int stage, const Domain& solution)
    {
      return false;
    }

    virtual void setOneStepMethod(std::shared_ptr<OneStepMethod>)
    {}

  };


  template<typename Domain_, typename Range_>
  struct LinearSolverComponent
    : public virtual TimestepAwareComponent<Domain_>
  {

    using Domain = Domain_;
    using Range  = Range_;
    using Field  = typename Domain::ElementType;
    using Real   = typename Dune::FieldTraits<Field>::real_type;

    virtual void setLinearizationPoint(const Domain& domain, bool keep_matrix)
    {}

  };

  template<typename Matrix_>
  class MatrixUsingLinearSolverComponent;

  template<typename Matrix_>
  class MatrixBasedLinearSolverComponent
    : public virtual LinearSolverComponent<typename Matrix_::Domain,typename Matrix_::Range>
  {

  public:

    using Matrix     = Matrix_;
    using Dependent  = MatrixUsingLinearSolverComponent<Matrix>;
    using Dependents = std::unordered_set<Dependent*>;

    virtual void updateMatrix(bool rebuild_pattern = false)
    {}

    virtual bool hasMatrix() const = 0;

    virtual const Matrix& matrix() const = 0;

    virtual const Matrix* matrixPointer() const
    {
      return nullptr;
    }

    void subscribe(Dependent& dependent)
    {
      auto& ti = typeid(dependent);
      _dependents.insert(&dependent);
      if (hasMatrix())
        dependent.matrixUpdated(true);
    }

    void unsubscribe(Dependent& dependent)
    {
      _dependents.erase(&dependent);
    }

    bool subscribed(const Dependent& dependent) const
    {
      return _dependents.count(&dependent) > 0;
    }

    const Dependents& dependents() const
    {
      return _dependents;
    }

  private:

    Dependents _dependents;

  };

  template<typename Matrix_>
  class MatrixUsingLinearSolverComponent
    : public virtual LinearSolverComponent<typename Matrix_::Domain,typename Matrix_::Range>
  {

  public:

    using Matrix         = Matrix_;
    using MatrixProvider = MatrixBasedLinearSolverComponent<Matrix>;

    virtual void matrixUpdated(bool pattern_updated)
    {}

    MatrixProvider& matrixProvider() const
    {
      return *_matrix_provider;
    }

    const Matrix& matrix() const
    {
      return matrixProvider().matrix();
    }

    ~MatrixUsingLinearSolverComponent()
    {
      if (_matrix_provider)
        _matrix_provider->unsubscribe(*this);
    }

  protected:

    MatrixUsingLinearSolverComponent(std::shared_ptr<MatrixProvider> matrix_provider)
      : _matrix_provider(std::move(matrix_provider))
    {}

    void setMatrixProvider(std::shared_ptr<MatrixProvider> matrix_provider)
    {
      if (_matrix_provider)
        _matrix_provider->unsubscribe(*this);
      _matrix_provider = std::move(matrix_provider);
      _matrix_provider->subscribe(*this);
    }

  private:

    std::shared_ptr<MatrixProvider> _matrix_provider;

  };


  template<typename Domain_, typename Range_>
  struct LinearSolver
    : public LinearSolverComponent<Domain_,Range_>
  {

    using Domain = Domain_;
    using Range  = Range_;
    using Real   = typename LinearSolverComponent<Domain_,Range_>::Real;

    virtual void solve(Domain&, Range& range) = 0;

    virtual void solve(Domain&, Range& range, Real reduction) = 0;

    virtual Real norm(const Range& range) const = 0;

  };

} // namespace Dune::PDELab::Backend

#endif // DUNE_PDELAB_BACKEND_COMMON_INTERFACE_H
