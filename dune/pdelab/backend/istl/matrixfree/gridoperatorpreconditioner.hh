// -*- tab-width: 2; indent-tabs-mode: nil -*-
// vi: set et ts=2 sw=2 sts=2:

#ifndef DUNE_PDELAB_BACKEND_ISTL_MATRIXFREE_GRIDOPERATORPRECONDITIONER_HH
#define DUNE_PDELAB_BACKEND_ISTL_MATRIXFREE_GRIDOPERATORPRECONDITIONER_HH

#include<dune/grid/common/rangegenerators.hh>
#include<dune/istl/preconditioner.hh>
#include<dune/istl/solvercategory.hh>

namespace Dune::PDELab
{

  /** \brief Turn a grid operator that represents a preconditioner solve into a ISTL preconditioner

      \tparam PrecGO Grid operator implementing the matrix-free preconditioner application
  */
  template<class PrecGO>
  class GridOperatorPreconditioner
    : public Dune::Preconditioner<typename PrecGO::Traits::Domain,typename PrecGO::Traits::Range>
  {
    using Domain = typename PrecGO::Traits::Domain;
    using Range = typename PrecGO::Traits::Range;

  public :
    static constexpr bool isLinear = PrecGO::LocalAssembler::isLinear();
    Dune::SolverCategory::Category category() const override
    {
      return Dune::SolverCategory::sequential;
    }

    GridOperatorPreconditioner(const PrecGO& precgo)
      : _precgo(precgo)
      , _u(nullptr)
    {}

    //! Set linearization point
    //! Must be called before apply() for nonlinear problems.
    void setLinearizationPoint(const Domain& u)
    {
      _u = &u;
    }

    //! prepare tensor product preconditioner if desired
    void pre(Domain& v, Range& d) override
    {
      if (not isLinear and _u == nullptr)
        DUNE_THROW(Dune::InvalidStateException, "You seem to apply a preconditioner to a linearized problem but haven't set the linearization point explicitly!");

      auto& lop = _precgo.localAssembler().localOperator();
      if (lop.requireSetup()){
        // We use the residual methods of the preconditioner grid operator to
        // do some set up work, if required (this could eg be the assembly and
        // inversion of the block diagonal). This was done through this
        // interface since the dimension match nicely and it avoids rewriting
        // assembler code like iterating over the grid, binding local function
        // spaces, etc.
        //
        // In the linear case the Jacobian does not depend on the current
        // solution (*_u in the nonlinear case). This means we can just pass a
        // dummy object of the right type, since it won't be used in the local
        // operator.
        //
        // Note:
        // - We do not have the current solution in the linear case so
        //   passing u is not possible
        // - Having access to the current solution here (for the linear case) is not
        //   difficult but needs changes in the solver backends and the linear problem
        //   solver.
        if (isLinear)
          _precgo.residual(d, d);
        else
          _precgo.residual(*_u, d);
        lop.setRequireSetup(false);
      }
    }

    void apply(Domain& v, const Range& d) override
    {
      if(isLinear)
        _precgo.jacobian_apply(d, v);
      else
        _precgo.jacobian_apply(*_u, d, v);
    }

    void post(Domain& v) override {}

  private :
    const PrecGO& _precgo;
    const Domain* _u;
  }; // end class GridOperatorPreconditioner

} // namespace Dune::PDELab
#endif
