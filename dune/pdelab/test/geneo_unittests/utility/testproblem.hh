#ifndef TESTPROBLEM_HH
#define TESTPROBLEM_HH

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/pdelab.hh>


namespace Utility
{

  // scoped enum indicating which boundary conditions to set up in the test
  // problem. The values are set in the respective class itself.
  //
  //   - BC::onlyDirichlet: The entire boundary satisfies a Dirichlet condition
  //   - BC::onlyNeumann: The entire boundary satisfies a Neumann condition
  //   - BC::DirichletAndNeumann: The upper and lower domain boundaries satisfy
  //     a Dirichlet condition and the left and right boundaries satisfy a
  //     Neumann condition
  enum class BC
  {
    onlyDirichlet,
    onlyNeumann,
    DirichletAndNeumann
  };


  // generic Poisson problem with homogeneous coefficient
  template<class GV, typename DF, typename RF>
  class PoissonTestProblem
  {
  private:

    using BCType = Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type;

    BC bc_;
    double bceps_; // threshold value for the distance to the boundary
    Dune::FieldVector<DF, GV::Grid::dimension> domainsize_;

  protected:

    RF coeffValue_;


  public:

    using Traits = Dune::PDELab::ConvectionDiffusionParameterTraits<GV, RF>;

    static constexpr bool permeabilityIsConstantPerCell() { return true; }


    PoissonTestProblem()
      : bc_{ BC::onlyDirichlet }
      , coeffValue_{ 1.0 }
      , bceps_{ 1e-6 }
      , domainsize_{ Dune::FieldVector<DF, GV::Grid::dimension>(1.0) }
    {}


    PoissonTestProblem(const GV& gv, BC bc, RF coeffValue)
      : bc_{ bc }
      , coeffValue_{ coeffValue }
      , bceps_{ 1e-6 }
      , domainsize_{ gv.grid().domainSize() }
    {}


    // permeability tensor
    typename Traits::PermTensorType
    A(const typename Traits::ElementType& e,
      const typename Traits::DomainType& x) const
    {
      typename Traits::PermTensorType I{{ coeffValue_, 0.0        },
                                        { 0.0        , coeffValue_}};
      return I;
    }


    // velocity field
    typename Traits::RangeType
    b(const typename Traits::ElementType& e,
      const typename Traits::DomainType& x) const
    {
      typename Traits::RangeType v(0.0);
      return v;
    }


    // sink term
    typename Traits::RangeFieldType
    c(const typename Traits::ElementType& e,
      const typename Traits::DomainType& x) const
    {
      return 0.0;
    }


    // source term
    typename Traits::RangeFieldType
    f(const typename Traits::ElementType& e,
      const typename Traits::DomainType& x) const
    {
      return 0.0;
    }


    // boundary condition type
    BCType bctype(const typename Traits::IntersectionType& is,
                  const typename Traits::IntersectionDomainType& x) const
    {
      switch (bc_)
      {
      case BC::onlyDirichlet :
        return BCType::Dirichlet;

      case BC::onlyNeumann :
        return BCType::Neumann;

      case BC::DirichletAndNeumann :
      {
        typename Traits::DomainType xglobal{ is.geometry().global(x) };

        if(xglobal[0] < bceps_ || xglobal[0] > domainsize_[0] - bceps_)
          return BCType::Neumann;
        else
          return BCType::Dirichlet;
      }
      default :
        DUNE_THROW(Dune::Exception,
                   "Unknown boundary condition type in testproblem.");
      }
    }


    // Dirichlet boundary condition value
    typename Traits::RangeFieldType
    g(const typename Traits::ElementType& e,
      const typename Traits::DomainType& x) const
    {
      typename Traits::DomainType xglobal{ e.geometry().global(x) };

      // 1.0 on the lower boundary, 0.0 everywhere else
      if (xglobal[1] < bceps_)
        return 1.0;
      else
        return 0.0;
    }


    // Neumann boundary condition
    typename Traits::RangeFieldType
    j(const typename Traits::IntersectionType& is,
      const typename Traits::IntersectionDomainType& x) const
    {
      typename Traits::DomainType xglobal{ is.geometry().global(x) };

      // 1.0 on the left boundary, 0.0 everywhere else
      if (xglobal[0] < bceps_)
        return 1.0;
      else
        return 0.0;
    }


    // outflow boundary condition
    typename Traits::RangeFieldType
    o(const typename Traits::IntersectionType& is,
      const typename Traits::IntersectionDomainType& x) const
    {
      return 0.0;
    }
  };


  // Poisson problem, where the coefficient is zero in outside the overlap
  // region. We do not need virtual functions here, as we only instantiate the
  // problems themselves once and virtual functions would result in
  // unneccessary performance overhead
  // TODO: Currently zero outside overlap for this subdomain. Needs to include
  //       the overlap of the other subdomains as well.
  template<class GV, typename DF, typename RF>
  class PoissonZeroInterior : public PoissonTestProblem<GV, DF, RF>
  {
  public:

    // use Traits of the parent class
    using Traits = typename PoissonTestProblem<GV, DF, RF>::Traits;

    PoissonZeroInterior()
      : PoissonTestProblem<GV, DF, RF>()
    {}


    PoissonZeroInterior(const GV& gv, BC bc, RF coeffValue)
      : PoissonTestProblem<GV, DF, RF>(gv, bc, coeffValue)
    {}


    // permeability tensor
    typename Traits::PermTensorType
    A(const typename Traits::ElementType& e,
      const typename Traits::DomainType& x) const
    {
      using PT = Dune::PartitionType;

      // zero for interior entities and the value from the testproblem for all
      // other entities
      switch(e.partitionType())
      {
      case PT::OverlapEntity :
      {
        typename Traits::PermTensorType I{{ 0.0, 0.0 },
                                          { 0.0, 0.0 }};
        return I;
      }
      default :
      {
        RF cVal{ PoissonTestProblem<GV, DF, RF>::coeffValue_ };
        typename Traits::PermTensorType I{{ cVal, 0.0 },
                                          { 0.0 , cVal}};
        return I;
      }
      }
    }
  };

} // namespace Utility

#endif
