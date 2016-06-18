// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_TEST_TESTPOISSON_PERIODIC_HH
#define DUNE_PDELAB_TEST_TESTPOISSON_PERIODIC_HH


#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/pdelab/boilerplate/pdelab.hh>
#include <dune/pdelab/localoperator/convectiondiffusionfem.hh>
#include <dune/pdelab/localoperator/convectiondiffusiondg.hh>

// Poisson problem definition
template<typename GV, typename RF>
class GenericEllipticProblem
{
  typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;

public:
  typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

  //! tensor diffusion coefficient
  typename Traits::PermTensorType
  A (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::PermTensorType I;
    for (std::size_t i=0; i<Traits::dimDomain; i++)
      for (std::size_t j=0; j<Traits::dimDomain; j++)
        I[i][j] = (i==j) ? 1 : 0;
    return I;
  }

  //! velocity field
  typename Traits::RangeType
  b (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::RangeType v(0.0);
    return v;
  }

  //! sink term
  typename Traits::RangeFieldType
  c (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return 0.0;
  }

  //! source term
  typename Traits::RangeFieldType
  f (const typename Traits::ElementType& e, const typename Traits::DomainType& xlocal) const
  {
    auto x = e.geometry().global(xlocal);
    return x[0] * std::sin(5.0*M_PI*x[1]) + std::exp(-((x[0]-0.5)*(x[0]-0.5) + (x[1]-0.5)*(x[1]-0.5)) / 0.02);
  }

  //! Boundary condition type function. Will not be evaluated on periodic boundary, so we simply set Dirichlet.
  BCType
  bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
  }

  //! Dirichlet boundary condition value
  typename Traits::RangeFieldType
  g (const typename Traits::ElementType& e, const typename Traits::DomainType& xlocal) const
  {
    return 0.0;
  }

  //! flux boundary condition
  typename Traits::RangeFieldType
  j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return 0.0;
  }

  //! outflow boundary condition
  typename Traits::RangeFieldType
  o (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return 0.0;
  }
};

// Problem solver
template<typename GridType, typename NumberType, int dim>
void poisson (GridType& grid)
{
#if DG == 0
    const Dune::PDELab::MeshType meshtype = Dune::PDELab::MeshType::conforming;
#endif
    const Dune::SolverCategory::Category solvertype = Dune::SolverCategory::overlapping;
    const Dune::GeometryType::BasicType elemtype = Dune::GeometryType::cube;

    // make problem parameters
    typedef GenericEllipticProblem<typename GridType::LeafGridView,NumberType> Problem;
    Problem problem;
    typedef Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<Problem> BCType;
    BCType bctype(grid.leafGridView(),problem);

    // make a finite element space
#if DG == 0
    typedef Dune::PDELab::CGSpace<GridType,NumberType,DEGREE,BCType,elemtype,meshtype,solvertype> FS;
    FS fs(grid,bctype);
#else
    typedef Dune::PDELab::DGQkSpace<GridType,NumberType,DEGREE,elemtype,solvertype> FS;
    FS fs(grid.leafGridView());
#endif

    // make a degree of freedom vector and initialize it with a function
    typedef typename FS::DOF V;
    V x(fs.getGFS(),0.0);

#if DG == 0
    typedef Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<Problem> G;
    G g(grid.leafGridView(),problem);
    Dune::PDELab::interpolate(g,fs.getGFS(),x);
#endif

    // assemble constraints
    fs.assembleConstraints(bctype);
    fs.setNonConstrainedDOFS(x,0.0);

    // assembler for finite elemenent problem
#if DG == 0
    typedef Dune::PDELab::ConvectionDiffusionFEM<Problem,typename FS::FEM> LOP;
    LOP lop(problem);
    typedef Dune::PDELab::GalerkinGlobalAssembler<FS,LOP,solvertype> ASSEMBLER;
    ASSEMBLER assembler(fs,lop,27);
#else
    typedef Dune::PDELab::ConvectionDiffusionDG<Problem,typename FS::FEM> LOP;
    LOP lop(problem,Dune::PDELab::ConvectionDiffusionDGMethod::SIPG,Dune::PDELab::ConvectionDiffusionDGWeights::weightsOn,2.0);
    typedef Dune::PDELab::GalerkinGlobalAssembler<FS,LOP,solvertype> ASSEMBLER;
    ASSEMBLER assembler(fs,lop,27);
#endif

    // make linear solver and solve problem
    typedef Dune::PDELab::ISTLSolverBackend_IterativeDefault<FS,ASSEMBLER,solvertype> SBE;
    SBE sbe(fs,assembler,5000,1);

    typedef Dune::PDELab::StationaryLinearProblemSolver<typename ASSEMBLER::GO,typename SBE::LS,V> SLP;
    SLP slp(*assembler,*sbe,x,1e-6);
    slp.apply();

    // // print statistics about nonzero values per row
    // typename ASSEMBLER::MAT m(assembler.getGO());
    // std::cout << m.patternStatistics() << std::endl;

    // output grid to VTK file
    Dune::SubsamplingVTKWriter<typename GridType::LeafGridView> vtkwriter(grid.leafGridView(),DEGREE-1);
    typename FS::DGF xdgf(fs.getGFS(),x);
    vtkwriter.addVertexData(std::make_shared<typename FS::VTKF>(xdgf,"x_h"));
    auto out_name = "poisson_periodic_" + std::to_string(dim) + "d_q" + std::to_string(DEGREE) + "_dg" + std::to_string(DG);
    vtkwriter.write(out_name, Dune::VTK::appendedraw);
}

#endif // DUNE_PDELAB_TEST_TESTPOISSON_PERIODIC_HH
