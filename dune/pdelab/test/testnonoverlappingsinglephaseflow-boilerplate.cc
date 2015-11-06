// -*- tab-width: 4; indent-tabs-mode: nil -*-
/** \file
    \brief Solve single phase flow problem in parallel on non-overlapping grids using conforming linear finite elements with boilerplate code
*/
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/pdelab/boilerplate/pdelab.hh>
#include <dune/pdelab/localoperator/convectiondiffusionfem.hh>

#include<dune/pdelab/gridfunctionspace/vtk.hh>

#include "testnonoverlappingsinglephaseflow-boilerplate-problem.hh"

int main(int argc, char **argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
    if(Dune::MPIHelper::isFake)
      std::cout<< "This is a sequential program." << std::endl;
    else
    {
      if(helper.rank()==0)
        std::cout << "parallel run on " << helper.size() << " process(es)" << std::endl;
    }

    // define parameters
    const unsigned int dim = 3;
    const unsigned int degree = 1;
    const Dune::PDELab::MeshType meshtype = Dune::PDELab::MeshType::conforming;
    const Dune::GeometryType::BasicType elemtype = Dune::GeometryType::simplex;
    const Dune::SolverCategory::Category solvertype = Dune::SolverCategory::nonoverlapping;
    typedef double NumberType;

    // make grid
    typedef Dune::UGGrid<dim> GM;
    typedef Dune::PDELab::UnstructuredGrid<GM> Grid;
    Grid grid(GRIDSDIR "/cube1045.msh",true,true);
    grid->loadBalance();

    std::cout << " after load balance /" << helper.rank() << "/ " << grid->size(0) << std::endl;

    // get view
    typedef GM::LeafGridView GV;
    using ES = Dune::PDELab::NonOverlappingEntitySet<GV>;
    ES es(grid->leafGridView());

    // make problem parameters
    typedef ParameterA<ES,NumberType> Problem;
    Problem problem(es);
    typedef Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<Problem> BCType;
    BCType bctype(es,problem);

    // make a finite element space
    typedef Dune::PDELab::CGSpace<GM,NumberType,degree,BCType,elemtype,meshtype,solvertype> FS;
    FS fs(*grid,bctype);

    // make a degree of freedom vector and initialize it with a function
    typedef FS::DOF X;
    X x(fs.getGFS(),0.0);
    typedef Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<Problem> G;
    G g(es,problem);
    Dune::PDELab::interpolate(g,fs.getGFS(),x);

    // assemble constraints
    fs.assembleConstraints(bctype);
    fs.setNonConstrainedDOFS(x,0.0);

    // assembler for finite elemenent problem
    typedef Dune::PDELab::ConvectionDiffusionFEM<Problem,typename FS::FEM> LOP;
    LOP lop(problem);
    typedef Dune::PDELab::GalerkinGlobalAssembler<FS,LOP,solvertype> ASSEMBLER;
    ASSEMBLER assembler(fs,lop,std::pow(2*degree+1,dim));

    // make linear solver and solve problem
    typedef Dune::PDELab::ISTLSolverBackend_IterativeDefault<FS,ASSEMBLER,solvertype> SBE;
    SBE sbe(fs,assembler,5000,1);
    typedef Dune::PDELab::StationaryLinearProblemSolver<ASSEMBLER::GO,SBE::LS,X> SLP;
    SLP slp(*assembler,*sbe,x,1e-6);
    slp.apply();

    // output grid to VTK file
    typedef typename ES::Traits::GridView GV;
    Dune::SubsamplingVTKWriter<GV> vtkwriter(es.gridView(),3);
    Dune::PDELab::addSolutionToVTKWriter(vtkwriter,fs.getGFS(),x);
    vtkwriter.write("advection_stationary",Dune::VTK::ascii);

    // done
    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
    return 1;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
    return 1;
  }
}
