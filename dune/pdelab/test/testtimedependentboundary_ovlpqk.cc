// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/** \file
    \brief Solve heat equation with solution $t$ (time).
    Implemented after discovery of a bug in the overlapping scheme
    with time-dependent boundary conditions (March 2019).
*/
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream> // console writing
#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
// C includes
#include <sys/stat.h>
// dune-grid includes
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/io/file/gmshreader.hh>

#include <dune/pdelab.hh>

#include <vector>
#include <map>
#include <string>

#include <dune/common/filledarray.hh>
#include <dune/common/std/make_array.hh>

    /** \brief Parameter class for solving the linear heat equation with time dependent boundary
     *
     * A parameter class for the linear heat equation
     * \f{align*}{
     *     \partial_t u + \Delta u &=& f \mbox{ in } (\Omega,[0,T]),  \\
     *                           u &=& g \mbox{ on } (\partial\Omega_D,[0,T]), \\
     *            \nabla u \cdot n &=& 0 \mbox{ on } (\partial\Omega_N,[0,T]), \\
     *                           u &=& g \mbox{ in } (\Omega,t=0)
     * \f}
     * Note:
     *  - This class is a copy of convectiondiffusionparameter.hh with different functions f,g, and added time.
     *  - Functions $A,b,c,j,o,bctype$ are default, $g(x,t)=t$, and $f(x,t)=1$.
     *  - Analytic solution to this equation is $u=g(x,t)=t \mbox{in} (\Omega,[0,T])$.
     *
     * \tparam GV The GridView type
     * \tparam RF The range field type
     */
    template<typename GV, typename RF>
    class ConvectionDiffusionModelProblem
    {
      typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;

    public:
      typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

      //! tensor diffusion constant per cell? return false if you want more than one evaluation of A per cell.
      static constexpr bool permeabilityIsConstantPerCell()
      {
        return true;
      }

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
      f (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
      {
        return 1.0;
      }

      //! boundary condition type function
      BCType
      bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
      {
        return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
      }

      //! Dirichlet boundary condition value
      typename Traits::RangeFieldType
      g (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
      {
        return t;
      }

      //! Neumann boundary condition
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

      void setTime(double t_)
      {
        t=t_;
      }

    private:
      double t=0;
    };

//===============================================================
// Problem setup and solution
//===============================================================

// generate a P1 function and output it
template<typename GV, typename FEM>
void driver (const GV& gv, const FEM& fem)
{
  // constants and types
  typedef typename FEM::Traits::FiniteElementType::Traits::
    LocalBasisType::Traits::RangeFieldType RF;

  // make model problem
  typedef ConvectionDiffusionModelProblem<GV,RF> Problem;
  Problem problem;
  problem.setTime(0.);
  // initial condition
  auto glambda = [&](const auto& e, const auto& x){return problem.g(e,x);};
  auto gf = Dune::PDELab::makeInstationaryGridFunctionFromCallable(gv,glambda,problem);
  auto blambda = [&](const auto& i, const auto& x){return problem.bctype(i,x);};
  auto bf = Dune::PDELab::makeBoundaryConditionFromCallable(gv,blambda);
  using Dune::PDELab::Backend::native;

  // make function space
  typedef Dune::PDELab::OverlappingConformingDirichletConstraints CON;
  typedef Dune::PDELab::ISTL::VectorBackend<> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
  GFS gfs(gv,fem);
  gfs.name("timeDependentBoundary");

  // Coefficient vector
  using ZS = Dune::PDELab::Backend::Vector<GFS, RF>;
  ZS zs(gfs);
  Dune::PDELab::interpolate(gf,gfs,zs);
  // Assemble constraints
  typedef typename GFS::template ConstraintsContainer<RF>::Type CC;
  CC cc;
  Dune::PDELab::constraints(bf,gfs,cc); // assemble constraints
  // // skip text output
  // std::cout << "constrained dofs=" << cc.size() << " of "
  //           << gfs.globalSize() << std::endl;

  // make local operator
  typedef Dune::PDELab::ConvectionDiffusionFEM<Problem,FEM> LOP;
  LOP lop(problem);
  typedef Dune::PDELab::L2 TLOP;
  TLOP tlop(6);

  // make grid operator
  typedef Dune::PDELab::ISTL::BCRSMatrixBackend<> MBE;
  MBE mbe(27); // 27 is too large / correct for all test cases, so should work fine
  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,
                                     MBE,
                                     double,double,double,
                                     CC,CC> GO0;
  GO0 go0(gfs,cc,gfs,cc,lop,mbe);
  typedef Dune::PDELab::GridOperator<GFS,GFS,TLOP,
                                     MBE,
                                     double,double,double,
                                     CC,CC> GO1;
  GO1 go1(gfs,cc,gfs,cc,tlop,mbe);
  typedef Dune::PDELab::OneStepGridOperator<GO0,GO1> GO;
  GO go(go0,go1);

  // Select a linear solver backend  PARALLEL
  typedef Dune::PDELab::ISTLBackend_OVLP_CG_SSORk<GFS,CC> LS;
  int verbose=0;
  // if (gfs.gridView().comm().rank()==0) verbose=1; // no output necessary if test goes well
  LS ls(gfs,cc,100,5,verbose);

  typedef Dune::PDELab::StationaryLinearProblemSolver<GO,LS,ZS> PDESOLVER;
  PDESOLVER solver(go,ls,1e-10,1e-99,verbose);

  // select and prepare time-stepping scheme
  Dune::PDELab::OneStepThetaParameter<RF> method(1.0);
  typedef Dune::PDELab::OneStepMethod<RF,GO,PDESOLVER,ZS,ZS> OSM;
  OSM  osm(method,go,solver);
  osm.setVerbosityLevel(2*verbose);

  // // Does not produce VTK (very boring if correct). Lines in while loop commented out too.
  // // prepare VTK writer for time-dempendent problem and write first file
  // int subsampling=1.;
  // using VTKWRITER = Dune::SubsamplingVTKWriter<GV>;
  // VTKWRITER vtkwriter(gv,Dune::refinementIntervals(subsampling));
  // std::string filenamef="heatEq_OVLP_timeDependentBoundaries";
  // struct stat stf;
  // if( stat( filenamef.c_str(), &stf ) != 0 )
  // {
  //   int statf = 0;
  //   statf = mkdir(filenamef.c_str(),S_IRWXU|S_IRWXG|S_IRWXO);
  //   if( statf != 0 && statf != -1)
  //     std::cout << "Error: Cannot create directory "
  //               << filenamef << std::endl;
  // }
  // using VTKSEQUENCEWRITER = Dune::VTKSequenceWriter<GV>;
  // VTKSEQUENCEWRITER vtkSequenceWriter(
  //   std::make_shared<VTKWRITER>(vtkwriter),filenamef,filenamef,"");
  // // add data field for all components of the space to the VTK writer
  // Dune::PDELab::addSolutionToVTKWriter(vtkSequenceWriter,gfs,zs);
  // vtkSequenceWriter.write(0.0,Dune::VTK::appendedraw);

  ZS zsnew(zs);
  RF time = 0.;
  RF timestep = 0.01;
  RF T = 1.;
  while(time<T-1e-8)
  {
    try
    {
      problem.setTime(time+timestep);
      osm.apply(time,timestep,zs,gf,zsnew);

      zs=zsnew;
      time+=timestep;

      // // write solution
      // vtkSequenceWriter.write(time,Dune::VTK::appendedraw);

    } // end try
    catch(Dune::Exception& e)
    {
      std::cout << e.what() << std::endl;
      // vtkSequenceWriter.write(time+T,Dune::VTK::appendedraw);
      throw;
    }
    catch(...) // in case of unexpected exception terminate
    {
      throw;
    }
  } // end while

  // Error check:
  //   since the result is (spatially-wise) a constant function,
  //   no base function integration is necessary.
  //   All DoFs should have value "time", $u(x,t)=t$.
  RF error{0};
  for (const auto& v : zs)
  {
    RF vv = v-time;
    error += vv*vv;
  }
  if (error > 1E-18) // error is still squared, scale the tolerance accordingly
  {
    std::cout << "Too big error " << error << " in testtimedependentboundary_ovlpqk.cc" << std::endl;
    std::abort();
  }
  // std::cout << "On rank " << gfs.gridView().comm().rank() << " the error is: " << error << "." << std::endl;
} // end driver

//===============================================================
// Main program with grid setup
//===============================================================
int main(int argc, char** argv)
{
  // try{
    // initialize MPI, finalize is done automatically on exit
    Dune::MPIHelper::instance(argc,argv);
    if(Dune::MPIHelper::isFake)
    {
      std::cout<< "This is a sequential program. Designed for testing a parallel program with overlapping grid." << std::endl;
      std::abort();
    }

    // read ini file
    constexpr int refinement = 0;
    constexpr int degree = 1;
    constexpr int dim = 2;
    constexpr int overlap = 2;
    typedef Dune::YaspGrid<dim> Grid;
    typedef Grid::ctype DF;
    Dune::FieldVector<DF,dim> L;
    L[0] = 1.;
    L[1] = 1.;
    std::array<int,dim> N;
    N[0] = 32;
    N[1] = 32;
    std::bitset<dim> periodic(false);
    std::shared_ptr<Grid> gridp = std::shared_ptr<Grid>(new Grid(L,N,periodic,overlap));
    gridp->refineOptions(false); // keep overlap in cells
    gridp->globalRefine(refinement);
    typedef Grid::LeafGridView GV;
    GV gv=gridp->leafGridView();
    typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,double,degree> FEM;
    FEM fem(gv);
    driver(gv,fem);
  // }
  // catch (Dune::Exception &e)
  // {
  //   std::cerr << "Dune reported error: " << e << std::endl;
  //   return 1;
  // }
  // catch (...)
  // {
  //   std::cerr << "Unknown exception thrown!" << std::endl;
  //   return 1;
  // }
}
