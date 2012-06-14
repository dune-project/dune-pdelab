// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include<iostream>
#include<vector>
#include<map>
#include<string>
#include<dune/common/mpihelper.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include <dune/common/shared_ptr.hh>
#include<dune/common/static_assert.hh>
#include<dune/grid/yaspgrid.hh>
#include<dune/istl/bvector.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/preconditioners.hh>
#include<dune/istl/io.hh>
#include<dune/istl/paamg/amg.hh>

#include"../finiteelementmap/p0fem.hh"
#include"../finiteelementmap/p12dfem.hh"
#include"../finiteelementmap/pk2dfem.hh"
#include"../finiteelementmap/q12dfem.hh"
#include"../finiteelementmap/q22dfem.hh"
#include"../finiteelementmap/q1fem.hh"
#include"../finiteelementmap/conformingconstraints.hh"
#include"../gridfunctionspace/gridfunctionspace.hh"
#include"../gridfunctionspace/gridfunctionspaceutilities.hh"
#include"../gridfunctionspace/interpolate.hh"
#include"../constraints/constraints.hh"
#include"../common/function.hh"
#include"../common/vtkexport.hh"
#include"../backend/istlvectorbackend.hh"
#include"../backend/istlmatrixbackend.hh"
#include"../gridoperator/gridoperator.hh"
#include"../backend/seqistlsolverbackend.hh"
#include"../localoperator/laplacedirichletp12d.hh"
#include"../localoperator/poisson.hh"

#include <dune/pdelab/common/benchmarkhelper.hh>
#include <dune/common/parametertreeparser.hh>

#include"../test/gridexamples.hh"

//===============================================================
//===============================================================
// Solve the Poisson equation
//           - \Delta u = f in \Omega,
//                    u = g on \partial\Omega_D
//  -\nabla u \cdot \nu = j on \partial\Omega_N
//===============================================================
//===============================================================

//===============================================================
// Define parameter functions f,g,j and \partial\Omega_D/N
//===============================================================

// function for defining the source term
template<typename GV, typename RF>
class F
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
                                                  F<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,F<GV,RF> > BaseT;

  F (const GV& gv) : BaseT(gv) {}
  inline void evaluateGlobal (const typename Traits::DomainType& x,
                              typename Traits::RangeType& y) const
  {
    if (x[0]>0.25 && x[0]<0.375 && x[1]>0.25 && x[1]<0.375)
      y = 50.0;
    else
      y = 0.0;
  }
};


// boundary grid function selecting boundary conditions
class ConstraintsParameters
  : public Dune::PDELab::DirichletConstraintsParameters
{

public:

  template<typename I>
  bool isDirichlet(const I & ig, const Dune::FieldVector<typename I::ctype, I::dimension-1> & x) const
  {
    Dune::FieldVector<typename I::ctype,I::dimension>
      xg = ig.geometry().global(x);

    if (xg[1]<1E-6 || xg[1]>1.0-1E-6)
      {
        return false;
      }
    if (xg[0]>1.0-1E-6 && xg[1]>0.5+1E-6)
      {
        return false;
      }
    return true;
  }

};


// function for Dirichlet boundary conditions and initialization
template<typename GV, typename RF>
class G
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
                                                  G<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,G<GV,RF> > BaseT;

  G (const GV& gv) : BaseT(gv) {}
  inline void evaluateGlobal (const typename Traits::DomainType& x,
                              typename Traits::RangeType& y) const
  {
    typename Traits::DomainType center;
    for (int i=0; i<GV::dimension; i++) center[i] = 0.5;
    center -= x;
    y = exp(-center.two_norm2());
  }
};

// function for defining the flux boundary condition
template<typename GV, typename RF>
class J
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
                                                  J<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,J<GV,RF> > BaseT;

  J (const GV& gv) : BaseT(gv) {}
  inline void evaluateGlobal (const typename Traits::DomainType& x,
                              typename Traits::RangeType& y) const
  {
    if (x[1]<1E-6 || x[1]>1.0-1E-6)
      {
        y = 0;
        return;
      }
    if (x[0]>1.0-1E-6 && x[1]>0.5+1E-6)
      {
        y = -5.0;
        return;
      }
  }
};

//===============================================================
// Problem setup and solution
//===============================================================

// generate a P1 function and output it
template<typename GV, typename FEM, typename CON, int q>
void poisson (const GV& gv, const FEM& fem, std::string filename, const bool solve, std::size_t runs)
{

  Dune::PDELab::BenchmarkHelper<> bh(filename,runs);

  for (std::size_t run = 0; run < runs; ++run)
    {
      bh.start_run(std::cout);
      // constants and types
      typedef typename GV::Grid::ctype DF;
      typedef typename FEM::Traits::FiniteElementType::Traits::
        LocalBasisType::Traits::RangeFieldType R;

      bh.start("global setup",std::cout);
      bh.start("GFS setup",std::cout);

      // make function space
      typedef Dune::PDELab::GridFunctionSpace<
        GV,
        FEM,
        CON,
        Dune::PDELab::ISTLVectorBackend<1>
        > GFS;
      GFS gfs(gv,fem);

      bh.end("GFS setup",std::cout);

      // include timing for ordering update to make
      // listings compatible with new infrastructure
      bh.start("ordering update",std::cout);
      bh.end("ordering update",std::cout);

      bh.end("global setup",std::cout);
      bh.start("constraints",std::cout);

      // make constraints map and initialize it from a function
      typedef typename GFS::template ConstraintsContainer<R>::Type C;
      C cg;
      cg.clear();
      ConstraintsParameters constraintsparameters;
      Dune::PDELab::constraints(constraintsparameters,gfs,cg);

      bh.end("constraints",std::cout);
      bh.start("LOP construction",std::cout);

      // make local operator
      typedef G<GV,R> GType;
      GType g(gv);
      typedef F<GV,R> FType;
      FType f(gv);
      typedef J<GV,R> JType;
      JType j(gv);
      typedef Dune::PDELab::Poisson<FType,ConstraintsParameters,JType,q> LOP;
      LOP lop(f,constraintsparameters,j);

      bh.end("LOP construction",std::cout);
      bh.start("GOP construction",std::cout);

      // make grid operator
      typedef Dune::PDELab::GridOperator<
        GFS,GFS,LOP,
        Dune::PDELab::ISTLBCRSMatrixBackend<1,1>,
        double,double,double,
        C,C> GridOperator;
      GridOperator gridoperator(gfs,cg,gfs,cg,lop);

      bh.end("GOP construction",std::cout);
      bh.start("vector creation",std::cout);

      // make coefficent Vector and initialize it from a function
      typedef typename GridOperator::Traits::Domain DV;

      DV x0(gfs);
      x0 = 0.0;

      bh.end("vector creation",std::cout);
      bh.start("interpolation",std::cout);

      Dune::PDELab::interpolate(g,gfs,x0);

      bh.end("interpolation",std::cout);
      bh.start("set free DOFs",std::cout);

      Dune::PDELab::set_nonconstrained_dofs(cg,0.0,x0);

      bh.end("set free DOFs",std::cout);
      bh.start("matrix creation",std::cout);

      // represent operator as a matrix
      typedef typename GridOperator::Traits::Jacobian M;
      M m(gridoperator);
      m = 0.0;

      bh.end("matrix creation",std::cout);
      bh.start("jacobian",std::cout);

      gridoperator.jacobian(x0,m);
      //  Dune::printmatrix(std::cout,m.base(),"global stiffness matrix","row",9,1);

      bh.end("jacobian",std::cout);

      // evaluate residual w.r.t initial guess
      typedef typename GridOperator::Traits::Range RV;
      RV r(gfs);
      r = 0.0;

      bh.start("residual",std::cout);

      gridoperator.residual(x0,r);

      bh.end("residual",std::cout);
      bh.start("solve",std::cout);

      DV x(gfs,0.0);

      if (solve)
        {

          // make ISTL solver
          Dune::MatrixAdapter<M,DV,RV> opa(m);
          typedef Dune::PDELab::OnTheFlyOperator<DV,RV,GridOperator> ISTLOnTheFlyOperator;
          //ISTLOnTheFlyOperator opb(gridoperator);
          Dune::SeqSSOR<M,DV,RV> ssor(m,1,1.0);
          Dune::SeqILU0<M,DV,RV> ilu0(m,1.0);
          Dune::Richardson<DV,RV> richardson(1.0);

          //   typedef Dune::Amg::CoarsenCriterion<Dune::Amg::SymmetricCriterion<M,
          //     Dune::Amg::FirstDiagonal> > Criterion;
          //   typedef Dune::SeqSSOR<M,V,V> Smoother;
          //   typedef typename Dune::Amg::SmootherTraits<Smoother>::Arguments SmootherArgs;
          //   SmootherArgs smootherArgs;
          //   smootherArgs.iterations = 2;
          //   int maxlevel = 20, coarsenTarget = 100;
          //   Criterion criterion(maxlevel, coarsenTarget);
          //   criterion.setMaxDistance(2);
          //   typedef Dune::Amg::AMG<Dune::MatrixAdapter<M,V,V>,V,Smoother> AMG;
          //   AMG amg(opa,criterion,smootherArgs,1,1);

          Dune::CGSolver<DV> solvera(opa,ilu0,1E-10,5000,1);
          // FIXME: Use ISTLOnTheFlyOperator in the second solver again
          Dune::CGSolver<DV> solverb(opa,richardson,1E-10,5000,1);
          Dune::InverseOperatorResult stat;

          // solve the jacobian system
          r *= -1.0; // need -residual
          solvera.apply(x,r,stat);
          x += x0;
        }

      bh.end("solve",std::cout);
      bh.start("IO",std::cout);

      // make discrete function object
      typedef Dune::PDELab::DiscreteGridFunction<GFS,DV> DGF;
      DGF dgf(gfs,x);

      // output grid function with VTKWriter
      Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
      vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(dgf,"solution"));
      vtkwriter.write(filename,Dune::VTK::ascii);

      bh.end("IO",std::cout);

      bh.end_run(std::cout);

    }
  bh.print(std::cout);

  std::stringstream timings_name;
  timings_name << filename << "_timings.txt";

  std::ofstream timings_file(timings_name.str());
  timings_file << filename << " " << runs << std::endl;
  bh.print(timings_file);
}

template<typename Grid>
std::string name(std::string grid_name, const Grid& grid, std::string method_name, int p, int q, bool solve)
{
  std::stringstream n;
  n << "poisson-old_" << method_name << p
    << "_" << grid_name
    << "_" << Grid::dimension << "D"
    << "_l" << grid.maxLevel()
    << "_q" << q
    << (solve ? "_solve" : "_nosolve");
  return n.str();
}

//===============================================================
// Main program with grid setup
//===============================================================

int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper::instance(argc, argv);

    Dune::ParameterTree params;

    if (argc == 1)
      {
        std::cerr << "no parameter file passed, defaulting to poisson-old.ini..." << std::endl;
        Dune::ParameterTreeParser::readINITree("poisson-old.ini",params);
      }
    else if (argc == 2)
      {
        std::cerr << "reading parameters from " << argv[1] << "..." << std::endl;
        Dune::ParameterTreeParser::readINITree(argv[1],params);
      }
    else
      {
        std::cerr << "Usage: " << argv[0] << " [parameter file]" << std::endl;
        return 64;
      }

    const std::size_t global_runs = params.get("global.runs",5);
    const bool global_solve = params.get("global.solve",false);

    // YaspGrid Q1 2D test
    if (params.hasSub("Yasp_Q1_2D") && params.get<bool>("Yasp_Q1_2D.enabled"))
      {
        const Dune::ParameterTree& lp = params.sub("Yasp_Q1_2D");

        const int refine = lp.get<int>("refine");
        const int runs = lp.get("runs",global_runs);
        const bool solve = lp.get("solve",global_solve);

        // make grid
        Dune::FieldVector<double,2> L(1.0);
        Dune::FieldVector<int,2> N(1);
        Dune::FieldVector<bool,2> B(false);
        Dune::YaspGrid<2> grid(L,N,B,0);
        grid.globalRefine(refine);

        // get view
        typedef Dune::YaspGrid<2>::LeafGridView GV;
        const GV& gv=grid.leafView();

        // make finite element map
        typedef GV::Grid::ctype DF;
        typedef Dune::PDELab::Q12DLocalFiniteElementMap<DF,double> FEM;
        FEM fem;

        // solve problem
        poisson<GV,FEM,Dune::PDELab::ConformingDirichletConstraints,2>(gv,fem,name("Yasp",grid,"Q",1,2,solve),solve,runs);
      }

    // YaspGrid Q2 2D test
    if (params.hasSub("Yasp_Q2_2D") && params.get<bool>("Yasp_Q2_2D.enabled"))
      {
        const Dune::ParameterTree& lp = params.sub("Yasp_Q2_2D");

        const int refine = lp.get<int>("refine");
        const int runs = lp.get("runs",global_runs);
        const bool solve = lp.get("solve",global_solve);

        // make grid
        Dune::FieldVector<double,2> L(1.0);
        Dune::FieldVector<int,2> N(1);
        Dune::FieldVector<bool,2> B(false);
        Dune::YaspGrid<2> grid(L,N,B,0);
        grid.globalRefine(refine);

        // get view
        typedef Dune::YaspGrid<2>::LeafGridView GV;
        const GV& gv=grid.leafView();

        // make finite element map
        typedef GV::Grid::ctype DF;
        typedef Dune::PDELab::Q22DLocalFiniteElementMap<DF,double> FEM;
        FEM fem;

        // solve problem
        poisson<GV,FEM,Dune::PDELab::ConformingDirichletConstraints,4>(gv,fem,name("Yasp",grid,"Q",2,2,solve),solve,runs);
      }

    // YaspGrid Q2 3D test
    if (params.hasSub("Yasp_Q2_3D") && params.get<bool>("Yasp_Q2_3D.enabled"))
      {
        const Dune::ParameterTree& lp = params.sub("Yasp_Q2_3D");

        const int refine = lp.get<int>("refine");
        const int runs = lp.get("runs",global_runs);
        const bool solve = lp.get("solve",global_solve);

        // make grid
        Dune::FieldVector<double,3> L(1.0);
        Dune::FieldVector<int,3> N(1);
        Dune::FieldVector<bool,3> B(false);
        Dune::YaspGrid<3> grid(L,N,B,0);
        grid.globalRefine(refine);

        // get view
        typedef Dune::YaspGrid<3>::LeafGridView GV;
        const GV& gv=grid.leafView();

        // make finite element map
        typedef GV::Grid::ctype DF;
        typedef Dune::PDELab::Q1LocalFiniteElementMap<DF,double,3> FEM;
        FEM fem;

        // solve problem
        poisson<GV,FEM,Dune::PDELab::ConformingDirichletConstraints,4>(gv,fem,name("Yasp",grid,"Q",2,2,solve),solve,runs);
      }

    // UG Pk 2D test
#if HAVE_UG
    if (params.hasSub("UG_Pk_2D") && params.get<bool>("UG_Pk_2D.enabled"))
      {
        const Dune::ParameterTree& lp = params.sub("UG_Pk_2D");

        const int refine = lp.get<int>("refine");
        const int runs = lp.get("runs",global_runs);
        const bool solve = lp.get("solve",global_solve);

        // make grid
        Dune::shared_ptr<Dune::UGGrid<2> > grid(TriangulatedUnitSquareMaker<Dune::UGGrid<2> >::create());
        grid->globalRefine(refine);

        // get view
        typedef Dune::UGGrid<2>::LeafGridView GV;
        const GV& gv=grid->leafView();

        // make finite element map
        typedef GV::Grid::ctype DF;
        typedef double R;
        const int k=3;
        const int q=2*k;
        typedef Dune::PDELab::Pk2DLocalFiniteElementMap<GV,DF,double,k> FEM;
        FEM fem(gv);

        // solve problem
        poisson<GV,FEM,Dune::PDELab::ConformingDirichletConstraints,q>(gv,fem,name("UG",*grid,"P",k,q,solve),solve,runs);

      }
#endif

#if HAVE_ALBERTA
    if (params.hasSub("Alberta_Pk_2D") && params.get<bool>("Alberta_Pk_2D.enabled"))
      {
        const Dune::ParameterTree& lp = params.sub("Alberta_Pk_2D");

        const int refine = lp.get<int>("refine");
        const int runs = lp.get("runs",global_runs);
        const bool solve = lp.get("solve",global_solve);

        // make grid
        AlbertaUnitSquare grid;
        grid.globalRefine(refine);

        // get view
        typedef AlbertaUnitSquare::LeafGridView GV;
        const GV& gv=grid.leafView();

        // make finite element map
        typedef GV::Grid::ctype DF;
        typedef double R;
        const int k=3;
        const int q=2*k;
        typedef Dune::PDELab::Pk2DLocalFiniteElementMap<GV,DF,double,k> FEM;
        FEM fem(gv);

        // solve problem
        poisson<GV,FEM,Dune::PDELab::ConformingDirichletConstraints,q>(gv,fem,name("Alberta",grid,"P",k,q,solve),solve,runs);
      }
#endif

#if HAVE_ALUGRID
    if (params.hasSub("ALU_Pk_2D") && params.get<bool>("ALU_Pk_2D.enabled"))
      {
        const Dune::ParameterTree& lp = params.sub("ALU_Pk_2D");

        const int refine = lp.get<int>("refine");
        const int runs = lp.get("runs",global_runs);
        const bool solve = lp.get("solve",global_solve);

        // make grid
        ALUUnitSquare grid;
        grid.globalRefine(refine);

        // get view
        typedef ALUUnitSquare::LeafGridView GV;
        const GV& gv=grid.leafView();

        // make finite element map
        typedef GV::Grid::ctype DF;
        typedef double R;
        const int k=3;
        const int q=2*k;
        typedef Dune::PDELab::Pk2DLocalFiniteElementMap<GV,DF,double,k> FEM;
        FEM fem(gv);

        // solve problem
        poisson<GV,FEM,Dune::PDELab::ConformingDirichletConstraints,q>(gv,fem,name("ALU",grid,"P",k,q,solve),solve,runs);
      }
#endif

    // test passed
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
