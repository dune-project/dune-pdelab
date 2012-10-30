// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include <vector>
#include <map>
#include <sstream>
#include <fstream>

#include<dune/common/mpihelper.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/static_assert.hh>
#include<dune/grid/yaspgrid.hh>
#include<dune/istl/bvector.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/preconditioners.hh>
#include<dune/istl/io.hh>
//#include<dune/istl/paamg/amg.hh>

#include"../finiteelementmap/p0fem.hh"
#include"../finiteelementmap/p12dfem.hh"
#include"../finiteelementmap/pk2dfem.hh"
#include"../gridfunctionspace/gridfunctionspace.hh"
#include"../gridfunctionspace/gridfunctionspaceutilities.hh"
#include"../gridfunctionspace/interpolate.hh"
#include"../constraints/constraints.hh"
#include"../common/function.hh"
#include"../common/vtkexport.hh"
#include"../localoperator/laplacedirichletccfv.hh"
#include"../backend/backendselector.hh"
#include"../backend/istlvectorbackend.hh"
#include"../backend/istlmatrixbackend.hh"
#include"../backend/seqistlsolverbackend.hh"

#include"../gridoperator/gridoperator.hh"

#include <dune/pdelab/gridfunctionspace/vtk.hh>
#include <dune/pdelab/common/benchmarkhelper.hh>
#include <dune/pdelab/ordering/singlecodimleafordering.hh>
#include <dune/common/parametertreeparser.hh>

#include "../test/gridexamples.hh"


// define some grid functions to interpolate from
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

// define some boundary grid functions to define boundary conditions
template<typename GV>
class B
  : public Dune::PDELab::BoundaryGridFunctionBase<Dune::PDELab::BoundaryGridFunctionTraits<GV,int,1,
                                                                                           Dune::FieldVector<int,1> >,
                                                  B<GV> >
{
  const GV& gv;

public:
  typedef Dune::PDELab::BoundaryGridFunctionTraits<GV,int,1,Dune::FieldVector<int,1> > Traits;
  typedef Dune::PDELab::BoundaryGridFunctionBase<Traits,B<GV> > BaseT;

  B (const GV& gv_) : gv(gv_) {}

  template<typename I>
  inline void evaluate (const Dune::PDELab::IntersectionGeometry<I>& ig,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    y = 1; // all is Dirichlet boundary
  }

  //! get a reference to the GridView
  inline const GV& getGridView () const
  {
    return gv;
  }
};


template<typename GV, typename Mapper>
void test (const GV& gv, const Mapper& mapper, std::string filename, const bool solve, std::size_t runs)
{

  Dune::PDELab::BenchmarkHelper<> bh(filename,runs);

  for (std::size_t run = 0; run < runs; ++run)
    {

      bh.start_run(std::cout);

      typedef typename GV::Grid::ctype DF;
      typedef double RF;
      const int dim = GV::dimension;

      // instantiate finite element maps
      Dune::GeometryType gt;
      gt.makeCube(dim);
      typedef Dune::PDELab::P0LocalFiniteElementMap<DF,RF,dim> FEM;
      FEM fem(gt); // works only for cubes

      bh.start("global setup",std::cout);
      bh.start("GFS setup",std::cout);

      // make function space
      typedef Dune::PDELab::GridFunctionSpace<
        GV,
        FEM,
        Dune::PDELab::NoConstraints,
        Dune::PDELab::ISTLVectorBackend<>,
        Mapper
        > GFS;
      GFS gfs(gv,fem);
      gfs.name("u");
      // Output FV data as VTK cell data set
      gfs.setDataSetType(GFS::Output::cellData);

      bh.end("GFS setup",std::cout);
      bh.start("ordering update",std::cout);

      gfs.ordering();

      bh.end("ordering update",std::cout);
      bh.end("global setup",std::cout);

      typedef G<GV,RF> GType;
      GType g(gv);

      // make grid function operator
      typedef Dune::PDELab::LaplaceDirichletCCFV<GType> LO;
      LO lo(g);

      typedef Dune::PDELab::GridOperator<
        GFS,GFS,LO,
        Dune::PDELab::ISTLMatrixBackend,
        RF,RF,RF> GO;
      GO go(gfs,gfs,lo);

      // make coefficent Vector and initialize it from a function
      typedef typename GO::Traits::Domain V;
      V x0(gfs);
      x0 = 0.0;

      bh.start("interpolation",std::cout);

      Dune::PDELab::interpolate(g,gfs,x0);

      bh.end("interpolation",std::cout);
      bh.start("matrix creation",std::cout);

      // represent operator as a matrix
      typedef typename GO::Traits::Jacobian M;
      M m(go);
      m = 0.0;

      bh.end("matrix creation",std::cout);
      bh.start("jacobian",std::cout);

      go.jacobian(x0,m);

      bh.end("jacobian",std::cout);

      // evaluate residual w.r.t initial guess
      V r(gfs);
      r = 0.0;

      bh.start("residual",std::cout);

      go.residual(x0,r);

      bh.end("residual",std::cout);
      bh.start("solve",std::cout);

      V x(gfs,0.0);

      if (solve)
        {

          typedef typename M::Container ISTLM;
          typedef typename V::Container ISTLV;

          // make ISTL solver
          Dune::MatrixAdapter<ISTLM,ISTLV,ISTLV> opa(m.base());
          //  typedef Dune::PDELab::OnTheFlyOperator<V,V,GOS> ISTLOnTheFlyOperator;
          //  ISTLOnTheFlyOperator opb(gos);
          //  Dune::SeqSSOR<M,V,V> ssor(m,1,1.0);
          Dune::SeqILU0<ISTLM,ISTLV,ISTLV> ilu0(m.base(),1.0);
          //  Dune::Richardson<V,V> richardson(1.0);
          Dune::CGSolver<ISTLV> solvera(opa,ilu0,1E-10,5000,0);
          //  Dune::CGSolver<V> solverb(opb,richardson,1E-10,5000,2);
          Dune::InverseOperatorResult stat;

          // solve the jacobian system
          r *= -1.0; // need -residual
          solvera.apply(x.base(),r.base(),stat);
          x += x0;

        }

      bh.end("solve",std::cout);
      bh.start("I/O",std::cout);

      // output grid function with VTKWriter
      Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::nonconforming);
      Dune::PDELab::add_solution_to_vtk_writer(vtkwriter,gfs,x);
      vtkwriter.write(filename,Dune::VTKOptions::ascii);

      bh.end("I/O",std::cout);

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
std::string name(std::string grid_name, const Grid& grid, std::string backend_name, bool solve)
{
  std::stringstream n;
  n << "laplacedirichletccfv"
    << "_" << grid_name
    << "_" << Grid::dimension << "D"
    << "_l" << grid.maxLevel()
    << "_" << backend_name
    << (solve ? "_solve" : "_nosolve");
  return n.str();
}


int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper::instance(argc, argv);

    Dune::ParameterTree params;

    if (argc == 1)
      {
        std::cerr << "no parameter file passed, defaulting to laplacedirichletccfv.ini..." << std::endl;
        Dune::ParameterTreeParser::readINITree("laplacedirichletccfv.ini",params);
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

    // 2D - full ordering / backend infrastructure
    if (params.hasSub("Yasp_2D_full") && params.get<bool>("Yasp_2D_full.enabled"))
      {
        const Dune::ParameterTree& lp = params.sub("Yasp_2D_full");

        const int refine = lp.get<int>("refine");
        const int runs = lp.get("runs",global_runs);
        const bool solve = lp.get("solve",global_solve);

        // need a grid in order to test grid functions
        Dune::FieldVector<double,2> L(1.0);
        Dune::FieldVector<int,2> N(1);
        Dune::FieldVector<bool,2> B(false);
        Dune::YaspGrid<2> grid(L,N,B,0);
        grid.globalRefine(refine);

        typedef Dune::PDELab::DefaultLeafOrderingTag OrderingTag;
        OrderingTag ordering_tag;

        test(grid.leafView(),ordering_tag,name("Yasp",grid,"full",solve),solve,runs);
      }

    // 2D - simplified ordering / backend infrastructure
    if (params.hasSub("Yasp_2D_simple") && params.get<bool>("Yasp_2D_simple.enabled"))
      {
        const Dune::ParameterTree& lp = params.sub("Yasp_2D_simple");

        const int refine = lp.get<int>("refine");
        const int runs = lp.get("runs",global_runs);
        const bool solve = lp.get("solve",global_solve);

        // need a grid in order to test grid functions
        Dune::FieldVector<double,2> L(1.0);
        Dune::FieldVector<int,2> N(1);
        Dune::FieldVector<bool,2> B(false);
        Dune::YaspGrid<2> grid(L,N,B,0);
        grid.globalRefine(refine);

        typedef Dune::PDELab::SingleCodimMapper Mapper;
        Mapper mapper;

        test(grid.leafView(),mapper,name("Yasp",grid,"simple",solve),solve,runs);
      }

    // 3D - full ordering / backend infrastructure
    if (params.hasSub("Yasp_3D_full") && params.get<bool>("Yasp_3D_full.enabled"))
      {
        const Dune::ParameterTree& lp = params.sub("Yasp_3D_full");

        const int refine = lp.get<int>("refine");
        const int runs = lp.get("runs",global_runs);
        const bool solve = lp.get("solve",global_solve);

        // need a grid in order to test grid functions
        Dune::FieldVector<double,3> L(1.0);
        Dune::FieldVector<int,3> N(1);
        Dune::FieldVector<bool,3> B(false);
        Dune::YaspGrid<3> grid(L,N,B,0);
        grid.globalRefine(refine);

        typedef Dune::PDELab::DefaultLeafOrderingTag OrderingTag;
        OrderingTag ordering_tag;

        test(grid.leafView(),ordering_tag,name("Yasp",grid,"full",solve),solve,runs);
      }

    // 3D - simplified ordering / backend infrastructure
    if (params.hasSub("Yasp_3D_simple") && params.get<bool>("Yasp_3D_simple.enabled"))
      {
        const Dune::ParameterTree& lp = params.sub("Yasp_3D_simple");

        const int refine = lp.get<int>("refine");
        const int runs = lp.get("runs",global_runs);
        const bool solve = lp.get("solve",global_solve);

        // need a grid in order to test grid functions
        Dune::FieldVector<double,3> L(1.0);
        Dune::FieldVector<int,3> N(1);
        Dune::FieldVector<bool,3> B(false);
        Dune::YaspGrid<3> grid(L,N,B,0);
        grid.globalRefine(refine);

        typedef Dune::PDELab::SingleCodimMapper Mapper;
        Mapper mapper;

        test(grid.leafView(),mapper,name("Yasp",grid,"simple",solve),solve,runs);
      }

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
