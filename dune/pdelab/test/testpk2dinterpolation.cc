// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cmath>
#include <string>
#include <sstream>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#if HAVE_ALBERTA
#include <dune/grid/albertagrid.hh>
#endif
#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif
#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif

#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/finiteelementmap/pkfem.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>

#include "gridexamples.hh"
#include "l2difference.hh"

template<typename GV, typename RF>
class U
  : public Dune::PDELab::AnalyticGridFunctionBase<
      Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
      U<GV,RF>
      >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,U<GV,RF> > Base;

  U (const GV& gv)
    : Base(gv)
  {}

  inline void evaluateGlobal (const typename Traits::DomainType& x,
                              typename Traits::RangeType& y) const
  {
    typename Traits::DomainType center(0.5);
    center -= x;
    y = exp(-3.0*center.two_norm2());
  }
};

template<typename GV, typename FEM>
double interpolationerror (const GV& gv, const FEM &fem)
{
  typedef typename FEM::Traits::FiniteElementType::Traits
    ::LocalBasisType::Traits::RangeFieldType R;  // range type

  typedef Dune::PDELab::GridFunctionSpace<GV, FEM> GFS;
  GFS gfs(gv,fem);                    // make grid function space

  using X = Dune::PDELab::Backend::Vector<GFS, R>;
  X x(gfs,0.0);                       // make coefficient vector

  U<GV,R> u(gv);                      // make analytic function object
  Dune::PDELab::interpolate(u,gfs,x); // make x interpolate u

  Dune::PDELab::DiscreteGridFunction<GFS, X> v(gfs,x);

  return l2difference(u,v,4);
}

template<int k, typename Grid>
void run_test(std::shared_ptr<Grid> grid, int &result, unsigned int maxelements, std::string name)
{
  std::cout << std::endl
            << "Testing P" << k << "2D interpolation with " << name << std::endl;

  typedef Dune::PDELab::PkLocalFiniteElementMap<typename Grid::LeafGridView, double, double, k> FEM;
  FEM fem(grid->leafGridView());

  std::cout << "interpolation level 0" << std::endl;
  double error0 = interpolationerror(grid->leafGridView(), fem);
  double h0 = std::pow(1/double(grid->leafGridView().size(0)), 1/double(Grid::dimension));
  std::cout << "interpolation error: "
            << std::setw(8) << grid->leafGridView().size(0) << " elements, h="
            << std::scientific << h0 << ", error="
            << std::scientific << error0 << std::endl;

  while((unsigned int)(grid->leafGridView().size(0)) < maxelements)
    grid->globalRefine(1);

  std::cout << "interpolation level " << grid->maxLevel() << std::endl;
  double errorf = interpolationerror(grid->leafGridView(), fem);
  double hf = std::pow(1/double(grid->leafGridView().size(0)), 1/double(Grid::dimension));
  std::cout << "interpolation error: "
            << std::setw(8) << grid->leafGridView().size(0) << " elements, h="
            << std::scientific << hf << ", error="
            << std::scientific << errorf << std::endl;

  double total_convergence = std::log(errorf/error0)/std::log(hf/h0);
  std::cout << "interpolation total convergence: "
            << std::scientific << total_convergence << std::endl;

  if(result != 1)
    result = 0;

  const double min_convergence[] = {1.7, 2.7, 3.7};

  if(total_convergence < min_convergence[k-1]) {
    std::cout << "Error: interpolation total convergence < " << min_convergence[k-1] << std::endl;
    result = 1;
  }
}

template<typename GridFactory>
void test(const GridFactory& gf, int &result, unsigned int maxelements, std::string name)
{
  run_test<1>(gf.create(),result,maxelements,name);
  run_test<2>(gf.create(),result,maxelements,name);
  run_test<3>(gf.create(),result,maxelements,name);
}

int main(int argc, char** argv)
{
  try{
    Dune::MPIHelper::instance(argc, argv);

    // default exitcode 77 (=skipped); returned in case none of the supported
    // Grids were found
    int result = 77;

#if HAVE_ALBERTA
#if (ALBERTA_DIM != 2)
#error ALBERTA_DIM is not set to 2 -- please check the Makefile.am
#endif
    test(UnitTriangleMaker          <Dune::AlbertaGrid<2, 2>    >(),
         result, 250000, "alberta-triangle");
    test(TriangulatedUnitSquareMaker<Dune::AlbertaGrid<2, 2>    >(),
         result, 250000, "alberta-square");
#endif // HAVE_ALBERTA

#if HAVE_DUNE_ALUGRID
    using ALUType = Dune::ALUGrid<2, 2, Dune::simplex, Dune::nonconforming>;
    auto alugrid = Dune::StructuredGridFactory<ALUType>::createSimplexGrid(Dune::FieldVector<ALUType::ctype, 2>(0.0), Dune::FieldVector<ALUType::ctype, 2>(1.0), Dune::make_array(1u, 1u));
    run_test<1>(alugrid, result, 25000, "alu-triangle");
    run_test<2>(alugrid, result, 25000, "alu-triangle");
    run_test<3>(alugrid, result, 25000, "alu-triangle");
#endif

#if HAVE_UG
    test(UnitTriangleMaker          <Dune::UGGrid<2>            >(),
         result, 250000, "ug-triangle");
    test(TriangulatedUnitSquareMaker<Dune::UGGrid<2>            >(),
         result, 250000, "ug-square");
#endif // HAVE_ALBERTA

    return result;
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
