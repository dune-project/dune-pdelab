// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifdef HAVE_CONFIG_H
#include "config.h"     
#endif

#include <cmath>
#include <string>
#include <sstream>

#include <dune/common/fvector.hh>
#include <dune/common/geometrytype.hh>
#include <dune/common/mpihelper.hh>
#include <dune/common/shared_ptr.hh>

#include <dune/grid/common/quadraturerules.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#ifdef HAVE_ALBERTA
#include <dune/grid/albertagrid.hh>
#endif
#ifdef HAVE_ALUGRID
#include <dune/grid/alugrid.hh>
#endif
#ifdef HAVE_UG
#include <dune/grid/uggrid.hh>
#endif

#include "../common/function.hh"
#include "../finiteelementmap/p12dfem.hh"
#include "../gridfunctionspace/gridfunctionspace.hh"
#include "../gridfunctionspace/gridfunctionspaceutilities.hh"
#include "../gridfunctionspace/interpolate.hh"

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
    ::LocalBasisType::Traits::DomainFieldType D; // domain type
  typedef typename FEM::Traits::FiniteElementType::Traits
    ::LocalBasisType::Traits::RangeFieldType R;  // range type

  typedef Dune::PDELab::GridFunctionSpace<GV, FEM> GFS;    
  GFS gfs(gv,fem);                    // make grid function space

  typedef typename GFS::template VectorContainer<R>::Type X;
  X x(gfs,0.0);                       // make coefficient vector

  U<GV,R> u(gv);                      // make analytic function object
  Dune::PDELab::interpolate(u,gfs,x); // make x interpolate u

  Dune::PDELab::DiscreteGridFunction<GFS, X> v(gfs,x);

  return l2difference(u,v,4);
}

template<typename Grid>
void test(Dune::shared_ptr<Grid> grid, int &result, unsigned int maxelements, std::string name)
{
  std::cout << std::endl
            << "Testing P12D interpolation with " << name << std::endl;

  typedef Dune::PDELab::P12DLocalFiniteElementMap<typename Grid::ctype, double> FEM;
  FEM fem;

  std::cout << "interpolation level 0" << std::endl;
  double error0 = interpolationerror(grid->leafView(), fem);
  double h0 = std::pow(1/double(grid->leafView().size(0)), 1/double(Grid::dimension));
  std::cout << "interpolation error: " 
            << std::setw(8) << grid->leafView().size(0) << " elements, h=" 
            << std::scientific << h0 << ", error="
            << std::scientific << error0 << std::endl;

  while((unsigned int)(grid->leafView().size(0)) < maxelements)
    grid->globalRefine(1);

  std::cout << "interpolation level " << grid->maxLevel() << std::endl;
  double errorf = interpolationerror(grid->leafView(), fem);
  double hf = std::pow(1/double(grid->leafView().size(0)), 1/double(Grid::dimension));
  std::cout << "interpolation error: " 
            << std::setw(8) << grid->leafView().size(0) << " elements, h=" 
            << std::scientific << hf << ", error="
            << std::scientific << errorf << std::endl;

  double total_convergence = std::log(errorf/error0)/std::log(hf/h0);
  std::cout << "interpolation total convergence: "
            << std::scientific << total_convergence << std::endl;

  if(result != 1)
    result = 0;

  if(total_convergence < 1.7) {
    std::cout << "Error: interpolation total convergence < 1.7" << std::endl;
    result = 1;
  }
}

int main(int argc, char** argv)
{
  try{
    Dune::MPIHelper::instance(argc, argv);

    // default exitcode 77 (=skipped); returned in case none of the supported
    // Grids were found
    int result = 77;

#ifdef HAVE_ALBERTA
#if (ALBERTA_DIM != 2)
#error ALBERTA_DIM is not set to 2 -- please check the Makefile.am
#endif
    test(UnitTriangleMaker          <Dune::AlbertaGrid<2, 2>    >::create(),
         result, 250000, "alberta-triangle");
    test(TriangulatedUnitSquareMaker<Dune::AlbertaGrid<2, 2>    >::create(),
         result, 250000, "alberta-square");
#endif

#ifdef HAVE_ALUGRID
    test(UnitTriangleMaker          <Dune::ALUSimplexGrid<2, 2> >::create(),
         result, 250000, "alu-triangle");
    test(TriangulatedUnitSquareMaker<Dune::ALUSimplexGrid<2, 2> >::create(),
         result, 250000, "alu-square");
#endif // HAVE_ALUGRID

#ifdef HAVE_UG
    test(UnitTriangleMaker          <Dune::UGGrid<2>            >::create(),
         result, 250000, "ug-triangle");
    test(TriangulatedUnitSquareMaker<Dune::UGGrid<2>            >::create(),
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
