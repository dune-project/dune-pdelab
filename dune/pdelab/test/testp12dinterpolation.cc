// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifdef HAVE_CONFIG_H
#include "config.h"     
#endif

#include <cmath>
#include <string>
#include <sstream>

#include <dune/common/fvector.hh>
#include <dune/common/geometrytype.hh>
#include <dune/common/smartpointer.hh>

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

#include "stddomains.hh"


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

template<typename GV, typename U, typename V> 
double l2difference (const GV & gv, const U& u, const V &v, int qorder=1)
{
  // constants and types
  const int dim = GV::Grid::dimension;
  typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
  typedef typename GV::Grid::ctype ct;
  
  // loop over grid view
  double sum = 0.0;
  for (ElementIterator eit = gv.template begin<0>();
       eit!=gv.template end<0>(); ++eit)
  {

    Dune::GeometryType gt = eit->geometry().type();
    const Dune::QuadratureRule<ct,dim>& 
      rule = Dune::QuadratureRules<ct,dim>::rule(gt,qorder);

    for (typename Dune::QuadratureRule<ct,dim>::const_iterator qit=rule.begin(); 
         qit!=rule.end(); ++qit)
    {
      // evaluate the given grid functions at integration point
      typename U::Traits::RangeType u_val;
      u.evaluate(*eit,qit->position(),u_val);

      typename V::Traits::RangeType v_val;
      v.evaluate(*eit,qit->position(),v_val);

      // accumulate error
      v_val -= u_val;
      sum += v_val.two_norm2()*qit->weight()*
        eit->geometry().integrationElement(qit->position());
    }
  }
  return std::sqrt(sum);
}

template<typename GV, typename FEM>
void interpolationerror (const GV& gv, const FEM &fem)
{
  typedef typename FEM::Traits::LocalFiniteElementType::Traits
    ::LocalBasisType::Traits::DomainFieldType D; // domain type
  typedef typename FEM::Traits::LocalFiniteElementType::Traits
    ::LocalBasisType::Traits::RangeFieldType R;  // range type

  typedef Dune::PDELab::GridFunctionSpace<GV, FEM> GFS;    
  GFS gfs(gv,fem);                    // make grid function space

  typedef typename GFS::template VectorContainer<R>::Type X;
  X x(gfs,0.0);                       // make coefficient vector

  U<GV,R> u(gv);                      // make analytic function object
  Dune::PDELab::interpolate(u,gfs,x); // make x interpolate u

  Dune::PDELab::DiscreteGridFunction<GFS, X> v(gfs,x);

  std::cout.precision(8);
  std::cout << "interpolation error: " 
            << std::setw(8) << gv.size(0) << " elements " 
            << std::scientific << l2difference(gv,u,v,4) << std::endl;
}

template<typename Grid>
void test(Dune::SmartPointer<Grid> grid, int &result, std::string name = "")
{
  if(name == "") name = grid->name();
  name = "p12dinterpolation-"+name;

  typedef Dune::PDELab::P12DLocalFiniteElementMap<typename Grid::ctype, double> FEM;
  FEM fem;

  for (int l=0; l<5; l++)
  {
    interpolationerror(grid->leafView(), fem);
    grid->globalRefine(1);
  }

  result = 0;
}

int main(int argc, char** argv)
{
  using Dune::PDELab::UnitTriangleMaker;
  using Dune::PDELab::TriangulatedUnitSquareMaker;

  try{
    //Maybe initialize Mpi
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

    // default exitcode 77 (=skipped); returned in case none of the supported
    // Grids were found
    int result = 77;

#ifdef HAVE_ALBERTA
#if (ALBERTA_DIM != 2)
#error ALBERTA_DIM is not set to 2 -- please check the Makefile.am
#endif
    test(UnitTriangleMaker          <Dune::AlbertaGrid<2, 2>    >::create(),
         result, "alberta-triangle");
    test(TriangulatedUnitSquareMaker<Dune::AlbertaGrid<2, 2>    >::create(),
         result, "alberta-square");
#endif

#ifdef HAVE_ALUGRID
    test(UnitTriangleMaker          <Dune::ALUSimplexGrid<2, 2> >::create(),
         result, "alu-triangle");
    test(TriangulatedUnitSquareMaker<Dune::ALUSimplexGrid<2, 2> >::create(),
         result, "alu-square");
#endif // HAVE_ALUGRID

#ifdef HAVE_UG
    test(UnitTriangleMaker          <Dune::UGGrid<2>            >::create(),
         result, "ug-triangle");
    test(TriangulatedUnitSquareMaker<Dune::UGGrid<2>            >::create(),
         result, "ug-square");
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
