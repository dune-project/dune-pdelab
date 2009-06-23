// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifdef HAVE_CONFIG_H
#include "config.h"     
#endif

#include <cmath>
#include <string>
#include <sstream>

#include <dune/common/float_cmp.hh>
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
#include "../finiteelementmap/edges03dfem.hh"
#include "../gridfunctionspace/gridfunctionspace.hh"
#include "../gridfunctionspace/gridfunctionspaceutilities.hh"
#include "../gridfunctionspace/interpolate.hh"

#include "gridexamples.hh"
#include "stddomains.hh"


template<typename GV, typename RF>
class U
  : public Dune::PDELab::AnalyticGridFunctionBase<
      Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,3>,
      U<GV,RF>
      >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,3> Traits;
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
double interpolationerror (const GV& gv, const FEM &fem, const std::string &name = "")
{
  typedef typename FEM::Traits::LocalFiniteElementType::Traits
    ::LocalBasisType::Traits::DomainFieldType D; // domain type
  typedef typename FEM::Traits::LocalFiniteElementType::Traits
    ::LocalBasisType::Traits::RangeFieldType R;  // range type

  typedef Dune::PDELab::GridFunctionSpace<GV, FEM> GFS;
  GFS gfs(gv,fem);                    // make grid function space

  typedef typename GFS::template VectorContainer<R>::Type X;
  X x(gfs,0.0);                       // make coefficient vector

  typedef U<GV,R> AFunc;
  AFunc u(gv);                      // make analytic function object
  Dune::PDELab::interpolateGlobal(u,gfs,x); // make x interpolate u

  typedef Dune::PDELab::DiscreteGridFunctionGlobal<GFS, X> IFunc;
  IFunc v(gfs,x);

  if(x.size() < 10)
    for(unsigned i = 0; i < x.size(); ++i)
      std::cout << "x[" << i << "] = " << x[i] << std::endl;

  if(name != "") {
    Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,1);  // plot result
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<AFunc>(u,"analytic"));
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<IFunc>(v,"interpolated"));
    vtkwriter.write(name,Dune::VTKOptions::ascii);
  }

  return l2difference(gv,u,v,4);
}

template<typename Grid>
void test(Dune::SmartPointer<Grid> grid, int &result, unsigned int maxelements, std::string name = "")
{
  if(name == "") name = grid->name();

  std::cout << std::endl
            << "Testing EdgeS03D interpolation with " << name << std::endl;

  name = "edges03dinterpolationglobal-" + name;

  typedef Dune::PDELab::EdgeS03DLocalFiniteElementMap<typename Grid::LeafGridView, double> FEM;

  std::cout << "interpolation level 0" << std::endl;
  double error0 = interpolationerror(grid->leafView(), FEM(grid->leafView()), name+"-coarse");
  double h0 = std::pow(1/double(grid->leafView().size(0)), 1/double(Grid::dimension));
  std::cout << "interpolation error: " 
            << std::setw(8) << grid->leafView().size(0) << " elements, h=" 
            << std::scientific << h0 << ", error="
            << std::scientific << error0 << std::endl;

  if(Dune::FloatCmp::eq(error0, 0.0)) {
    std::cerr << "Error: The analytic function was perfectly interpolated." << std::endl
              << "Error: This makes the interpolation convergence test meaningless." << std::endl
              << "Error: Please change this test program to use an analyting function which cannot be" << std::endl
              << "Error: represented exacly by this basis, i.e. something containing exp()" << std::endl;
    result = 1;
    return;
  }

  while((unsigned int)(grid->leafView().size(0)) < maxelements)
    grid->globalRefine(1);

  std::cout << "interpolation level " << grid->maxLevel() << std::endl;
  double errorf = interpolationerror(grid->leafView(), FEM(grid->leafView()), name+"-fine");
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

  if(total_convergence < 0.85) {
    std::cout << "Error: interpolation total convergence < 0.85" << std::endl;
    result = 1;
  }
}

int main(int argc, char** argv)
{
  using Dune::PDELab::UnitTetrahedronMaker;
  using Dune::PDELab::TriangulatedUnitCubeMaker;

  try{
    // default exitcode 77 (=skipped); returned in case none of the supported
    // Grids were found
    int result = 77;

#ifdef HAVE_ALBERTA
#if (ALBERTA_DIM != 3)
#error ALBERTA_DIM is not set to 3 -- please check the Makefile.am
#endif
    test(UnitTetrahedronMaker         <Dune::AlbertaGrid<3, 3>    >::create(),
         result, 50000, "alberta-tetrahedron");
    test(KuhnTriangulatedUnitCubeMaker<Dune::AlbertaGrid<3, 3>    >::create(),
          result, 50000, "alberta-cube");
#endif

#ifdef HAVE_ALUGRID
    test(UnitTetrahedronMaker         <Dune::ALUSimplexGrid<3, 3> >::create(),
         result, 10000, "alu-tetrahedron");
    test(KuhnTriangulatedUnitCubeMaker<Dune::ALUSimplexGrid<3, 3> >::create(),
         result, 10000, "alu-cube");
#endif // HAVE_ALUGRID

#ifdef HAVE_UG
    test(UnitTetrahedronMaker         <Dune::UGGrid<3>            >::create(),
         result, 10000, "ug-tetrahedron");
    test(KuhnTriangulatedUnitCubeMaker<Dune::UGGrid<3>            >::create(),
         result, 10000, "ug-cube");
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
