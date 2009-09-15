// -*- tab-width: 4; indent-tabs-mode: nil -*-

#ifdef HAVE_CONFIG_H
#include "config.h"     
#endif

#include <string>
#include <sstream>

#include <dune/common/fvector.hh>
#include <dune/common/mpihelper.hh>
#include <dune/common/smartpointer.hh>

#ifdef HAVE_ALBERTA
#include <dune/grid/albertagrid.hh>
#endif
#ifdef HAVE_ALUGRID
#include <dune/grid/alugrid.hh>
#endif
#ifdef HAVE_UG
#include <dune/grid/uggrid.hh>
#endif

#include "../finiteelementmap/edges03dfem.hh"

#include "gridexamples.hh"

const double eps = 1e-10;

const unsigned refine_limit = 100;

template<typename T, int n>
std::string fmt(const Dune::FieldVector<T, n>& v) {
  std::ostringstream s;
  if(n > 0)
    s << v[0];
  for(unsigned i = 1; i < n; ++i)
    s << ", " << v[i];
  return s.str();
}

template<typename U, typename V>
std::string fmt(const std::pair<U, V>& p) {
  std::ostringstream s;
  s << p.first << ", " << p.second;
  return s.str();
}

template<typename GV>
bool testFEM(const GV& gv, const std::string indent = "")
{
  typedef typename GV::Grid::ctype DF; // domain field type
  static const unsigned dimDomain = GV::dimension;
  typedef double RF;                   // range field type
  typedef Dune::PDELab::EdgeS03DLocalFiniteElementMap<GV,RF> FEM;
  typedef typename FEM::Traits::LocalFiniteElementType::Traits
    ::LocalBasisType LocalBasis;
  typedef typename LocalBasis::Traits::DomainType D;
  typedef typename LocalBasis::Traits::RangeType R;
  typedef typename GV::template Codim<0>::Iterator ElemIterator;
  typedef typename GV::template Codim<0>::Geometry Geometry;
  typedef Dune::GenericReferenceElements<DF, dimDomain> ReferenceElements;
  typedef Dune::GenericReferenceElement<DF, dimDomain> ReferenceElement;

  bool success = true;
  FEM fem(gv);   // maps entity to finite element

  // a map containing a for each edge a list of (element index, tangential
  // value) pairs.  An edge is identified by a the vertex indices of its end
  // points, where first vertex index < second vertex index.
  std::map<std::pair<unsigned, unsigned>, std::list<std::pair<unsigned, RF> > > values;

  std::vector<R> result;
  std::vector<D> edgeCenters;
  std::vector<D> edgeTangents;
  std::vector<std::pair<unsigned, unsigned> > edges;
  std::vector<std::pair<unsigned, unsigned> > edgesLocal;

  const ElemIterator end = gv.template end<0>();
  for(ElemIterator it = gv.template begin<0>(); it != end; ++it) {
    bool element_printed = false;
    const unsigned nEdges = it->template count<dimDomain-1>();
    const Geometry& geo = it->geometry();
    const ReferenceElement& ref = ReferenceElements::general(geo.type());
    const LocalBasis& localBasis = fem.find(*it).localBasis();
    const unsigned elemIndex = gv.indexSet().index(*it);
    
    edgeCenters.resize(nEdges);
    edgeTangents.resize(nEdges);
    edges.resize(nEdges);
    edgesLocal.resize(nEdges);
    for(unsigned i = 0; i < nEdges; ++i) {
      unsigned l0 = ref.subEntity(i,dimDomain-1, 0,dimDomain);
      unsigned l1 = ref.subEntity(i,dimDomain-1, 1,dimDomain);
      unsigned g0 = gv.indexSet().subIndex(*it, l0, dimDomain);
      unsigned g1 = gv.indexSet().subIndex(*it, l1, dimDomain);
      if(g0 > g1) {
        std::swap(g0, g1);
        std::swap(l0, l1);
      }
      D center = geo.corner(l0);  center += geo.corner(l1);
      center /= 2;
      D tangent = geo.corner(l1); tangent -= geo.corner(l0);
      tangent /= tangent.two_norm();

      edges[i] = std::make_pair(g0, g1);
      edgesLocal[i] = std::make_pair(l0, l1);
      edgeCenters[i] = center;
      edgeTangents[i] = tangent;
    }

    for(unsigned i = 0; i < nEdges; ++i) {
      bool edge_printed = false;
      localBasis.evaluateFunctionGlobal(edgeCenters[i], result, geo);
      for(unsigned j = 0; j < nEdges; ++j) {
        RF tComp = edgeTangents[i]*result[j];
        if(j == i)
          values[edges[i]].push_back(std::make_pair(elemIndex,tComp));
        else
          if(std::abs(tComp) > eps) {
            if(!element_printed) {
              std::cout << indent
                        << "Element " << elemIndex
                        << " @(" << fmt(geo.global(ref.position(0,0))) << ")"
                        << std::endl;
              element_printed = true;
            }
            if(!edge_printed) {
              std::cout << indent << "  "
                        << "Edge " << i
                        << " local vertices (" << fmt(edgesLocal[i]) << ")"
                        << " global vertices (" << fmt(edges[i]) << ")"
                        << " center (" << fmt(edgeCenters[i]) << ")"
                        << " tangent (" << fmt(edgeTangents[i]) << ")"
                        << std::endl;
              edge_printed = true;
            }
            std::cout << indent << "    "
                      << "Shape function " << j << ":"
                      << " Tangential component is non-zero"
                      << " (" << tComp << ")"
                      << std::endl;
            success = false;
          }
      }
    }
  }

  return success;
}

template<typename Grid>
void test(Dune::SmartPointer<Grid> grid, int &result, std::string name = "")
{
  if(name == "") name = grid->name();
  std::cout << "*** Checking " << name << " ***\n" << std::endl;

  for(;
      static_cast<unsigned>(grid->size(0)) < refine_limit;
      grid->globalRefine(1)) {
    std::cout << "  Level " << grid->maxLevel()
              << " with " << grid->size(0) << " elements" << std::endl;
    bool success = testFEM(grid->leafView(), "    ");
    if(!success)
      result = 1;
    else if(result == 77)
      result = 0;
    std::cout << std::endl;
  }
}

int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

    // default exitcode 77 (=skipped); returned in case none of the supported
    // Grids were found
    int result = 77;

#ifdef HAVE_ALBERTA
#if (ALBERTA_DIM != 3)
#error ALBERTA_DIM is not set to 3 -- please check the Makefile.am
#endif
    test(UnitTetrahedronMaker         <Dune::AlbertaGrid<3, 3>    >::create(),
         result, "alberta-tetrahedron");
    test(KuhnTriangulatedUnitCubeMaker<Dune::AlbertaGrid<3, 3>    >::create(),
         result, "alberta-cube");
#endif

#ifdef HAVE_ALUGRID
    test(UnitTetrahedronMaker         <Dune::ALUSimplexGrid<3, 3> >::create(),
         result, "alu-tetrahedron");
    test(TriangulatedUnitCubeMaker    <Dune::ALUSimplexGrid<3, 3> >::create(),
         result, "alu-cube");
#endif // HAVE_ALUGRID

#ifdef HAVE_UG
    test(UnitTetrahedronMaker         <Dune::UGGrid<3>            >::create(),
         result, "ug-tetrahedron");
    test(TriangulatedUnitCubeMaker    <Dune::UGGrid<3>            >::create(),
         result, "ug-cube");
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
