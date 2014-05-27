// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/pdelab/finiteelementmap/brezzidouglasmarinifem.hh>
#include <dune/pdelab/finiteelementmap/bdm1simplex2dfem.hh>
#include <dune/pdelab/finiteelementmap/bdm1cube2dfem.hh>

#include "gridexamples.hh"

// Run unit tests for given FEM
template<typename GV, typename FEM, typename BaseFEM>
void test(GV gv, const FEM& fem, const BaseFEM& baseFEM)
{
  typename GV::template Codim<0>::Iterator it = gv.template begin<0>();
  const typename FEM::Traits::FiniteElement& fe = fem.find(*it);
  const typename BaseFEM::Traits::FiniteElement& bfe = baseFEM.find(*it);
  dune_static_assert((Dune::is_same<typename FEM::Traits::FiniteElement,typename BaseFEM::Traits::FiniteElement>::value),
                     "Implementation error in BrezziDouglasMariniLocalFiniteElementMap: picked wrong base FEM");
  if (fe.localBasis().size() != bfe.localBasis().size())
    DUNE_THROW(Dune::InvalidStateException,"finite elements should be of same size");
}


int main(int argc, char** argv)
{
  try{

    Dune::MPIHelper::instance(argc,argv);

    // 2D cube tests
    {
      // make grid
      Dune::FieldVector<double,2> L(1.0);
      Dune::array<int,2> N(Dune::fill_array<int,2>(1));
      Dune::YaspGrid<2> grid(L,N);
      grid.globalRefine(3);

      // get view
      typedef Dune::YaspGrid<2>::LeafGridView GV;
      const GV& gv=grid.leafGridView();

      typedef GV::Grid::ctype DF;
      typedef double RF;

      typedef Dune::PDELab::BrezziDouglasMariniLocalFiniteElementMap<GV,DF,RF,1> BDM1FEM;
      typedef Dune::PDELab::BDM1Cube2DLocalFiniteElementMap<GV,DF,RF> BDM1BASEFEM;
      BDM1FEM bdm1_fem(gv);
      BDM1BASEFEM bdm1_base_fem(gv);
      test(gv,bdm1_fem,bdm1_base_fem);

    }

#if HAVE_ALUGRID

    {
      // 2D simplex tests
      ALUUnitSquare grid;
      grid.globalRefine(4);

      // get view
      typedef ALUUnitSquare::LeafGridView GV;
      const GV& gv=grid.leafGridView();

      typedef GV::Grid::ctype DF;
      typedef double RF;

      typedef Dune::PDELab::BrezziDouglasMariniLocalFiniteElementMap<GV,DF,RF,1> BDM1FEM;
      typedef Dune::PDELab::BDM1Simplex2DLocalFiniteElementMap<GV,DF,RF> BDM1BASEFEM;

      BDM1FEM bdm1_fem(gv);
      BDM1BASEFEM bdm1_base_fem(gv);
      test(gv,bdm1_fem,bdm1_base_fem);

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
