// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include <memory>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/powergridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/compositegridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/vectorgridfunctionspace.hh>
#include <dune/pdelab/ordering/interleavedordering.hh>
#include <dune/pdelab/ordering/chunkedblockordering.hh>
#include <dune/pdelab/ordering/permutedordering.hh>

// test chunked block ordering with a taylor hood space

int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper::instance(argc, argv);

    std::size_t chunk_size = argc > 1 ? atoi(argv[1]) : 5;

    // 2D
    {
      std::cout << "2D tests" << std::endl;
      // need a grid in order to test grid functions
      // typedef Dune::YaspGrid<2> Grid;
      Dune::FieldVector<double,2> L(1.0);
      std::array<int,2> N{{2,2}};
      typedef Dune::YaspGrid<2> Grid;
      Grid grid(L,N);

      typedef Grid::LeafGridView GV;
      GV gv = grid.leafGridView();

      typedef Dune::PDELab::QkLocalFiniteElementMap<GV,double,double,2> V_FEM;
      typedef Dune::PDELab::QkLocalFiniteElementMap<GV,double,double,1> P_FEM;
      V_FEM v_fem(gv);
      P_FEM p_fem(gv);

      typedef Dune::PDELab::ISTL::VectorBackend<> V_Component_Backend;
      typedef Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed,2> V_Backend;
      typedef Dune::PDELab::ISTL::VectorBackend<> P_Backend;
      typedef Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::bcrs> TH_Backend;

      typedef Dune::PDELab::DefaultLeafOrderingTag V_Component_Ordering;
      typedef Dune::PDELab::ordering::Chunked<Dune::PDELab::EntityBlockedOrderingTag> V_Ordering;
      typedef Dune::PDELab::DefaultLeafOrderingTag P_Ordering;
      typedef Dune::PDELab::LexicographicOrderingTag TH_Ordering;

      typedef Dune::PDELab::VectorGridFunctionSpace<
        GV,
        V_FEM,
        2,
        V_Backend,
        V_Component_Backend,
        Dune::PDELab::NoConstraints,
        V_Ordering,
        V_Component_Ordering
        > V_GFS;
      V_GFS v_gfs(gv,v_fem,V_Backend(),V_Component_Backend(),V_Ordering(chunk_size),V_Component_Ordering());

      typedef Dune::PDELab::GridFunctionSpace<
        GV,
        P_FEM,
        Dune::PDELab::NoConstraints,
        P_Backend,
        P_Ordering
        > P_GFS;
      P_GFS p_gfs(gv,p_fem,P_Backend(),P_Ordering());

      typedef Dune::PDELab::CompositeGridFunctionSpace<
        TH_Backend,
        TH_Ordering,
        V_GFS,
        P_GFS
        > TH_GFS;
      TH_GFS th_gfs(TH_Backend(),TH_Ordering(),v_gfs,p_gfs);

      th_gfs.update();

      std::cout << th_gfs.ordering().size() << " " << th_gfs.ordering().blockCount() << std::endl;
      std::cout << th_gfs.ordering().template child<0>().size() << " " << th_gfs.ordering().template child<0>().blockCount() << std::endl;

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
