// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/fvector.hh>
#include <dune/common/typetraits.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/test/gridexamples.hh>


// Include snippets for the different FEMs

#include "rtbdmfem.hh"
#include "opbfem.hh"
#include "pkfem.hh"
#include "rannacherturekfem.hh"


// Run unit tests for given FEM
template<typename FEM, typename GV, typename Constraints, typename VBE>
void test_fem(const FEM& fem, GV gv, const Constraints& constraints, const VBE& vbe)
{

  typedef Dune::PDELab::GridFunctionSpace<
    GV,
    FEM,
    Constraints,
    VBE
    > GFS;

  GFS gfs(gv,fem,constraints,vbe);

  std::cout << gfs.ordering().size() << std::endl;

}


template<typename RF, typename Constraints, typename VBE>
void test_2d_cube(const Constraints& constraints, const VBE& vbe)
{

  // make grid
  Dune::FieldVector<double,2> L(1.0);
  std::array<int,2> N;
  std::fill(N.begin(),N.end(),1);
  Dune::YaspGrid<2> grid(L,N);
  grid.globalRefine(3);

  // get view
  typedef Dune::YaspGrid<2>::LeafGridView GV;
  auto gv=grid.leafGridView();

  typedef GV::Grid::ctype DF;

  typedef typename FEM_FACTORY::template FEM<GV,DF,RF,Dune::GeometryType::cube>::pointer PFEM;

  PFEM pfem = FEM_FACTORY::template create<GV,DF,RF,Dune::GeometryType::cube>(gv);

  test_fem(*pfem,gv,constraints,vbe);

}


template<typename RF, typename Constraints, typename VBE>
void test_3d_cube(const Constraints& constraints, const VBE& vbe)
{

  // make grid
  Dune::FieldVector<double,3> L(1.0);
  std::array<int,3> N;
  std::fill(N.begin(),N.end(),1);

  Dune::YaspGrid<3> grid(L,N);
  grid.globalRefine(3);

  // get view
  typedef Dune::YaspGrid<3>::LeafGridView GV;
  const GV& gv=grid.leafGridView();

  typedef GV::Grid::ctype DF;

  typedef typename FEM_FACTORY::template FEM<GV,DF,RF,Dune::GeometryType::cube>::pointer PFEM;

  PFEM pfem = FEM_FACTORY::template create<GV,DF,RF,Dune::GeometryType::cube>(gv);

  test_fem(*pfem,gv,constraints,vbe);

}


template<typename RF, typename Constraints, typename VBE>
void test_2d_simplex(const Constraints& constraints, const VBE& vbe)
{

#if HAVE_DUNE_ALUGRID

    {
      // make grid
      using ALUType = Dune::ALUGrid<2, 2, Dune::simplex, Dune::nonconforming>;
      auto alugrid = Dune::StructuredGridFactory<ALUType>::createSimplexGrid(Dune::FieldVector<ALUType::ctype, 2>(0.0), Dune::FieldVector<ALUType::ctype, 2>(1.0), Dune::make_array(1u, 1u));
      alugrid->globalRefine(3);

      // get view
      typedef ALUType::LeafGridView GV;
      auto gv = alugrid->leafGridView();

      typedef GV::Grid::ctype DF;

      typedef typename FEM_FACTORY::template FEM<GV,DF,RF,Dune::GeometryType::simplex>::pointer PFEM;

      PFEM pfem = FEM_FACTORY::template create<GV,DF,RF,Dune::GeometryType::simplex>(gv);

      test_fem(*pfem,gv,constraints,vbe);
    }

#elif HAVE_UG

    {
      // make grid
      typedef Dune::UGGrid<2> Grid;
      std::shared_ptr<Grid> gridptr = TriangulatedUnitSquareMaker<Grid>::create();
      Grid& grid = *gridptr;
      grid.globalRefine(3);

      // get view
      typedef Grid::LeafGridView GV;
      auto gv=grid.leafGridView();

      typedef GV::Grid::ctype DF;

      typedef typename FEM_FACTORY::template FEM<GV,DF,RF,Dune::GeometryType::simplex>::pointer PFEM;

      PFEM pfem = FEM_FACTORY::template create<GV,DF,RF,Dune::GeometryType::simplex>(gv);

      test_fem(*pfem,gv,constraints,vbe);
    }

#else

#warning Could not find supported 2D simplex grid, 2D simplex tests will be skipped.
    std::cerr << "Warning: 2D simplex tests were skipped." << std::endl;

#endif

}


template<typename RF, typename Constraints, typename VBE>
void test_3d_simplex(const Constraints& constraints, const VBE& vbe)
{

#if HAVE_DUNE_ALUGRID

    {
      using ALUType = Dune::ALUGrid<3, 3, Dune::simplex, Dune::nonconforming>;
      auto alugrid = Dune::StructuredGridFactory<ALUType>::createSimplexGrid(Dune::FieldVector<ALUType::ctype, 3>(0.0), Dune::FieldVector<ALUType::ctype, 3>(1.0), Dune::make_array(1u, 1u, 1u));
      alugrid->globalRefine(3);

      // get view
      typedef ALUType::LeafGridView GV;
      auto gv = alugrid->leafGridView();

      typedef GV::Grid::ctype DF;

      typedef typename FEM_FACTORY::template FEM<GV,DF,RF,Dune::GeometryType::simplex>::pointer PFEM;

      PFEM pfem = FEM_FACTORY::template create<GV,DF,RF,Dune::GeometryType::simplex>(gv);

      test_fem(*pfem,gv,constraints,vbe);
    }

#elif HAVE_UG

    {
      // make grid
      typedef Dune::UGGrid<3> Grid;
      std::shared_ptr<Grid> gridptr = TriangulatedUnitCubeMaker<Grid>::create();
      Grid& grid = *gridptr;
      grid.globalRefine(3);

      // get view
      typedef Grid::LeafGridView GV;
      auto gv=grid.leafGridView();

      typedef GV::Grid::ctype DF;

      typedef typename FEM_FACTORY::template FEM<GV,DF,RF,Dune::GeometryType::simplex>::pointer PFEM;

      PFEM pfem = FEM_FACTORY::template create<GV,DF,RF,Dune::GeometryType::simplex>(gv);

      test_fem(*pfem,gv,constraints,vbe);
    }

#else

#warning Could not find supported 3D simplex grid, 3D simplex tests will be skipped.
    std::cerr << "Warning: 3D simplex tests were skipped." << std::endl;

#endif

}


int main(int argc, char** argv)
{
  try{

    Dune::MPIHelper::instance(argc,argv);

    typedef Dune::PDELab::NoConstraints Constraints;
    Constraints constraints DUNE_UNUSED;

    typedef Dune::PDELab::ISTL::VectorBackend<> VBE;
    VBE vbe DUNE_UNUSED;

    typedef double RF DUNE_UNUSED;

#if FEM_DIM == 1

#ifdef FEM_CUBE
    test_1d_cube<RF>(constraints,vbe);
#endif

#ifdef FEM_SIMPLEX
    test_1d_simplex<RF>(constraints,vbe);
#endif

#elif FEM_DIM == 2

#ifdef FEM_CUBE
    test_2d_cube<RF>(constraints,vbe);
#endif

#ifdef FEM_SIMPLEX
    test_2d_simplex<RF>(constraints,vbe);
#endif

#elif FEM_DIM == 3

#ifdef FEM_CUBE
    test_3d_cube<RF>(constraints,vbe);
#endif

#ifdef FEM_SIMPLEX
    test_3d_simplex<RF>(constraints,vbe);
#endif

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
