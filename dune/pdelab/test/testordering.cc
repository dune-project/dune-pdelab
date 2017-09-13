// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include <memory>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/alugrid/grid.hh>
#include <dune/grid/common/scsgmapper.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/pdelab/finiteelementmap/p0fem.hh>
#include <dune/pdelab/finiteelementmap/pkfem.hh>
#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/pdelab/finiteelementmap/variablemonomfem.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/powergridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/compositegridfunctionspace.hh>
#include <dune/pdelab/ordering/interleavedordering.hh>

#include <dune/pdelab/backend/istl.hh>

#include <dune/pdelab/constraints/common/constraints.hh>

template<typename GFS>
void check_ordering(const GFS& gfs)
{
    const typename GFS::Ordering& ordering = gfs.ordering();

    Dune::PDELab::LocalFunctionSpace<GFS> lfs(gfs);

    typedef typename GFS::Traits::GridView GV;

    for (typename GV::template Codim<0>::Iterator it = gfs.gridView().template begin<0>();
         it != gfs.gridView().template end<0>(); ++it)
    {
      lfs.bind(*it);

      std::vector<typename GFS::Ordering::Traits::DOFIndex> vdi(lfs.size());
      std::vector<typename GFS::Ordering::Traits::ContainerIndex> vci(lfs.size());
      for (unsigned i = 0; i < lfs.size(); ++i)
      {
        vdi[i] = lfs.dofIndex(i);
      }
      ordering.map_lfs_indices(vdi.begin(),vdi.end(),vci.begin());

      for (unsigned i = 0; i < lfs.size(); ++i)
      {
        const typename GFS::Ordering::Traits::DOFIndex& di = lfs.dofIndex(i);
        typename GFS::Ordering::Traits::ContainerIndex ci;
        ordering.mapIndex(di.view(),ci);
        std::cout << di << "    " << vci[i] << "    " << ci << std::endl;
      }
    }

    using V = Dune::PDELab::Backend::Vector<GFS,double>;
    V x(gfs);
    x = 0.0;
    std::cout << std::endl;
}


// test function trees
template<int dim, bool cube>
struct test;

template<>
struct test<2,true> {
  template<class GV>
  static void testleafgridfunction(const GV& gv)
  {
    // instantiate finite element maps
    auto gt = Dune::GeometryTypes::quadrilateral;
    typedef Dune::PDELab::P0LocalFiniteElementMap<float,double,GV::dimension> P0FEM;
    P0FEM p0fem(gt);
    typedef Dune::PDELab::QkLocalFiniteElementMap<GV,float,double,1> Q12DFEM;
    Q12DFEM q12dfem(gv);
    typedef Dune::PDELab::QkLocalFiniteElementMap<GV,float,double,2> Q22DFEM;
    Q22DFEM q22dfem(gv);
    typedef Dune::SingleCodimSingleGeomTypeMapper<GV, 0> CellMapper;
    CellMapper cellmapper(gv);
    typedef Dune::PDELab::VariableMonomLocalFiniteElementMap<CellMapper,float,double,GV::dimension> MonomFEM;
    MonomFEM monomfem(cellmapper,2);

    typedef Dune::PDELab::NoConstraints CON;

    typedef Dune::PDELab::ISTL::VectorBackend<> VBE;

    // make a grid function space
    typedef Dune::PDELab::GridFunctionSpace<GV,P0FEM,CON,VBE> P0GFS;
    P0GFS p0gfs(gv,p0fem);
    typedef Dune::PDELab::GridFunctionSpace<GV,Q12DFEM,CON,VBE> GFS1;
    GFS1 gfs1(gv,q12dfem);
    typedef Dune::PDELab::GridFunctionSpace<GV,Q22DFEM,CON,VBE> GFS2;
    GFS2 gfs2(gv,q22dfem);
    typedef Dune::PDELab::GridFunctionSpace<GV,MonomFEM,CON,VBE> GFS3;
    GFS3 gfs3(gv,monomfem);

    typedef Dune::PDELab::PowerGridFunctionSpace<GFS1,3,VBE,Dune::PDELab::InterleavedOrderingTag> P1GFS;

    {
      // keep this block around to test the vector-based constructor of the ordering tag
      std::vector<std::size_t> p1_gfs_block_sizes(3);
      std::fill(p1_gfs_block_sizes.begin(),p1_gfs_block_sizes.end(),1);
      P1GFS p1gfs(gfs1,gfs1,gfs1,VBE(),p1_gfs_block_sizes);
    }

    P1GFS p1gfs(gfs1,gfs1,gfs1,VBE(),{{1,1,1}});

    typedef Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed,6> NVBE;

    typedef Dune::PDELab::PowerGridFunctionSpace<P1GFS,2,NVBE,Dune::PDELab::InterleavedOrderingTag> PGFS;
    std::vector<std::size_t> p_gfs_block_sizes(2);
    std::fill(p_gfs_block_sizes.begin(),p_gfs_block_sizes.end(),3);
    PGFS pgfs(p1gfs,p1gfs,NVBE(),p_gfs_block_sizes);

    typedef Dune::PDELab::PowerGridFunctionSpace<GFS1,3,VBE,Dune::PDELab::EntityBlockedOrderingTag> P2GFS;
    P2GFS p2gfs(gfs1,gfs1,gfs1);

    typedef Dune::PDELab::PowerGridFunctionSpace<P2GFS,2,NVBE,Dune::PDELab::EntityBlockedOrderingTag> P3GFS;
    P3GFS p3gfs(p2gfs,p2gfs,NVBE());


    // check_ordering(gfs1);
    check_ordering(gfs3);
    return;
    check_ordering(pgfs);
    check_ordering(p3gfs);

    typedef Dune::PDELab::PowerGridFunctionSpace<GFS1,3,VBE,Dune::PDELab::EntityBlockedOrderingTag> EBPGFS1;
    EBPGFS1 ebpgfs1(gfs1);

    check_ordering(ebpgfs1);

    typedef Dune::PDELab::CompositeGridFunctionSpace<
      VBE,
      Dune::PDELab::EntityBlockedOrderingTag,
      GFS1,
      EBPGFS1,
      GFS1
      > EBCGFS1;

    EBCGFS1 ebcgfs1(gfs1,ebpgfs1,gfs1);

    check_ordering(ebcgfs1);
  }
};

template<>
struct test<2,false> {
  template<class GV>
  static void testleafgridfunction(const GV& gv)
  {
    // instantiate finite element maps
    typedef Dune::SingleCodimSingleGeomTypeMapper<GV, 0> CellMapper;
    CellMapper cellmapper(gv);
    typedef Dune::PDELab::VariableMonomLocalFiniteElementMap<CellMapper,float,double,GV::dimension> MonomFEM;
    MonomFEM monomfem(cellmapper,gv.indexSet().types(0)[0],2);

    typedef Dune::PDELab::NoConstraints CON;

    typedef Dune::PDELab::ISTL::VectorBackend<> VBE;

    // make a grid function space
    typedef Dune::PDELab::GridFunctionSpace<GV,MonomFEM,CON,VBE> GFS3;
    GFS3 gfs3(gv,monomfem);
    check_ordering(gfs3);
  }
};

template<>
struct test<3,true> {
  template<class GV>
  static void testleafgridfunction(const GV& gv)
  {
    // instantiate finite element maps
    auto gt = Dune::GeometryTypes::hexahedron;
    typedef Dune::PDELab::P0LocalFiniteElementMap<float,double,GV::dimension> P0FEM;
    P0FEM p0fem(gt);
    typedef Dune::PDELab::QkLocalFiniteElementMap<GV,float,double,1> Q1FEM;
    Q1FEM q1fem(gv);

    // make a grid function space
    typedef Dune::PDELab::GridFunctionSpace<GV,P0FEM> P0GFS;
    P0GFS p0gfs(gv,p0fem);
// Doesn't work, we need a grid with triangular elemets for that
//  typedef Dune::PDELab::GridFunctionSpace<GV,P1FEM> P1GFS;
//  P1GFS p1gfs(gv,p1fem);
    typedef Dune::PDELab::GridFunctionSpace<GV,Q1FEM> Q1GFS;
    Q1GFS q1gfs(gv,q1fem);
  }
};

template<bool cube, class GV>
void testleafgridfunction(const GV& gv)
{
  test<GV::dimension,cube>::testleafgridfunction(gv);
}


int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper::instance(argc, argv);

    // 2D
    {
      std::cout << "2D tests" << std::endl;
      // need a grid in order to test grid functions
      // typedef Dune::YaspGrid<2> Grid;
      typedef Dune::ALUGrid<2,2,Dune::cube,Dune::nonconforming> Grid;
      Dune::FieldVector<double,2> l(0.0);
      Dune::FieldVector<double,2> u(1.0);
      std::array<unsigned int,2> N = {{1,1}};
      std::shared_ptr<Grid> grid = Dune::StructuredGridFactory<Grid>::createCubeGrid(l,u,N);
      grid->globalRefine(1);

      std::cout << Dune::GlobalGeometryTypeIndex::index(grid->leafGridView().template begin<0>()->type()) << std::endl;
      testleafgridfunction<true>(grid->leafGridView());
    }

    {
      std::cout << "2D tests" << std::endl;
      // need a grid in order to test grid functions
      // typedef Dune::YaspGrid<2> Grid;
      typedef Dune::ALUGrid<2,2,Dune::simplex,Dune::conforming> Grid;
      Dune::FieldVector<double,2> l(0.0);
      Dune::FieldVector<double,2> u(1.0);
      std::array<unsigned int,2> N = {{1,1}};
      std::shared_ptr<Grid> grid = Dune::StructuredGridFactory<Grid>::createSimplexGrid(l,u,N);
      grid->globalRefine(1);

      std::cout << Dune::GlobalGeometryTypeIndex::index(grid->leafGridView().template begin<0>()->type()) << std::endl;
      testleafgridfunction<false>(grid->leafGridView());
    }

    // 3D
    {
      std::cout << "3D tests" << std::endl;
      // need a grid in order to test grid functions
      Dune::FieldVector<double,3> L(1.0);
      std::array<int,3> N(Dune::fill_array<int,3>(1));
      Dune::YaspGrid<3> grid(L,N);
      grid.globalRefine(1);

      testleafgridfunction<true>(grid.leafGridView());
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
