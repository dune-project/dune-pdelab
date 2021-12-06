// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include <memory>

#include <dune/common/exceptions.hh>
#include <dune/common/filledarray.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/yaspgrid.hh>
#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif
#include <dune/grid/common/scsgmapper.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/pdelab.hh>

template<typename GFS>
void check_ordering(const GFS& gfs)
{
    const auto& ordering = gfs.ordering();

    Dune::PDELab::LocalFunctionSpace lfs{gfs};
    Dune::PDELab::LFSIndexCache lfs_cache{lfs};

    for (const auto& element : elements(gfs.entitySet()))
    {
      // lfs computes dof indices
      lfs.bind(element);
      // cache computes container indices
      lfs_cache.update();

      for (unsigned i = 0; i < lfs.size(); ++i) {
        // get i-th local dof index
        auto di = lfs.dofIndex(i);
        // map dof index to a container index (from ordering)
        auto ci_map = ordering.mapIndex(di);
        // map local index to a container index (from cache)
        auto ci_cache_i = lfs_cache.containerIndex(i);
        // map dof index to a container index (from cache)
        auto ci_cache_di = lfs_cache.containerIndex(di);
        if (ci_map != ci_cache_i or ci_map != ci_cache_di)
          DUNE_THROW(Dune::RangeError, "Container index mappings do not match");

        // check that all size suffixes fit the current container index
        auto size_suffix = ci_map;
        while (size_suffix.size() != 0) {
          // get outer container index in the suffix
          auto block_index = size_suffix.back();
          // calculate the size for a container that would hold such block
          size_suffix.pop_back();
          auto size = ordering.size(size_suffix);
          // the index should always fit into the size
          if (not (size > block_index))
            DUNE_THROW(Dune::RangeError, "Size `" << size << "` for CI suffix '"
                                   << size_suffix
                                   << "' is not big enough to hold sub-block '"
                                   << block_index << "'");
        }
      }
    }

    using V = Dune::PDELab::Backend::Vector<GFS,double>;
    V x(gfs);
    x = 0.0;
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
    MonomFEM monomfem(cellmapper, gt, 2);

    typedef Dune::PDELab::NoConstraints CON;

    typedef Dune::PDELab::ISTL::VectorBackend<> VBE;

    // make a grid function space
    typedef Dune::PDELab::GridFunctionSpace<GV,P0FEM,CON,VBE> P0GFS;
    P0GFS p0gfs(gv,p0fem);
    typedef Dune::PDELab::GridFunctionSpace<GV,Q12DFEM,CON,VBE> GFS1;
    {
      GFS1 gfs1(gv,q12dfem);
      check_ordering(gfs1);
    }
    GFS1 gfs1(gv,q12dfem);
    typedef Dune::PDELab::GridFunctionSpace<GV,Q22DFEM,CON,VBE> GFS2;
    GFS2 gfs2(gv,q22dfem);
    typedef Dune::PDELab::GridFunctionSpace<GV,MonomFEM,CON,VBE> GFS3;
    {
      GFS3 gfs3(gv,monomfem);
      check_ordering(gfs3);
    }
    GFS3 gfs3(gv,monomfem);

    typedef Dune::PDELab::PowerGridFunctionSpace<GFS1,3,VBE,Dune::PDELab::InterleavedOrderingTag> P1GFS;

    {
      // keep this block around to test the vector-based constructor of the ordering tag
      std::vector<std::size_t> p1_gfs_block_sizes(3);
      std::fill(p1_gfs_block_sizes.begin(),p1_gfs_block_sizes.end(),1);
      P1GFS p1gfs(gfs1,gfs1,gfs1,VBE(),p1_gfs_block_sizes);
    }

    P1GFS p1gfs(gfs1,gfs1,gfs1,VBE(),{{1,1,1}});

    typedef Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed> NVBE;

    using EGFS  = Dune::PDELab::GridFunctionSpace<GV,Q12DFEM,CON,VBE>;
    EGFS egfs(gv,q12dfem);
    using EPGFS = Dune::PDELab::PowerGridFunctionSpace<EGFS,3,VBE,Dune::PDELab::InterleavedOrderingTag>;
    {
      EPGFS epgfs(egfs,egfs,egfs,VBE(),{{1,1,1}});
      check_ordering(epgfs);
    }

    EPGFS epgfs(egfs,egfs,egfs,VBE(),{{1,1,1}});
    typedef Dune::PDELab::PowerGridFunctionSpace<EPGFS,2,NVBE,Dune::PDELab::InterleavedOrderingTag> PGFS;
    std::vector<std::size_t> p_gfs_block_sizes(2);
    std::fill(p_gfs_block_sizes.begin(),p_gfs_block_sizes.end(),3);
    PGFS pgfs(epgfs,epgfs,NVBE(),p_gfs_block_sizes);


    typedef Dune::PDELab::PowerGridFunctionSpace<GFS1,3,VBE,Dune::PDELab::EntityBlockedOrderingTag> P2GFS;
    P2GFS p2gfs(gfs1,gfs1,gfs1);

    typedef Dune::PDELab::PowerGridFunctionSpace<P2GFS,2,NVBE,Dune::PDELab::EntityBlockedOrderingTag> P3GFS;
    P3GFS p3gfs(p2gfs,p2gfs,NVBE());

    check_ordering(pgfs);
    //check_ordering(p3gfs);

    typedef Dune::PDELab::PowerGridFunctionSpace<GFS1,3,VBE,Dune::PDELab::EntityBlockedOrderingTag> EBPGFS1;
    {
      EBPGFS1 ebpgfs1(gfs1);
      check_ordering(ebpgfs1);
    }
    EBPGFS1 ebpgfs1(gfs1);

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
    check_ordering(p0gfs);
// Doesn't work, we need a grid with triangular elemets for that
//  typedef Dune::PDELab::GridFunctionSpace<GV,P1FEM> P1GFS;
//  P1GFS p1gfs(gv,p1fem);
    typedef Dune::PDELab::GridFunctionSpace<GV,Q1FEM> Q1GFS;
    Q1GFS q1gfs(gv,q1fem);
    check_ordering(q1gfs);
  }
};

template<bool cube, class GV>
void testleafgridfunction(const GV& gv)
{
  test<GV::dimension,cube>::testleafgridfunction(gv);
}


int main(int argc, char** argv)
{
  //Maybe initialize Mpi
  Dune::MPIHelper::instance(argc, argv);

#if HAVE_DUNE_ALUGRID
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
#endif


  // 2D
  {
    std::cout << "2D tests" << std::endl;
    // need a grid in order to test grid functions
    Dune::FieldVector<double,2> L(1.0);
    std::array<int,2> N(Dune::filledArray<2,int>(1));
    Dune::YaspGrid<2> grid(L,N);
    grid.globalRefine(1);

    testleafgridfunction<true>(grid.leafGridView());
  }

  // 3D
  {
    std::cout << "3D tests" << std::endl;
    // need a grid in order to test grid functions
    Dune::FieldVector<double,3> L(1.0);
    std::array<int,3> N(Dune::filledArray<3,int>(1));
    Dune::YaspGrid<3> grid(L,N);
    grid.globalRefine(1);

    testleafgridfunction<true>(grid.leafGridView());
  }

  // test passed
  return 0;
}
