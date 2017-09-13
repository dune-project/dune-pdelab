// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/pdelab/finiteelementmap/p0fem.hh>
#include <dune/pdelab/finiteelementmap/pkfem.hh>
#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/powergridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/compositegridfunctionspace.hh>
#include <dune/pdelab/ordering/interleavedordering.hh>
#include <dune/pdelab/ordering/permutedordering.hh>

#include <dune/pdelab/backend/istl.hh>

template<typename GFS>
void check_ordering_reference(const GFS& gfs)
{
    const typename GFS::Ordering& ordering = gfs.ordering();

    Dune::PDELab::LocalFunctionSpace<GFS> lfs(gfs);

    typedef typename GFS::Traits::GridView GV;

    typename GV::template Codim<0>::Iterator it = gfs.gridView().template begin<0>();

    lfs.bind(*it);

    for (unsigned i = 0; i < lfs.size(); ++i)
      {
        const typename GFS::Ordering::Traits::DOFIndex& di = lfs.dofIndex(i);
        typename GFS::Ordering::Traits::ContainerIndex ci;
        ordering.mapIndex(di.view(),ci);
        std::cout << di << "    " << ci << std::endl;
      }
    using V = Dune::PDELab::Backend::Vector<GFS,double>;
    V x(gfs);
    x = 0.0;
    std::cout << std::endl;
}

template<typename GFS>
void check_ordering(GFS gfs)
{
    gfs.update();
    check_ordering_reference(gfs);
}

template<typename GFS>
void check_permuted_ordering(GFS gfs)
{
  gfs.update();
  check_ordering_reference(gfs);
  std::size_t i = gfs.blockCount();
  auto& permutation = gfs.orderingTag().permutation();
  permutation.resize(i);
  std::generate(permutation.begin(),permutation.end(),[&](){ return --i; });
  check_ordering_reference(gfs);
}

template<typename GFS>
void check_double_permuted_ordering(GFS gfs)
{
  gfs.update();
  std::size_t i = gfs.blockCount();
  auto& permutation = gfs.orderingTag().template permuted<1>().permutation();
  permutation.resize(i);
  std::generate(permutation.begin(),permutation.end(),[&](){ return --i; });
  gfs.orderingTag().template permuted<2>().permutation() = permutation;
  check_ordering_reference(gfs);
}


template<class GV>
static void testpermutedordering(const GV& gv)
{
  // instantiate finite element maps
  auto gt = Dune::GeometryTypes::quadrilateral;
  typedef Dune::PDELab::P0LocalFiniteElementMap<float,double,GV::dimension> P0FEM;
  P0FEM p0fem(gt);
  typedef Dune::PDELab::QkLocalFiniteElementMap<GV,float,double,1> Q12DFEM;
  Q12DFEM q12dfem(gv);
  typedef Dune::PDELab::QkLocalFiniteElementMap<GV,float,double,2> Q22DFEM;
  Q22DFEM q22dfem(gv);

  typedef Dune::PDELab::NoConstraints CON;

  typedef Dune::PDELab::ISTL::VectorBackend<> VBE;

  // make a grid function space
  typedef Dune::PDELab::GridFunctionSpace<GV,P0FEM,CON,VBE> P0GFS;
  P0GFS p0gfs(gv,p0fem);
  typedef Dune::PDELab::GridFunctionSpace<GV,Q12DFEM,CON,VBE> GFS1;
  GFS1 gfs1(gv,q12dfem);
  typedef Dune::PDELab::GridFunctionSpace<GV,Q22DFEM,CON,VBE> GFS2;
  GFS2 gfs2(gv,q22dfem);

  typedef Dune::PDELab::PowerGridFunctionSpace<GFS1,3,VBE,Dune::PDELab::InterleavedOrderingTag> P1GFS;

  {
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


  typedef Dune::PDELab::PowerGridFunctionSpace<GFS1,3,VBE,Dune::PDELab::EntityBlockedOrderingTag> EBPGFS1;
  EBPGFS1 ebpgfs1(gfs1);

  check_ordering(ebpgfs1);

  // Permuting a PowerGFS with lexicographic ordering

  typedef Dune::PDELab::PowerGridFunctionSpace<
    GFS1,
    3,
    VBE,
    Dune::PDELab::ordering::Permuted<Dune::PDELab::LexicographicOrderingTag>
    > PLexPGFS1;
  PLexPGFS1 plexpgfs1(gfs1);

  check_permuted_ordering(plexpgfs1);

  // Double permutation of a PowerGFS with lexicographic ordering (should yield original order)

  typedef Dune::PDELab::PowerGridFunctionSpace<
    GFS1,
    3,
    VBE,
    Dune::PDELab::ordering::Permuted<Dune::PDELab::ordering::Permuted<Dune::PDELab::LexicographicOrderingTag> >
    > PPLexPGFS1;
  PPLexPGFS1 pplexpgfs1(gfs1);

  check_double_permuted_ordering(pplexpgfs1);

  // Permuting a PowerGFS with entity-blocked ordering

  typedef Dune::PDELab::PowerGridFunctionSpace<
    GFS1,
    3,
    VBE,
    Dune::PDELab::ordering::Permuted<Dune::PDELab::EntityBlockedOrderingTag>
    > PEBPGFS1;
  PEBPGFS1 pebpgfs1(gfs1);

  check_permuted_ordering(pebpgfs1);

  // Double permutation of a PowerGFS with entity-blocked ordering (should yield original order)

  typedef Dune::PDELab::PowerGridFunctionSpace<
    GFS1,
    3,
    VBE,
    Dune::PDELab::ordering::Permuted<Dune::PDELab::ordering::Permuted<Dune::PDELab::EntityBlockedOrderingTag> >
    > PPEBPGFS1;
  PPEBPGFS1 ppebpgfs1(gfs1);

  check_double_permuted_ordering(ppebpgfs1);

  // Permuting a CompositeGFS with entity-blocked ordering

  typedef Dune::PDELab::CompositeGridFunctionSpace<
    VBE,
    Dune::PDELab::EntityBlockedOrderingTag,
    GFS1,
    EBPGFS1,
    GFS1
    > EBCGFS1;

  EBCGFS1 ebcgfs1(gfs1,ebpgfs1,gfs1);

  check_ordering(ebcgfs1);

  // Double permutation of a CompositeGFS with entity-blocked ordering (should yield original order)

  typedef Dune::PDELab::CompositeGridFunctionSpace<
    VBE,
    Dune::PDELab::ordering::Permuted<Dune::PDELab::EntityBlockedOrderingTag>,
    GFS1,
    EBPGFS1,
    GFS1
    > PEBCGFS1;

  PEBCGFS1 pebcgfs1(gfs1,ebpgfs1,gfs1);

  check_permuted_ordering(pebcgfs1);

  typedef Dune::PDELab::CompositeGridFunctionSpace<
    VBE,
    Dune::PDELab::ordering::Permuted<Dune::PDELab::ordering::Permuted<Dune::PDELab::EntityBlockedOrderingTag> >,
    GFS1,
    EBPGFS1,
    GFS1
    > PPEBCGFS1;

  PPEBCGFS1 ppebcgfs1(gfs1,ebpgfs1,gfs1);

  check_double_permuted_ordering(ppebcgfs1);



  typedef Dune::PDELab::CompositeGridFunctionSpace<
    VBE,
    Dune::PDELab::LexicographicOrderingTag,
    GFS1,
    EBPGFS1,
    GFS1
    > LexCGFS1;

  LexCGFS1 lexcgfs1(gfs1,ebpgfs1,gfs1);

  check_ordering(lexcgfs1);


  typedef Dune::PDELab::CompositeGridFunctionSpace<
    VBE,
    Dune::PDELab::ordering::Permuted<Dune::PDELab::LexicographicOrderingTag>,
    GFS1,
    EBPGFS1,
    GFS1
    > PLexCGFS1;

  PLexCGFS1 plexcgfs1(gfs1,ebpgfs1,gfs1);

  check_permuted_ordering(plexcgfs1);

  typedef Dune::PDELab::CompositeGridFunctionSpace<
    VBE,
    Dune::PDELab::ordering::Permuted<Dune::PDELab::ordering::Permuted<Dune::PDELab::LexicographicOrderingTag> >,
    GFS1,
    EBPGFS1,
    GFS1
    > PPLexCGFS1;

  PPLexCGFS1 pplexcgfs1(gfs1,ebpgfs1,gfs1);

  check_double_permuted_ordering(pplexcgfs1);

  // we need to manually trigger the update here
  pplexcgfs1.update();

  std::cout << pplexcgfs1.size() << " "
            << pplexcgfs1.blockCount() << " "
            << pplexcgfs1.maxLocalSize() << " "
            << pplexcgfs1.globalSize() << std::endl;

  std::cout << pplexcgfs1.template child<0>().size() << " "
            << pplexcgfs1.template child<0>().blockCount() << " "
            << pplexcgfs1.template child<0>().maxLocalSize() << " "
            << pplexcgfs1.template child<0>().globalSize() << std::endl;

  std::cout << pplexcgfs1.template child<1>().size() << " "
            << pplexcgfs1.template child<1>().blockCount() << " "
            << pplexcgfs1.template child<1>().maxLocalSize() << " "
            << pplexcgfs1.template child<1>().globalSize() << std::endl;

  std::cout << pplexcgfs1.template child<2>().size() << " "
            << pplexcgfs1.template child<2>().blockCount() << " "
            << pplexcgfs1.template child<2>().maxLocalSize() << " "
            << pplexcgfs1.template child<2>().globalSize() << std::endl;
}


int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper::instance(argc, argv);

    // 2D
    {
      std::cout << "2D tests" << std::endl;
      // need a grid in order to test Orderings
      Dune::FieldVector<double,2> L(1.0);
      std::array<int,2> N(Dune::fill_array<int,2>(1));
      Dune::YaspGrid<2> grid(L,N);
      grid.globalRefine(1);

      testpermutedordering(grid.leafGridView());
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
