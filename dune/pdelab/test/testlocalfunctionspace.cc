// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include <vector>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/localvector.hh>
#include <dune/pdelab/gridfunctionspace/subspace.hh>


// test function trees
template<class GV>
void test (const GV& gv)
{
  // instantiate finite element maps
  typedef Dune::PDELab::QkLocalFiniteElementMap<GV,float,double,2> Q22DFEM;
  Q22DFEM q22dfem(gv);
  typedef Dune::PDELab::QkLocalFiniteElementMap<GV,float,double,1> Q12DFEM;
  Q12DFEM q12dfem(gv);

  // make a grid function space
  typedef Dune::PDELab::GridFunctionSpace<GV,Q22DFEM> Q2GFS;
  Q2GFS q2gfs(gv,q22dfem);
  typedef Dune::PDELab::GridFunctionSpace<GV,Q12DFEM> Q1GFS;
  Q1GFS q1gfs(gv,q12dfem);

  // power grid function space
  typedef Dune::PDELab::PowerGridFunctionSpace<Q2GFS,2,
    Dune::PDELab::ISTL::VectorBackend<>, Dune::PDELab::LexicographicOrderingTag> PowerGFS;
  PowerGFS powergfs(q2gfs);

  // composite grid function space
  typedef Dune::PDELab::CompositeGridFunctionSpace<Dune::PDELab::ISTL::VectorBackend<>,
      Dune::PDELab::LexicographicOrderingTag,PowerGFS,Q1GFS> CompositeGFS;
  CompositeGFS compositegfs(powergfs,q1gfs);

  // make coefficent Vectors - we need to make copies of the spaces because we stuck
  // them in a hierarchy
  using V = Dune::PDELab::Backend::Vector<Q2GFS,double>;
  Q2GFS q2gfs2(gv,q22dfem);
  V x(q2gfs2);
  x = 0.0;
  using VP = Dune::PDELab::Backend::Vector<PowerGFS,double>;
  Q2GFS q2gfs_pc(gv,q22dfem);
  PowerGFS powergfs2(q2gfs_pc);
  VP xp(powergfs2);
  xp = 0.0;

  // make local function space object
  typedef Dune::PDELab::AnySpaceTag Tag;
  typedef typename Dune::PDELab::LocalFunctionSpace<Q2GFS> Q2LFS;
  Q2LFS q2lfs(q2gfs2);
  typedef Dune::PDELab::LFSIndexCache<Q2LFS> Q2LFSCache;
  Q2LFSCache q2lfsCache(q2lfs);
  typedef typename V::template ConstLocalView<Q2LFSCache> VView;
  VView x_view(x);
  Dune::PDELab::LocalVector<double, Tag> xl(q2lfs.maxSize());

  typedef typename Dune::PDELab::LocalFunctionSpace<PowerGFS> PowerLFS;
  PowerLFS powerlfs(powergfs2);
  typedef Dune::PDELab::LFSIndexCache<PowerLFS> PowerLFSCache;
  PowerLFSCache powerlfsCache(powerlfs);
  typedef typename VP::template ConstLocalView<PowerLFSCache> VPView;
  VPView xp_view(xp);
  Dune::PDELab::LocalVector<double, Tag> xlp(powerlfs.maxSize());

  typedef typename Dune::PDELab::LocalFunctionSpace<CompositeGFS> CompositeLFS;
  CompositeLFS compositelfs(compositegfs);
  //typedef Dune::PDELab::LFSIndexCache<CompositeLFS> CompositeLFSCache;
  //CompositeLFSCache compositelfsCache(compositelfs);
  //  std::vector<double> xlc(compositelfs.maxSize());

  typedef Dune::TypeTree::TreePath<1> Path1;
  typedef Dune::PDELab::GridFunctionSubSpace<CompositeGFS, Path1> SubGFS1;
  typedef Dune::PDELab::LocalFunctionSpace<SubGFS1> SubLFS1;
  SubGFS1 subgfs1(compositegfs);
  SubLFS1 sublfs1(subgfs1);

  // loop over elements
  typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
  for (ElementIterator it = gv.template begin<0>();
       it!=gv.template end<0>(); ++it)
    {
      q2lfs.bind(*it);
      q2lfs.debug();
      q2lfsCache.update();
      x_view.bind(q2lfsCache);
      x_view.read(xl);
      x_view.unbind();
      assert(q2lfs.size() ==
          q2lfs.localVectorSize());

      powerlfs.bind(*it);
      powerlfs.debug();
      powerlfsCache.update();
      xp_view.bind(powerlfsCache);
      xp_view.read(xlp);
      assert(powerlfs.size() ==
          powerlfs.localVectorSize());
      assert(powerlfs.localVectorSize() ==
          powerlfs.template child<0>().localVectorSize());
      assert(powerlfs.localVectorSize() ==
          powerlfs.template child<1>().localVectorSize());

      compositelfs.bind(*it);
      compositelfs.debug();
      assert(compositelfs.size() ==
          compositelfs.localVectorSize());
      assert(compositelfs.localVectorSize() ==
          compositelfs.template child<0>().localVectorSize());
      assert(compositelfs.localVectorSize() ==
          compositelfs.template child<0>().template child<0>().localVectorSize());
      assert(compositelfs.localVectorSize() ==
          compositelfs.template child<0>().template child<1>().localVectorSize());
      assert(compositelfs.localVectorSize() ==
          compositelfs.template child<1>().localVectorSize());

      // check LFS<SubSpace<CompositeSpace,1>> == LFS<CompositeSpace>::Child<1>
      sublfs1.bind(*it);
      assert(compositelfs.template child<1>().size() == sublfs1.size());
      for (unsigned int i=0; i<compositelfs.template child<1>().size(); i++)
      {
        if ( compositelfs.template child<1>().dofIndex(i) != sublfs1.dofIndex(i) )
        {
          std::cout << "DOF " << i << ": \n"
                    << "LFS::CHILD\t" << compositelfs.template child<1>().dofIndex(i) << "\n"
                    << "LFS<SUB>\t" << sublfs1.dofIndex(i) << "\n";
        }
        assert( compositelfs.template child<1>().dofIndex(i) == sublfs1.dofIndex(i) );
      }

    }
}

int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper::instance(argc, argv);

    // need a grid in order to test grid functions
    Dune::FieldVector<double,2> L(1.0);
    std::array<int,2> N(Dune::fill_array<int,2>(1));
    Dune::YaspGrid<2> grid(L,N);
    grid.globalRefine(1);

    test(grid.leafGridView());

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
