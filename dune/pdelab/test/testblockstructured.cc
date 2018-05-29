// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <vector>

#include <dune/grid/yaspgrid.hh>
#include <dune/pdelab/finiteelementmap/blockstructuredqkfem.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/blockstructured/localfunctionspace.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/backend/istl/bcrsmatrixbackend.hh>
#include <dune/pdelab/backend/simple/vector.hh>
#include <dune/pdelab/backend/istl/vector.hh>
#include <dune/pdelab/gridoperator/blockstructured.hh>
#include <dune/pdelab/constraints/conforming.hh>

template<typename GV>
struct TestData {
  using Q1_FEM = Dune::PDELab::BlockstructuredQkLocalFiniteElementMap<GV, double, double, 1, 4>;
  using Q2_FEM = Dune::PDELab::BlockstructuredQkLocalFiniteElementMap<GV, double, double, 2, 4>;
  using Backend = Dune::PDELab::ISTL::VectorBackend<>;
  using Q1_LeafGFS = Dune::PDELab::GridFunctionSpace<GV, Q1_FEM, Dune::PDELab::NoConstraints, Backend>;
  using Q2_LeafGFS = Dune::PDELab::GridFunctionSpace<GV, Q2_FEM, Dune::PDELab::NoConstraints, Backend>;
  using PowerGFS = Dune::PDELab::PowerGridFunctionSpace<Q1_LeafGFS, 2, Backend, Dune::PDELab::LexicographicOrderingTag>;
  using CompositeGFS = Dune::PDELab::CompositeGridFunctionSpace<Backend, Dune::PDELab::LexicographicOrderingTag, PowerGFS, Q2_LeafGFS>;

  Q1_FEM fem;
  Q2_FEM _q2_fem;
  std::shared_ptr<Q1_LeafGFS> pLeafGFS;
  Q1_LeafGFS _q1_gfs;
  Q2_LeafGFS _q2_gfs;
  PowerGFS _pgfs;
  std::shared_ptr<CompositeGFS> pCompositeGFS;

  TestData(const GV& gv)
      : fem(gv), _q2_fem(gv), pLeafGFS(std::make_shared<Q1_LeafGFS>(gv, fem)), _q1_gfs(gv, fem), _q2_gfs(gv, _q2_fem),
        _pgfs(_q1_gfs), pCompositeGFS(std::make_shared<CompositeGFS>(_pgfs, _q2_gfs)) {}
};

class LocalOperator
    : public Dune::PDELab::LocalOperatorDefaultFlags {
public:
  enum {
    doAlphaVolume = true
  };

  template<typename LFSU, typename R, typename LFSV, typename X, typename EG>
  void alpha_volume(const EG &eg, const LFSU &lfsu, const X &x, const LFSV &lfsv, R &r) const {
    using namespace Dune::Indices;

    auto lfsu_0 = child(lfsu, _0);
    auto lfsu_0_0 = child(lfsu_0, _0);
    auto lfsu_0_1 = child(lfsu_0, _1);

    auto lfsu_1 = child(lfsu, _1);

    for (std::size_t i = 0; i < lfsu_0_0.size(); ++i) {
      r.accumulate(lfsu_0_0, i, 0);
    }
    for (std::size_t i = 0; i < lfsu_0_1.size(); ++i) {
      r.accumulate(lfsu_0_1, i, 1);
    }
    for (std::size_t i = 0; i < lfsu_1.size(); ++i) {
      r.accumulate(lfsu_1, i, 2);
    }
  }
};

template<typename LFS, typename PGFS>
auto createAndBindLFS(PGFS pgfs) {
  const auto element = *pgfs->gridView().template begin<0>();
  auto plfs = std::make_shared<LFS>(pgfs);
  plfs->bind(element);
  return plfs;
}

template<typename PGFS>
auto setupBlockstructuredLFS(PGFS pgfs) {
  using LFS = Dune::Blockstructured::LocalFunctionSpace<typename PGFS::element_type>;
  return createAndBindLFS<LFS>(pgfs);
}

template<typename PGFS>
auto setupPDELabLFS(PGFS pgfs) {
  using LFS = Dune::PDELab::LocalFunctionSpace<typename PGFS::element_type>;
  return createAndBindLFS<LFS>(pgfs);
}

template<typename TestData>
void testBlockstructuredLeafLFS(const TestData &td) {
  auto plfs = setupBlockstructuredLFS(td.pLeafGFS);

  auto lfs_indices = *plfs->_subentityWiseDOFs_ptr;
  assert(lfs_indices.size() == 1);

  auto pcompareLFS = setupPDELabLFS(td.pLeafGFS);

  auto coeffs = plfs->finiteElement().localCoefficients();
  for (std::size_t i = 0; i < pcompareLFS->size(); ++i) {
    auto localKey = coeffs.localKey(i);
    if(localKey.index() == 0) {
      assert(lfs_indices[0].index(localKey.subEntity(), localKey.codim()) == pcompareLFS->dofIndex(i));
    }
  }
}

template<typename TestData>
void testBlockstructuredTreeLFS(const TestData &td) {
  auto lfs_ptr = setupBlockstructuredLFS(td.pCompositeGFS);

  const auto& lfs_indices = *lfs_ptr->_subentityWiseDOFs_ptr;
  assert(lfs_indices.size() == 3);

  auto compareLFS_ptr = setupPDELabLFS(td.pCompositeGFS);

  Dune::TypeTree::forEachLeafNode(*compareLFS_ptr, [lfs_ptr, lfs_indices] (auto& Node, auto& TreePath)
  {
    auto coeffs = Node.finiteElement().localCoefficients();

    int leaf = Dune::TypeTree::child(*lfs_ptr, TreePath).offsetLeafs;

    for (std::size_t i = 0; i < coeffs.size(); ++i) {
      auto localKey = coeffs.localKey(i);

      if(localKey.index() == 0)
        assert(lfs_indices[leaf].index(localKey.subEntity(), localKey.codim()) == Node.dofIndex(i));
    }
  });
}

template<typename LFSC, typename LFS>
auto createAndUpdateLFSC(const LFS& lfs)
{
  auto plfsc = std::make_shared<LFSC>(lfs, Dune::PDELab::EmptyTransformation(), false);
  plfsc->update();
  return plfsc;
}

template<typename LFS>
auto setupBlockstructuredLFSC(const LFS &lfs){
  using LFSC = Dune::Blockstructured::LFSIndexCache<LFS, Dune::PDELab::EmptyTransformation>;
  return createAndUpdateLFSC<LFSC>(lfs);
}

template<typename LFS>
auto setupPDELabLFSC(const LFS& lfs){
  using LFSC = Dune::PDELab::LFSIndexCache<LFS, Dune::PDELab::EmptyTransformation>;
  return createAndUpdateLFSC<LFSC>(lfs);
}

template<typename TestData>
void testBlockstructuredLFSC(const TestData &td) {
  auto plfs = setupBlockstructuredLFS(td.pCompositeGFS);
  auto plfsc = setupBlockstructuredLFSC(*plfs);

  auto compare_plfs = setupPDELabLFS(td.pCompositeGFS);
  auto compare_plfsc = setupPDELabLFSC(*compare_plfs);

  Dune::TypeTree::forEachLeafNode(*compare_plfs, [compare_plfsc, plfs, plfsc] (auto& Node, auto& TreePath)
  {
    auto coeffs = Node.finiteElement().localCoefficients();

    int leaf = Dune::TypeTree::child(*plfs, TreePath).offsetLeafs;

    for (std::size_t i = 0; i < coeffs.size(); ++i) {
      auto localKey = coeffs.localKey(i);

      if(localKey.index() == 0) {
        assert(
            plfsc->containerIndex(leaf, localKey.subEntity(), localKey.codim()) ==
            compare_plfsc->containerIndex(Node.offset + i));
      }

      assert(plfsc->localIndex(leaf, localKey.subEntity(), localKey.codim(), localKey.index()) == Node.offset + i);
    }
  });
}

template<typename TestData>
void testBlockstructuredUncachedVectorView(const TestData &td) {
  auto plfs = setupBlockstructuredLFS(td.pCompositeGFS);
  auto plfsc = setupBlockstructuredLFSC(*plfs);

  std::vector<int> local_write_to(plfsc->size(), 0);
  std::vector<int> local_read_from(plfsc->size());

  for (auto &i: local_read_from)
    i = std::rand();

  using Container = Dune::PDELab::Simple::VectorContainer<typename TestData::CompositeGFS, std::vector<int>>;
  Container container(*td.pCompositeGFS, 0);

  Dune::PDELab::UncachedVectorView<Container, typename decltype(plfsc)::element_type> vectorView;

  vectorView.attach(container);
  vectorView.bind(*plfsc);

  vectorView.add(local_read_from);
  vectorView.read(local_write_to);

  for (std::size_t i = 0; i < vectorView.size(); ++i)
    assert(local_write_to[i] == local_read_from[i]);
}

template<typename TestData>
void testBlockstructuredGridOperator(const TestData &td) {
  using MatrixBackend = Dune::PDELab::ISTL::BCRSMatrixBackend<>;
  using GO = Dune::Blockstructured::BlockstructuredGridOperator<
      typename TestData::CompositeGFS, typename TestData::CompositeGFS,
      LocalOperator, MatrixBackend, double, double, double>;
  using GO_Compare = Dune::PDELab::GridOperator<
      typename TestData::CompositeGFS, typename TestData::CompositeGFS,
      LocalOperator, MatrixBackend, double, double, double>;

  LocalOperator lop;
  MatrixBackend mb(1);

  GO go(*td.pCompositeGFS, *td.pCompositeGFS, lop, mb);
  GO_Compare go_compare(*td.pCompositeGFS, *td.pCompositeGFS, lop, mb);

  Dune::PDELab::Backend::Vector<typename TestData::CompositeGFS, double> x(*td.pCompositeGFS, 0);
  Dune::PDELab::Backend::Vector<typename TestData::CompositeGFS, double> r(*td.pCompositeGFS, 0);
  Dune::PDELab::Backend::Vector<typename TestData::CompositeGFS, double> r_compare(*td.pCompositeGFS, 0);

  go.residual(x, r);
  go_compare.residual(x, r_compare);

  for (std::size_t i = 0; i < r.N(); ++i) {
    assert((*r.storage())[i] == (*r_compare.storage())[i]);
  }
}

template<typename TestData>
void runAllTests(const TestData& td){
  testBlockstructuredLeafLFS(td);
  testBlockstructuredTreeLFS(td);

  testBlockstructuredLFSC(td);

  testBlockstructuredUncachedVectorView(td);

  testBlockstructuredGridOperator(td);
}

int main(int argc, char **argv) {
  try {
    Dune::MPIHelper::instance(argc, argv);

    using Grid2d = Dune::YaspGrid<2>;
    using Grid3d = Dune::YaspGrid<3>;

    Grid2d grid2({1, 1}, {1, 1});
    Grid3d grid3({1, 1, 1}, {1, 1, 1});

    TestData<typename Grid2d::LeafGridView> td2(grid2.leafGridView());
    TestData<typename Grid3d::LeafGridView> td3(grid3.leafGridView());

    runAllTests(td2);
    runAllTests(td3);
  }
  catch (Dune::Exception &e) {
    std::cerr << "Dune reported error: " << e << std::endl;
    return 1;
  }
  catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
    return 1;
  }
}