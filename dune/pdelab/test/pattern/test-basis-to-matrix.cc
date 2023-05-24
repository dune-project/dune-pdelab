#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/pdelab/pattern/basis_to_pattern.hh>
#include <dune/pdelab/pattern/pattern_to_matrix.hh>
#include <dune/pdelab/pattern/sparsity_pattern.hh>

#include <dune/pdelab/common/tree_traversal.hh>

#include <dune/pdelab/operator/local_assembly/archetype.hh>

#include <dune/istl/io.hh>

#include <gtest/gtest.h>

#include <random>
#include <algorithm>

#include "../basis/fixture.hh"

// declare a dummy class that forwards it's template argument as the fixture to test
template <typename Fixture> class ForwardFixture : public Fixture {};

// declare the templated fixture
TYPED_TEST_SUITE_P(ForwardFixture);

void makeFullPattern(const Dune::PDELab::Concept::LocalBasis auto& ltrial,
                     const Dune::PDELab::Concept::LocalBasis auto& ltest,
                     auto& pattern)
{
  Dune::PDELab::forEachLeafNode(ltest.tree(), [&](auto& ltest_leaf, auto){
    Dune::PDELab::forEachLeafNode(ltrial.tree(), [&](auto& ltrial_leaf, auto){
      for (size_t i=0; i<ltest_leaf.size(); ++i)
        for (size_t j=0; j<ltrial_leaf.size(); ++j)
          pattern.addLink(ltest_leaf,i,ltrial_leaf,j);
    });
  });
}


class FullVolumePattern : public Dune::PDELab::LocalAssembly::Archetype
{
public:
  static constexpr bool localAssembleDoVolume() {return true;}

  // define sparsity pattern of operator representation
  static void localAssemblePatternVolume(const Dune::PDELab::Concept::LocalBasis auto& ltrial,
                                         const Dune::PDELab::Concept::LocalBasis auto& ltest,
                                         auto& pattern)
  {
    makeFullPattern(ltrial, ltest, pattern);
  }
};


TYPED_TEST_P(ForwardFixture, Test2DQ1VolumePattern) {

  auto basis = this->makeFixtureBasis(this->_grid->leafGridView());
  using Basis = decltype(basis);
  using Pattern = Dune::PDELab::LeafSparsePattern<Basis, Basis>;
  Pattern pattern{basis, {}, basis, {}, basis.localView().maxSize()*4};

  FullVolumePattern lop;

  Dune::BCRSMatrix<int> matrix;
  Dune::PDELab::basisToPattern(lop, pattern);
  Dune::PDELab::patternToMatrix(pattern, matrix);
  matrix = 0;

  std::fstream file(this->basisName()+"_volume_pattern.svg", std::ios_base::out);
  Dune::writeSVGMatrix(matrix, file);
}


class FullSkeletonPattern : public Dune::PDELab::LocalAssembly::Archetype
{
public:
  static constexpr bool localAssembleDoSkeleton() {return true;}

  // define sparsity pattern connecting self and neighbor dofs
  static void localAssemblePatternSkeleton(const Dune::Concept::Intersection auto& intersection,
                                           const Dune::PDELab::Concept::LocalBasis auto& ltrial_s, const Dune::PDELab::Concept::LocalBasis auto& ltest_s,
                                           const Dune::PDELab::Concept::LocalBasis auto& ltrial_n, const Dune::PDELab::Concept::LocalBasis auto& ltest_n,
                                           auto& pattern_ss, auto& pattern_sn,
                                           auto& pattern_ns, auto& pattern_nn)
  {
    makeFullPattern(ltrial_s, ltest_n, pattern_ns);
    makeFullPattern(ltrial_n, ltest_s, pattern_sn);
  }
};


TYPED_TEST_P(ForwardFixture, Test2DQ1SkeletonPattern) {
  auto basis = this->makeFixtureBasis(this->_grid->leafGridView());
  using Basis = decltype(basis);
  using Pattern = Dune::PDELab::LeafSparsePattern<Basis, Basis>;
  Pattern pattern{basis, {}, basis, {}, basis.localView().maxSize()*4 };

  FullSkeletonPattern lop;

  Dune::BCRSMatrix<int> matrix;
  Dune::PDELab::basisToPattern(lop, pattern);
  Dune::PDELab::patternToMatrix(pattern, matrix);
  matrix = 0;

  std::fstream file( this->basisName()+"_skeleton_pattern.svg", std::ios_base::out);
  Dune::writeSVGMatrix(matrix, file);
}

class FullBoundaryPattern : public Dune::PDELab::LocalAssembly::Archetype
{
public:
  static constexpr bool localAssembleDoBoundary() {return true;}

  // define sparsity pattern connecting dofs on boundary elements
  static void localAssemblePatternBoundary(const Dune::Concept::Intersection auto& intersection,
                                           const Dune::PDELab::Concept::LocalBasis auto& ltrial_s,
                                           const Dune::PDELab::Concept::LocalBasis auto& ltest_s,
                                           auto& pattern_ss)
  {
    makeFullPattern(ltrial_s, ltrial_s, pattern_ss);
  }
};

TYPED_TEST_P(ForwardFixture, Test2DQ1BoundaryPattern) {
  auto basis = this->makeFixtureBasis(this->_grid->leafGridView());
  using Basis = decltype(basis);
  using Pattern = Dune::PDELab::LeafSparsePattern<Basis, Basis>;
  Pattern pattern{basis, {}, basis, {}, basis.localView().maxSize()*4 };

  FullBoundaryPattern lop;

  Dune::BCRSMatrix<int> matrix;
  Dune::PDELab::basisToPattern(lop, pattern);
  Dune::PDELab::patternToMatrix(pattern, matrix);
  matrix = 0;

  std::fstream file(this->basisName()+"_boundary_pattern.svg", std::ios_base::out);
  Dune::writeSVGMatrix(matrix, file);
}

REGISTER_TYPED_TEST_SUITE_P(ForwardFixture,
    Test2DQ1VolumePattern,
    Test2DQ1SkeletonPattern,
    Test2DQ1BoundaryPattern);

// make a list of basis fictures with flat orderings
using FlatBasisFixtures = ::testing::Types<
  Basis1DQ1FlatFixture,
  Basis2DQ1FlatFixture,
  Basis2DQ2FlatFixture,
  Basis2DQ1DgFlatFixture,
  Basis2DQ2DgFlatFixture,
  Basis2DQ3DgFlatFixture,
  Basis2DQ2FlatArrayFlatByEntityFixture,
  Basis2DQ2FlatArrayFlatLexicographicFixture,
  Basis2DQ2FlatVectorFlatByEntityFixture,
  Basis2DQ2FlatVectorFlatLexicographicFixture,
  Basis2DQ2FlatTupleFlatByEntityFixture,
  Basis2DQ2FlatTupleFlatLexicographicFixture
>;

// instantiate the test for flat orderings
INSTANTIATE_TYPED_TEST_SUITE_P(Pattern, ForwardFixture, FlatBasisFixtures);


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  Dune::MPIHelper::instance(argc, argv);
  return RUN_ALL_TESTS();
}
