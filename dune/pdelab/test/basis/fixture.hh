#ifndef DUNE_PDELAB_TEST_FIXTURE_BASIS_HH
#define DUNE_PDELAB_TEST_FIXTURE_BASIS_HH

#include <dune/pdelab/basis/basis.hh>

#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/pdelab/finiteelementmap/pkfem.hh>
#include <dune/pdelab/finiteelementmap/qkdg.hh>
#include <dune/pdelab/finiteelementmap/p0fem.hh>

#include "../fixture-grid.hh"

enum class MergingType {
  ByEntity,
  Lexicographic
};

std::string to_string(MergingType mode) {
  if (mode == MergingType::ByEntity)
    return "ByEntity";
  else if (mode == MergingType::Lexicographic)
    return "Lexicographic";
  else
    DUNE_THROW(Dune::RangeError, "");
}


enum class NodeType {
  Array,
  Vector,
  Tuple
};

std::string to_string(NodeType node_type) {
  if (node_type == NodeType::Array)
    return "Array";
  else if (node_type == NodeType::Array)
    return "Vector";
  else
    return "Tuple";
}


template<std::size_t degree, bool Blocked, class EntitySet>
auto makeQkPreBasis(const EntitySet& entity_set) {

  using QkFEM = Dune::PDELab::QkLocalFiniteElementMap<EntitySet, double, double, degree>;
  auto qkfem = std::make_shared<QkFEM>(entity_set);

  auto strategy = [&]{
    if constexpr (Blocked)
      return Dune::PDELab::Strategy::blockedByEntity(entity_set);
    else
      return Dune::PDELab::Strategy::flatByEntity(entity_set);
  }();
  return makePreBasis(strategy, qkfem);
}

template<std::size_t k, bool Blocked, class EntitySet>
auto makeQkDgBasis(const EntitySet& entity_set) {

  using DkDGFEM = Dune::PDELab::QkDGLocalFiniteElementMap<typename EntitySet::GridView::ctype,double,k,EntitySet::dimension>;
  auto qkfem = std::make_shared<DkDGFEM>();

  auto strategy = [&]{
    if constexpr (Blocked)
      return Dune::PDELab::Strategy::blockedByEntity(entity_set);
    else
      return Dune::PDELab::Strategy::flatByEntity(entity_set);
  }();
  return makePreBasis(strategy, qkfem);
}

template<MergingType merging_type, NodeType node_type, bool blocked, class PreBasis>
auto makeCompositePreBasis(const PreBasis& space) {
  auto strategy = [&]{
    if constexpr (merging_type == MergingType::ByEntity){
      auto entity_set = space.mergingStrategy().entitySet();
      if constexpr (blocked)
        return Dune::PDELab::Strategy::blockedByEntity(entity_set);
      else
        return Dune::PDELab::Strategy::flatByEntity(entity_set);
    } else {
      if constexpr (blocked)
        return Dune::PDELab::Strategy::blockedLexicographic();
      else
        return Dune::PDELab::Strategy::flatLexicographic();
    }
  }();
  if constexpr (node_type == NodeType::Array)
    return composite(strategy, std::array{space, space, space});
  else if constexpr (node_type == NodeType::Vector)
    return composite(strategy, std::vector{space, space, space});
  else
    return composite(strategy, std::tuple{space, space, space});
}


template<std::size_t dim, std::size_t degree, bool Blocked>
class QkFixture : public StructuredGridFixture<dim> {
 protected:
  static auto makeFixtureBasis(const auto& entity_set) {
    return makeBasis(entity_set, makeQkPreBasis<degree,Blocked>(entity_set));
  }
  static auto spaceName() {
    using std::to_string;
    return "Basis" + to_string(dim) + "DQ" + to_string(degree) + (Blocked ? "Blocked" : "Flat");
  }
};

template<std::size_t dim, std::size_t degree, bool Blocked>
class QkDgFixture : public StructuredGridFixture<dim> {
 protected:
  static auto makeFixtureBasis(const auto& entity_set) {
    return makeBasis(entity_set, makeQkDgBasis<degree,Blocked>(entity_set));
  }
  static auto spaceName() {
    using std::to_string;
    return "Basis" + to_string(dim) + "DQ" + to_string(degree) + "Dg" + (Blocked ? "Blocked" : "Flat");
  }
};

template<std::size_t dim, std::size_t degree, bool leaf_blocked, MergingType composite_merging, NodeType composite_node, bool composite_blocked>
class CompositeQkFixture : public StructuredGridFixture<dim> {
 protected:
  static auto makeFixtureBasis(const auto& entity_set) {
    auto space_qk = makeQkPreBasis<degree,leaf_blocked>(entity_set);
    return makeBasis(entity_set, makeCompositePreBasis<composite_merging, composite_node, composite_blocked>(space_qk));
  }

  static auto spaceName() {
    using std::to_string;
    return "Basis" + to_string(dim) + "DQ" + to_string(degree)
      + (leaf_blocked ? "Blocked" : "Flat")
      + to_string(composite_node)
      + (composite_blocked ? "Blocked" : "Flat")
      + to_string(composite_merging);
  }
};

template<std::size_t dim, std::size_t degree, bool leaf_blocked, MergingType composite_merging, NodeType composite_node, bool composite_blocked>
class CompositeQkDgFixture : public StructuredGridFixture<dim> {
 protected:
  static auto makeFixtureBasis(const auto& entity_set) {
    auto space_qk = makeQkDgBasis<degree,leaf_blocked>(entity_set);
    return makeBasis(entity_set, makeCompositePreBasis<composite_merging, composite_node, composite_blocked>(space_qk));
  }

  static auto spaceName() {
    using std::to_string;
    return "Basis" + to_string(dim) + "DQ" + to_string(degree)
      + "Dg"
      + (leaf_blocked ? "Blocked" : "Flat")
      + to_string(composite_node)
      + (composite_blocked ? "Blocked" : "Flat")
      + to_string(composite_merging);
  }
};

// these fixtures add a method `makeFixtureBasis` that generates an spaces out of an entity set

// q1/2 spaces
using Basis1DQ1FlatFixture = QkFixture<1,1,false>;
using Basis1DQ1BlockedFixture = QkFixture<1,1,true>;

using Basis2DQ1FlatFixture = QkFixture<2,1,false>;
using Basis2DQ2FlatFixture = QkFixture<2,2,false>;
using Basis2DQ1BlockedFixture = QkFixture<2,1,true>;
using Basis2DQ2BlockedFixture = QkFixture<2,2,true>;

// dg spaces
using Basis2DQ1DgFlatFixture = QkDgFixture<2,1,false>;
using Basis2DQ2DgFlatFixture = QkDgFixture<2,2,false>;
using Basis2DQ3DgFlatFixture = QkDgFixture<2,3,false>;
using Basis2DQ1DgBlockedFixture = QkDgFixture<2,1,true>;
using Basis2DQ2DgBlockedFixture = QkDgFixture<2,2,true>;
using Basis2DQ3DgBlockedFixture = QkDgFixture<2,3,true>;

// composite q12spaces array
using Basis2DQ2FlatArrayFlatByEntityFixture = CompositeQkFixture<2,2,false,MergingType::ByEntity,NodeType::Array,false>;
using Basis2DQ2FlatArrayBlockedByEntityFixture = CompositeQkFixture<2,2,false,MergingType::ByEntity,NodeType::Array,true>;
using Basis2DQ2FlatArrayFlatLexicographicFixture = CompositeQkFixture<2,2,false,MergingType::Lexicographic,NodeType::Array,false>;
using Basis2DQ2FlatArrayBlockedLexicographicFixture = CompositeQkFixture<2,2,false,MergingType::Lexicographic,NodeType::Array,true>;

using Basis2DQ2BlockedArrayFlatByEntityFixture = CompositeQkFixture<2,2,true,MergingType::ByEntity,NodeType::Array,false>;
using Basis2DQ2BlockedArrayBlockedByEntityFixture = CompositeQkFixture<2,2,true,MergingType::ByEntity,NodeType::Array,true>;
using Basis2DQ2BlockedArrayFlatLexicographicFixture = CompositeQkFixture<2,2,true,MergingType::Lexicographic,NodeType::Array,false>;
using Basis2DQ2BlockedArrayBlockedLexicographicFixture = CompositeQkFixture<2,2,true,MergingType::Lexicographic,NodeType::Array,true>;

// composite q12spaces vector
using Basis2DQ2FlatVectorFlatByEntityFixture = CompositeQkFixture<2,2,false,MergingType::ByEntity,NodeType::Vector,false>;
using Basis2DQ2FlatVectorBlockedByEntityFixture = CompositeQkFixture<2,2,false,MergingType::ByEntity,NodeType::Vector,true>;
using Basis2DQ2FlatVectorFlatLexicographicFixture = CompositeQkFixture<2,2,false,MergingType::Lexicographic,NodeType::Vector,false>;
using Basis2DQ2FlatVectorBlockedLexicographicFixture = CompositeQkFixture<2,2,false,MergingType::Lexicographic,NodeType::Vector,true>;

using Basis2DQ2BlockedVectorFlatByEntityFixture = CompositeQkFixture<2,2,true,MergingType::ByEntity,NodeType::Vector,false>;
using Basis2DQ2BlockedVectorBlockedByEntityFixture = CompositeQkFixture<2,2,true,MergingType::ByEntity,NodeType::Vector,true>;
using Basis2DQ2BlockedVectorFlatLexicographicFixture = CompositeQkFixture<2,2,true,MergingType::Lexicographic,NodeType::Vector,false>;
using Basis2DQ2BlockedVectorBlockedLexicographicFixture = CompositeQkFixture<2,2,true,MergingType::Lexicographic,NodeType::Vector,true>;

// composite q12spaces tuple
using Basis2DQ2FlatTupleFlatByEntityFixture = CompositeQkFixture<2,2,false,MergingType::ByEntity,NodeType::Tuple,false>;
using Basis2DQ2FlatTupleBlockedByEntityFixture = CompositeQkFixture<2,2,false,MergingType::ByEntity,NodeType::Tuple,true>;
using Basis2DQ2FlatTupleFlatLexicographicFixture = CompositeQkFixture<2,2,false,MergingType::Lexicographic,NodeType::Tuple,false>;
using Basis2DQ2FlatTupleBlockedLexicographicFixture = CompositeQkFixture<2,2,false,MergingType::Lexicographic,NodeType::Tuple,true>;

using Basis2DQ2BlockedTupleFlatByEntityFixture = CompositeQkFixture<2,2,true,MergingType::ByEntity,NodeType::Tuple,false>;
using Basis2DQ2BlockedTupleBlockedByEntityFixture = CompositeQkFixture<2,2,true,MergingType::ByEntity,NodeType::Tuple,true>;
using Basis2DQ2BlockedTupleFlatLexicographicFixture = CompositeQkFixture<2,2,true,MergingType::Lexicographic,NodeType::Tuple,false>;
using Basis2DQ2BlockedTupleBlockedLexicographicFixture = CompositeQkFixture<2,2,true,MergingType::Lexicographic,NodeType::Tuple,true>;


// list of spaces to test
using LeafBasisFixtures = ::testing::Types<
  Basis1DQ1FlatFixture,
  Basis1DQ1BlockedFixture,
  Basis2DQ1FlatFixture,
  Basis2DQ2FlatFixture,
  Basis2DQ1BlockedFixture,
  Basis2DQ2BlockedFixture,
  Basis2DQ1DgFlatFixture,
  Basis2DQ2DgFlatFixture,
  Basis2DQ3DgFlatFixture,
  Basis2DQ1DgBlockedFixture,
  Basis2DQ2DgBlockedFixture,
  Basis2DQ3DgBlockedFixture
>;

using ArrayBasisFixtures = ::testing::Types<
  Basis2DQ2FlatArrayFlatByEntityFixture,
  Basis2DQ2FlatArrayBlockedByEntityFixture,
  Basis2DQ2FlatArrayFlatLexicographicFixture,
  Basis2DQ2FlatArrayBlockedLexicographicFixture,
  Basis2DQ2BlockedArrayFlatByEntityFixture,
  Basis2DQ2BlockedArrayBlockedByEntityFixture,
  Basis2DQ2BlockedArrayFlatLexicographicFixture,
  Basis2DQ2BlockedArrayBlockedLexicographicFixture
>;

using VectorBasisFixtures = ::testing::Types<
  Basis2DQ2FlatVectorFlatByEntityFixture,
  Basis2DQ2FlatVectorBlockedByEntityFixture,
  Basis2DQ2FlatVectorFlatLexicographicFixture,
  Basis2DQ2FlatVectorBlockedLexicographicFixture,
  Basis2DQ2BlockedVectorFlatByEntityFixture,
  Basis2DQ2BlockedVectorBlockedByEntityFixture,
  Basis2DQ2BlockedVectorFlatLexicographicFixture,
  Basis2DQ2BlockedVectorBlockedLexicographicFixture
>;

using TupleBasisFixtures = ::testing::Types<
  Basis2DQ2FlatTupleFlatByEntityFixture,
  Basis2DQ2FlatTupleBlockedByEntityFixture,
  Basis2DQ2FlatTupleFlatLexicographicFixture,
  Basis2DQ2FlatTupleBlockedLexicographicFixture,
  Basis2DQ2BlockedTupleFlatByEntityFixture,
  Basis2DQ2BlockedTupleBlockedByEntityFixture,
  Basis2DQ2BlockedTupleFlatLexicographicFixture,
  Basis2DQ2BlockedTupleBlockedLexicographicFixture
>;

// Setting all the spaces into one big list blows up memory,
// so we need to split test into different exceutables
#ifdef BASIS_FIXTURES
using BasisFixtures = BASIS_FIXTURES;
#endif

#endif // DUNE_PDELAB_TEST_FIXTURE_BASIS_HH
