#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "../../grid_fixture.hh"

#include <dune/pdelab/basis/ordering/entity.hh>
#include <dune/pdelab/basis/prebasis/leaf.hh>
#include <dune/pdelab/basis/prebasis/composite.hh>

#include <dune/pdelab/basis/merging_strategy.hh>

#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/pdelab/finiteelementmap/variablemonomfem.hh>

#include <dune/pdelab/constraints/noconstraints.hh>

#include <dune/typetree/treepath.hh>

#include <dune/grid/common/scsgmapper.hh>

#include <unordered_set>

template<class MergingStragetgy, std::size_t k>
auto makeQkPreBasis(const MergingStragetgy& strategy, Dune::index_constant<k>) {
  using Q1FEM = Dune::PDELab::QkLocalFiniteElementMap<typename MergingStragetgy::EntitySet, double, double, k>;
  auto qkfem = std::make_shared<Q1FEM>(strategy.entitySet());
  return makePreBasis(strategy, qkfem);
}

template<class Q1EntityOrdering>
void testQ1EntityOrdering(Q1EntityOrdering& entity_ordering) {
  static_assert(not Q1EntityOrdering::containerBlocked());
  static_assert(Q1EntityOrdering::prioryFixedSize());
  EXPECT_TRUE(entity_ordering.fixedSize());

  static_assert(Q1EntityOrdering::maxContainerDepth() == 1);
  EXPECT_TRUE(entity_ordering.containsGeometry(Dune::GeometryTypes::vertex));
  EXPECT_FALSE(entity_ordering.containsGeometry(Dune::GeometryTypes::line));
  EXPECT_FALSE(entity_ordering.containsGeometry(Dune::GeometryTypes::quadrilateral));
  EXPECT_FALSE(entity_ordering.containsGeometry(Dune::GeometryTypes::triangle));
  EXPECT_FALSE(entity_ordering.containsCodim(0));
  EXPECT_FALSE(entity_ordering.containsCodim(1));
  EXPECT_TRUE(entity_ordering.containsCodim(2));
  EXPECT_TRUE(entity_ordering.singleCodim());
  EXPECT_FALSE(entity_ordering.disjointCodimClosure());
  EXPECT_EQ(entity_ordering.maxLocalCount(), 4);

  auto gt_index = Dune::GlobalGeometryTypeIndex::index(Dune::GeometryTypes::vertex);
  EXPECT_EQ(entity_ordering.blockCount(gt_index), 1);

  for(const auto& vertex : vertices(entity_ordering.entitySet())) {
    auto entity_index = entity_ordering.entitySet().indexSet().index(vertex);
    EXPECT_EQ(entity_ordering.blockCount(gt_index, entity_index), 1);
    auto tpath = Dune::TypeTree::treePath();
    auto mi = entity_ordering.firstContainerIndex(tpath, gt_index, entity_index);
    EXPECT_EQ(mi.size(), 1);
    EXPECT_EQ(mi[0], 0);
    EXPECT_EQ(entity_ordering.containerSize(mi, gt_index, entity_index), 0);
  }

  EXPECT_EQ(entity_ordering.blockCount(), entity_ordering.entitySet().indexSet().size(2));
  EXPECT_EQ(entity_ordering.dimension(), entity_ordering.blockCount());
}

TEST_F(StructuredGridFixture2D, TestQ1EntityOrdering) {
  using GridView = typename Dune::YaspGrid<2>::LeafGridView;
  GridView entity_set{_grid->leafGridView()};
  {
    auto flat_strategy = Dune::PDELab::Strategy::flatByEntity(entity_set);
    auto q1_space = makeQkPreBasis(flat_strategy, Dune::Indices::_1);
    using Space = decltype(q1_space);
    using EntityOrdering = Dune::PDELab::Impl::LeafEntityOrdering<Space>;
    auto entity_ordering = std::make_shared<EntityOrdering>(q1_space);
    entity_ordering->update();
    testQ1EntityOrdering(*entity_ordering);

    using VectorEntityOrdering = Dune::PDELab::Impl::VectorEntityOrdering<decltype(flat_strategy), EntityOrdering>;

    VectorEntityOrdering p_entity_ordering{std::vector{entity_ordering,entity_ordering}, flat_strategy};
    p_entity_ordering.update();
    testQ1EntityOrdering(p_entity_ordering.child(0));
    EXPECT_EQ(p_entity_ordering.childStorage(0), p_entity_ordering.childStorage(1));

    EXPECT_EQ(p_entity_ordering.degree(), 2);
    static_assert(not VectorEntityOrdering::containerBlocked());
    static_assert(VectorEntityOrdering::prioryFixedSize());
    EXPECT_TRUE(p_entity_ordering.fixedSize());

    static_assert(VectorEntityOrdering::maxContainerDepth() == 1);
    EXPECT_TRUE(p_entity_ordering.containsGeometry(Dune::GeometryTypes::vertex));
    EXPECT_FALSE(p_entity_ordering.containsGeometry(Dune::GeometryTypes::line));
    EXPECT_FALSE(p_entity_ordering.containsGeometry(Dune::GeometryTypes::quadrilateral));
    EXPECT_FALSE(p_entity_ordering.containsGeometry(Dune::GeometryTypes::triangle));
    EXPECT_FALSE(p_entity_ordering.containsCodim(0));
    EXPECT_FALSE(p_entity_ordering.containsCodim(1));
    EXPECT_TRUE(p_entity_ordering.containsCodim(2));
    EXPECT_TRUE(p_entity_ordering.singleCodim());
    EXPECT_FALSE(p_entity_ordering.disjointCodimClosure());
    EXPECT_EQ(p_entity_ordering.maxLocalCount(), 2*4);

    auto gt_index = Dune::GlobalGeometryTypeIndex::index(Dune::GeometryTypes::vertex);
    EXPECT_EQ(p_entity_ordering.blockCount(gt_index), 2);

    for(const auto& vertex : vertices(p_entity_ordering.entitySet())) {
      auto entity_index = p_entity_ordering.entitySet().indexSet().index(vertex);
      EXPECT_EQ(p_entity_ordering.blockCount(gt_index, entity_index), 2);
      auto tpath0 = Dune::TypeTree::treePath(0);
      auto mi0 = p_entity_ordering.firstContainerIndex(tpath0, gt_index, entity_index);
      EXPECT_EQ(mi0.size(), 1);
      EXPECT_EQ(mi0[0], 0);
      EXPECT_EQ(p_entity_ordering.containerSize(mi0, gt_index, entity_index), 0);

      auto tpath1 = Dune::TypeTree::treePath(1);
      auto mi1 = p_entity_ordering.firstContainerIndex(tpath1, gt_index, entity_index);
      EXPECT_EQ(mi1.size(), 1);
      EXPECT_EQ(mi1[0], 1);
      EXPECT_EQ(p_entity_ordering.containerSize(mi1, gt_index, entity_index), 0);
    }

    EXPECT_EQ(p_entity_ordering.blockCount(), 2*p_entity_ordering.entitySet().indexSet().size(2));
    EXPECT_EQ(p_entity_ordering.dimension(), p_entity_ordering.blockCount());
  }

  {
    auto blocked_strategy = Dune::PDELab::Strategy::blockedByEntity(entity_set);
    auto q1_space = makeQkPreBasis(blocked_strategy, Dune::Indices::_1);
    using Space = decltype(q1_space);
    using EntityOrdering = Dune::PDELab::Impl::LeafEntityOrdering<Space>;
    auto entity_ordering = std::make_shared<EntityOrdering>(q1_space);
    entity_ordering->update();
    testQ1EntityOrdering(*entity_ordering);

    using ArrayEntityOrdering = Dune::PDELab::Impl::ArrayEntityOrdering<decltype(blocked_strategy), EntityOrdering, 2>;

    ArrayEntityOrdering p_entity_ordering{std::array{entity_ordering,entity_ordering}, blocked_strategy};
    p_entity_ordering.update();
    testQ1EntityOrdering(p_entity_ordering.child(0));
    EXPECT_EQ(p_entity_ordering.childStorage(0), p_entity_ordering.childStorage(1));

    static_assert(ArrayEntityOrdering::degree() == 2);
    static_assert(ArrayEntityOrdering::containerBlocked());
    static_assert(ArrayEntityOrdering::prioryFixedSize());
    EXPECT_TRUE(p_entity_ordering.fixedSize());

    static_assert(ArrayEntityOrdering::maxContainerDepth() == 2);
    EXPECT_TRUE(p_entity_ordering.containsGeometry(Dune::GeometryTypes::vertex));
    EXPECT_FALSE(p_entity_ordering.containsGeometry(Dune::GeometryTypes::line));
    EXPECT_FALSE(p_entity_ordering.containsGeometry(Dune::GeometryTypes::quadrilateral));
    EXPECT_FALSE(p_entity_ordering.containsGeometry(Dune::GeometryTypes::triangle));
    EXPECT_FALSE(p_entity_ordering.containsCodim(0));
    EXPECT_FALSE(p_entity_ordering.containsCodim(1));
    EXPECT_TRUE(p_entity_ordering.containsCodim(2));
    EXPECT_TRUE(p_entity_ordering.singleCodim());
    EXPECT_FALSE(p_entity_ordering.disjointCodimClosure());
    EXPECT_EQ(p_entity_ordering.maxLocalCount(), 2*4);

    auto gt_index = Dune::GlobalGeometryTypeIndex::index(Dune::GeometryTypes::vertex);

    auto gt_count = p_entity_ordering.blockCount(gt_index);
    EXPECT_EQ(gt_count, 2);
    static_assert(std::same_as<decltype(gt_count), std::integral_constant<std::size_t,2>>);

    for(const auto& vertex : vertices(p_entity_ordering.entitySet())) {
      auto entity_index = p_entity_ordering.entitySet().indexSet().index(vertex);
      auto count = p_entity_ordering.blockCount(gt_index, entity_index);
      EXPECT_EQ(count, 2);
      static_assert(std::same_as<decltype(count), std::integral_constant<std::size_t,2>>);

      using namespace Dune::Indices;
      auto tpath0 = Dune::TypeTree::treePath(_0);
      auto mi0 = reverse(p_entity_ordering.firstContainerIndex(tpath0, gt_index, entity_index));
      EXPECT_EQ(mi0.size(), 2);
      EXPECT_EQ(mi0[0], 0);
      EXPECT_EQ(mi0[1], 0);
      EXPECT_EQ(p_entity_ordering.containerSize(mi0, gt_index, entity_index), 0);
      EXPECT_EQ(p_entity_ordering.containerSize(pop_front(mi0), gt_index, entity_index), 1);

      auto tpath1 = Dune::TypeTree::treePath(1);
      auto mi1 = reverse(p_entity_ordering.firstContainerIndex(tpath1, gt_index, entity_index));
      EXPECT_EQ(mi1.size(), 2);
      EXPECT_EQ(mi1[0], 0);
      EXPECT_EQ(mi1[1], 1);
      EXPECT_EQ(p_entity_ordering.containerSize(mi1, gt_index, entity_index), 0);
      EXPECT_EQ(p_entity_ordering.containerSize(pop_front(mi1), gt_index, entity_index), 1);

      EXPECT_EQ(p_entity_ordering.containerSize(Dune::TypeTree::treePath(), gt_index, entity_index), 2);
    }

    auto count = p_entity_ordering.blockCount();
    EXPECT_EQ(count, 2);
    static_assert(std::same_as<decltype(count), std::integral_constant<std::size_t,2>>);
    EXPECT_EQ(p_entity_ordering.dimension(), 2*p_entity_ordering.entitySet().indexSet().size(2));
  }
}

template<class Q2EntityOrdering>
void testQ2EntityOrdering(Q2EntityOrdering& entity_ordering) {
  static_assert(not Q2EntityOrdering::containerBlocked());
  static_assert(Q2EntityOrdering::prioryFixedSize());
  EXPECT_TRUE(entity_ordering.fixedSize());

  static_assert(Q2EntityOrdering::maxContainerDepth() == 1);
  EXPECT_TRUE(entity_ordering.containsGeometry(Dune::GeometryTypes::vertex));
  EXPECT_TRUE(entity_ordering.containsGeometry(Dune::GeometryTypes::line));
  EXPECT_TRUE(entity_ordering.containsGeometry(Dune::GeometryTypes::quadrilateral));
  EXPECT_FALSE(entity_ordering.containsGeometry(Dune::GeometryTypes::triangle));
  EXPECT_TRUE(entity_ordering.containsCodim(0));
  EXPECT_TRUE(entity_ordering.containsCodim(1));
  EXPECT_TRUE(entity_ordering.containsCodim(2));
  EXPECT_FALSE(entity_ordering.singleCodim());
  EXPECT_FALSE(entity_ordering.disjointCodimClosure());
  EXPECT_EQ(entity_ordering.maxLocalCount(), 9);

  auto gt_index_vertex = Dune::GlobalGeometryTypeIndex::index(Dune::GeometryTypes::vertex);
  EXPECT_EQ(entity_ordering.blockCount(gt_index_vertex), 1);

  for(const auto& vertex : vertices(entity_ordering.entitySet())) {
    auto entity_index = entity_ordering.entitySet().indexSet().index(vertex);
    EXPECT_EQ(entity_ordering.blockCount(gt_index_vertex, entity_index), 1);
    auto tpath = Dune::TypeTree::treePath();
    auto mi = entity_ordering.firstContainerIndex(tpath, gt_index_vertex, entity_index);
    EXPECT_EQ(mi.size(), 1);
    EXPECT_EQ(mi[0], 0);
    EXPECT_EQ(entity_ordering.containerSize(mi, gt_index_vertex, entity_index), 0);
  }

  auto gt_index_edge = Dune::GlobalGeometryTypeIndex::index(Dune::GeometryTypes::line);
  EXPECT_EQ(entity_ordering.blockCount(gt_index_edge), 1);

  for(const auto& edge : edges(entity_ordering.entitySet())) {
    auto entity_index = entity_ordering.entitySet().indexSet().index(edge);
    EXPECT_EQ(entity_ordering.blockCount(gt_index_edge, entity_index), 1);
    auto tpath = Dune::TypeTree::treePath();
    auto mi = entity_ordering.firstContainerIndex(tpath, gt_index_edge, entity_index);
    EXPECT_EQ(mi.size(), 1);
    EXPECT_EQ(mi[0], 0);
    EXPECT_EQ(entity_ordering.containerSize(mi, gt_index_edge, entity_index), 0);
  }

  auto gt_index_quad = Dune::GlobalGeometryTypeIndex::index(Dune::GeometryTypes::quadrilateral);
  EXPECT_EQ(entity_ordering.blockCount(gt_index_quad), 1);

  for(const auto& edge : edges(entity_ordering.entitySet())) {
    auto entity_index = entity_ordering.entitySet().indexSet().index(edge);
    EXPECT_EQ(entity_ordering.blockCount(gt_index_quad, entity_index), 1);
    auto tpath = Dune::TypeTree::treePath();
    auto mi = entity_ordering.firstContainerIndex(tpath, gt_index_edge, entity_index);
    EXPECT_EQ(mi.size(), 1);
    EXPECT_EQ(mi[0], 0);
    EXPECT_EQ(entity_ordering.containerSize(mi, gt_index_edge, entity_index), 0);
  }

  EXPECT_EQ(entity_ordering.blockCount(),
    entity_ordering.entitySet().indexSet().size(2) +
    entity_ordering.entitySet().indexSet().size(1) +
    entity_ordering.entitySet().indexSet().size(0)
  );
  EXPECT_EQ(entity_ordering.dimension(), entity_ordering.blockCount());
}


TEST_F(StructuredGridFixture2D, TestQ2EntityOrdering) {
  using GridView = typename Dune::YaspGrid<2>::LeafGridView;
  GridView entity_set{_grid->leafGridView()};

  auto flat_strategy = Dune::PDELab::Strategy::flatByEntity(entity_set);
  auto q2_pb = makeQkPreBasis(flat_strategy, Dune::Indices::_2);
  Dune::PDELab::Impl::LeafEntityOrdering entity_ordering{q2_pb};

  entity_ordering.update();
  testQ2EntityOrdering(entity_ordering);
}


template<class FixedMonomialEntityOrdering>
void testFixedMonomialEntityOrdering(FixedMonomialEntityOrdering& entity_ordering) {
  static_assert(not FixedMonomialEntityOrdering::containerBlocked());
  static_assert(not FixedMonomialEntityOrdering::prioryFixedSize());
  EXPECT_TRUE(entity_ordering.fixedSize());

  static_assert(FixedMonomialEntityOrdering::maxContainerDepth() == 1);
  EXPECT_FALSE(entity_ordering.containsGeometry(Dune::GeometryTypes::vertex));
  EXPECT_FALSE(entity_ordering.containsGeometry(Dune::GeometryTypes::line));
  EXPECT_TRUE(entity_ordering.containsGeometry(Dune::GeometryTypes::quadrilateral));
  EXPECT_FALSE(entity_ordering.containsGeometry(Dune::GeometryTypes::triangle));
  EXPECT_TRUE(entity_ordering.containsCodim(0));
  EXPECT_FALSE(entity_ordering.containsCodim(1));
  EXPECT_FALSE(entity_ordering.containsCodim(2));
  EXPECT_TRUE(entity_ordering.singleCodim());
  EXPECT_TRUE(entity_ordering.disjointCodimClosure());
  EXPECT_EQ(entity_ordering.maxLocalCount(), 6);

  auto gt_index = Dune::GlobalGeometryTypeIndex::index(Dune::GeometryTypes::quadrilateral);
  EXPECT_EQ(entity_ordering.blockCount(gt_index), 6);

  for(const auto& element : elements(entity_ordering.entitySet())) {
    auto entity_index = entity_ordering.entitySet().indexSet().index(element);
    EXPECT_EQ(entity_ordering.blockCount(gt_index, entity_index), 6);
    auto tpath = Dune::TypeTree::treePath();
    auto mi = entity_ordering.firstContainerIndex(tpath, gt_index, entity_index);
    EXPECT_EQ(mi.size(), 1);
    EXPECT_EQ(mi[0], 0);
    EXPECT_EQ(entity_ordering.containerSize(mi, gt_index, entity_index), 0);
  }

  EXPECT_EQ(entity_ordering.blockCount(),
    entity_ordering.entitySet().indexSet().size(0) * entity_ordering.blockCount(gt_index)
  );
  EXPECT_EQ(entity_ordering.dimension(), entity_ordering.blockCount());
}


TEST_F(StructuredGridFixture2D, TestFixedMonomialEntityOrdering) {
  using GridView = typename Grid::LeafGridView;
  GridView entity_set{_grid->leafGridView()};

  using CellMapper = Dune::SingleCodimSingleGeomTypeMapper<GridView, 0> ;
  CellMapper cellmapper(entity_set);
  using MonomFEM = Dune::PDELab::VariableMonomLocalFiniteElementMap<CellMapper, float, double, GridView::dimension>;
  auto monomfem = std::make_shared<MonomFEM>(cellmapper, Dune::GeometryTypes::quadrilateral, 2);

  auto strategy = Dune::PDELab::Strategy::flatByEntity(entity_set);
  auto monom_pb = makePreBasis(strategy, monomfem);

  using Space = decltype(monom_pb);
  using EntityOrdering = Dune::PDELab::Impl::LeafEntityOrdering<Space>;
  auto entity_ordering = std::make_shared<EntityOrdering>(monom_pb);

  entity_ordering->update();

  testFixedMonomialEntityOrdering(*entity_ordering);

  using ArrayEntityOrdering = Dune::PDELab::Impl::ArrayEntityOrdering<decltype(strategy), EntityOrdering, 2>;

  ArrayEntityOrdering p_entity_ordering{std::array{entity_ordering,entity_ordering}, strategy};
  p_entity_ordering.update();
  testFixedMonomialEntityOrdering(p_entity_ordering.child(0));
  EXPECT_EQ(p_entity_ordering.childStorage(0), p_entity_ordering.childStorage(1));

  static_assert(not p_entity_ordering.containerBlocked());
  static_assert(not p_entity_ordering.prioryFixedSize());
  EXPECT_TRUE(p_entity_ordering.fixedSize());

  static_assert(p_entity_ordering.maxContainerDepth() == 1);
  EXPECT_FALSE(p_entity_ordering.containsGeometry(Dune::GeometryTypes::vertex));
  EXPECT_FALSE(p_entity_ordering.containsGeometry(Dune::GeometryTypes::line));
  EXPECT_TRUE(p_entity_ordering.containsGeometry(Dune::GeometryTypes::quadrilateral));
  EXPECT_FALSE(p_entity_ordering.containsGeometry(Dune::GeometryTypes::triangle));
  EXPECT_TRUE(p_entity_ordering.containsCodim(0));
  EXPECT_FALSE(p_entity_ordering.containsCodim(1));
  EXPECT_FALSE(p_entity_ordering.containsCodim(2));
  EXPECT_TRUE(p_entity_ordering.singleCodim());
  EXPECT_TRUE(p_entity_ordering.disjointCodimClosure());
  EXPECT_EQ(p_entity_ordering.maxLocalCount(), 2*6);

  auto gt_index = Dune::GlobalGeometryTypeIndex::index(Dune::GeometryTypes::quadrilateral);
  EXPECT_EQ(p_entity_ordering.blockCount(gt_index), 2*6);

  for(const auto& element : elements(p_entity_ordering.entitySet())) {
    auto entity_index = p_entity_ordering.entitySet().indexSet().index(element);
    EXPECT_EQ(p_entity_ordering.blockCount(gt_index, entity_index), 2*6);

    auto tpath0 = Dune::TypeTree::treePath(0);
    auto mi0 = p_entity_ordering.firstContainerIndex(tpath0, gt_index, entity_index);
    EXPECT_EQ(mi0.size(), 1);
    EXPECT_EQ(mi0[0], 0);
    EXPECT_EQ(p_entity_ordering.containerSize(mi0, gt_index, entity_index), 0);

    auto tpath1 = Dune::TypeTree::treePath(1);
    auto mi1 = p_entity_ordering.firstContainerIndex(tpath1, gt_index, entity_index);
    EXPECT_EQ(mi1.size(), 1);
    EXPECT_EQ(mi1[0], 6);
    EXPECT_EQ(p_entity_ordering.containerSize(mi1, gt_index, entity_index), 0);
  }

  EXPECT_EQ(p_entity_ordering.blockCount(),
    p_entity_ordering.entitySet().indexSet().size(0) * p_entity_ordering.blockCount(gt_index)
  );
  EXPECT_EQ(p_entity_ordering.dimension(), p_entity_ordering.blockCount());
}

template<class VariableMonomialEntityOrdering>
void testVariableMonomialEntityOrdering(VariableMonomialEntityOrdering& entity_ordering) {
  static_assert(not VariableMonomialEntityOrdering::containerBlocked());
  static_assert(not VariableMonomialEntityOrdering::prioryFixedSize());
  EXPECT_FALSE(entity_ordering.fixedSize());

  static_assert(VariableMonomialEntityOrdering::maxContainerDepth() == 1);
  EXPECT_FALSE(entity_ordering.containsGeometry(Dune::GeometryTypes::vertex));
  EXPECT_FALSE(entity_ordering.containsGeometry(Dune::GeometryTypes::line));
  EXPECT_TRUE(entity_ordering.containsGeometry(Dune::GeometryTypes::quadrilateral));
  EXPECT_FALSE(entity_ordering.containsGeometry(Dune::GeometryTypes::triangle));
  EXPECT_TRUE(entity_ordering.containsCodim(0));
  EXPECT_FALSE(entity_ordering.containsCodim(1));
  EXPECT_FALSE(entity_ordering.containsCodim(2));
  EXPECT_TRUE(entity_ordering.singleCodim());
  EXPECT_TRUE(entity_ordering.disjointCodimClosure());
  EXPECT_EQ(entity_ordering.maxLocalCount(), 6);

  auto gt_index = Dune::GlobalGeometryTypeIndex::index(Dune::GeometryTypes::quadrilateral);

  std::size_t order = 0;
  for(const auto& element : elements(entity_ordering.entitySet())) {
    auto& fe = entity_ordering.space().finiteElementMap().getFEM(order++ % 3);
    auto entity_index = entity_ordering.entitySet().indexSet().index(element);
    EXPECT_EQ(entity_ordering.blockCount(gt_index, entity_index), fe.size());
    auto tpath = Dune::TypeTree::treePath();
    auto mi = entity_ordering.firstContainerIndex(tpath, gt_index, entity_index);
    EXPECT_EQ(mi.size(), 1);
    EXPECT_EQ(mi[0], 0);
    EXPECT_EQ(entity_ordering.containerSize(mi, gt_index, entity_index), 0);
  }

  EXPECT_EQ(entity_ordering.entitySet().indexSet().size(0), 4);

  EXPECT_EQ(entity_ordering.blockCount(),
    1 + 3 + 6 + 1
  );
  EXPECT_EQ(entity_ordering.dimension(), entity_ordering.blockCount());
}


TEST(TestEntityOrderings, TestVariableMonomialEntityOrdering) {
  Dune::FieldVector<double,2> L(1.0);
  std::array<int,2> N;
  std::fill(begin(N), end(N), 2); // for this test we require exactly 4 entities
  Dune::YaspGrid<2> grid{L,N};
  using GridView = typename Dune::YaspGrid<2>::LeafGridView;
  GridView entity_set{grid.leafGridView()};

  using CellMapper = Dune::SingleCodimSingleGeomTypeMapper<GridView, 0> ;
  CellMapper cellmapper(entity_set);
  using MonomFEM = Dune::PDELab::VariableMonomLocalFiniteElementMap<CellMapper, float, double, GridView::dimension>;
  auto monomfem = std::make_shared<MonomFEM>(cellmapper, Dune::GeometryTypes::quadrilateral, 2);

  std::size_t order = 0;
  for (const auto &e : elements(entity_set))
    monomfem->setOrder(e, order++ % 3);

  auto strategy = Dune::PDELab::Strategy::flatByEntity(entity_set);
  auto monom_pb = makePreBasis(strategy, monomfem);

  using Space = decltype(monom_pb);
  using EntityOrdering = Dune::PDELab::Impl::LeafEntityOrdering<Space>;
  auto entity_ordering = std::make_shared<EntityOrdering>(monom_pb);

  entity_ordering->update();

  testVariableMonomialEntityOrdering(*entity_ordering);

  using namespace Dune::Indices;
  using TupleEntityOrdering = Dune::PDELab::Impl::TupleEntityOrdering<decltype(strategy), EntityOrdering, EntityOrdering>;

  TupleEntityOrdering p_entity_ordering{std::tuple{entity_ordering,entity_ordering}, strategy};
  p_entity_ordering.update();
  testVariableMonomialEntityOrdering(p_entity_ordering.child(_0));
  testVariableMonomialEntityOrdering(p_entity_ordering.child(_1));

  static_assert(not TupleEntityOrdering::containerBlocked());
  static_assert(not TupleEntityOrdering::prioryFixedSize());
  EXPECT_FALSE(p_entity_ordering.fixedSize());

  static_assert(TupleEntityOrdering::maxContainerDepth() == 1);
  EXPECT_FALSE(p_entity_ordering.containsGeometry(Dune::GeometryTypes::vertex));
  EXPECT_FALSE(p_entity_ordering.containsGeometry(Dune::GeometryTypes::line));
  EXPECT_TRUE(p_entity_ordering.containsGeometry(Dune::GeometryTypes::quadrilateral));
  EXPECT_FALSE(p_entity_ordering.containsGeometry(Dune::GeometryTypes::triangle));
  EXPECT_TRUE(p_entity_ordering.containsCodim(0));
  EXPECT_FALSE(p_entity_ordering.containsCodim(1));
  EXPECT_FALSE(p_entity_ordering.containsCodim(2));
  EXPECT_TRUE(p_entity_ordering.singleCodim());
  EXPECT_TRUE(p_entity_ordering.disjointCodimClosure());
  EXPECT_EQ(p_entity_ordering.maxLocalCount(), 2*6);

  auto gt_index = Dune::GlobalGeometryTypeIndex::index(Dune::GeometryTypes::quadrilateral);

  order = 0;
  for(const auto& element : elements(p_entity_ordering.entitySet())) {
    auto& fe = p_entity_ordering.child(_0).space().finiteElementMap().getFEM(order++ % 3);
    auto entity_index = p_entity_ordering.entitySet().indexSet().index(element);
    EXPECT_EQ(p_entity_ordering.blockCount(gt_index, entity_index), 2*fe.size());

    EXPECT_EQ(p_entity_ordering.child(_0).blockCount(gt_index, entity_index), fe.size());
    auto tpath0 = Dune::TypeTree::treePath(_0);
    auto mi0 = p_entity_ordering.firstContainerIndex(tpath0, gt_index, entity_index);
    EXPECT_EQ(mi0.size(), 1);
    EXPECT_EQ(mi0[0], 0);
    EXPECT_EQ(p_entity_ordering.containerSize(mi0, gt_index, entity_index), 0);

    EXPECT_EQ(p_entity_ordering.child(_1).blockCount(gt_index, entity_index), fe.size());
    auto tpath1 = Dune::TypeTree::treePath(_1);
    auto mi1 = p_entity_ordering.firstContainerIndex(tpath1, gt_index, entity_index);
    EXPECT_EQ(mi1.size(), 1);
    EXPECT_EQ(mi1[0], fe.size());
    EXPECT_EQ(p_entity_ordering.containerSize(mi1, gt_index, entity_index), 0);
  }

  EXPECT_EQ(p_entity_ordering.dimension(),
    2*(1 + 3 + 6 + 1)
  );
  EXPECT_EQ(p_entity_ordering.blockCount(), p_entity_ordering.dimension());
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  Dune::MPIHelper::instance(argc, argv);
  return RUN_ALL_TESTS();
}
