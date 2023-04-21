#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/pdelab/basis/constraints/container_affine.hh>
#include <dune/pdelab/basis/constraints/container.hh>
#include <dune/pdelab/basis/constraints/dirichlet.hh>
#include <dune/pdelab/common/multiindex.hh>

#include <dune/functions/gridfunctions/analyticgridviewfunction.hh>

#include <dune/grid/test/basicunitcube.hh>

#include <gtest/gtest.h>

#include <random>

#include "../fixture.hh"

TEST(TestConstraintsContainer, FixedConstraints) {

  Dune::PDELab::Impl::AffineConstraintsContainerBase<double, Dune::TypeTree::HybridTreePath<int>> constraints;

  using Dune::TypeTree::treePath;

  //     | 1    0     0    0  0 |      |  0 |
  //     | 0    0     0    0  0 |      | 10 |
  // C = | 0.5  0.5   0    0  0 |  b = |  0 |
  //     | 0.5  0.25  0.5  0  0 |      |  5 |
  //     | 0    0     0    0  1 |      |  0 |

  constraints.addLinearConstraint(treePath(2), std::array{std::pair{treePath(0), 0.5},
                                                          std::pair{treePath(1), 0.5}});
  constraints.addAffineConstraint(treePath(3), 5., std::array{std::pair{treePath(0), 0.5},
                                                              std::pair{treePath(1), 0.25},
                                                              std::pair{treePath(2), 0.5}});
  constraints.addTranslationConstraint(treePath(1), 10.);

  constraints.globalCompress();

  std::vector<double> x_init{1., 2., 3., 4., 5.};

  EXPECT_FALSE( constraints.isConstrained( treePath(0) ) );
  EXPECT_TRUE(  constraints.isConstrained( treePath(1) ) );
  EXPECT_TRUE(  constraints.isConstrained( treePath(2) ) );
  EXPECT_TRUE(  constraints.isConstrained( treePath(3) ) );
  EXPECT_FALSE( constraints.isConstrained( treePath(4) ) );

  auto x = x_init;

  // x = Cx
  constraints.applyLinearConstraints(x = x_init);

  std::vector x_linear{
    x_init[0],
    0.,
    x_init[0]*0.5+x_init[1]*0.5,
    x_init[0]*0.5+x_init[1]*0.25+x_init[2]*0.5,
    x_init[4]
  };
  for (std::size_t i = 0; i != x.size(); ++i)
    EXPECT_FLOAT_EQ(x[i], x_linear[i]) << "  i:=" << i;

  // x = C^Tx
  constraints.applyLinearTransposedConstraints(x = x_init);

  std::vector x_linear_t{
    x_init[0]+x_init[2]*0.5+x_init[3]*0.5,
    x_init[2]*0.5+x_init[3]*0.25,
    x_init[3]*0.5,
    0.,
    x_init[4]
  };
  for (std::size_t i = 0; i != x.size(); ++i)
    EXPECT_FLOAT_EQ(x[i], x_linear_t[i]) << "  i:=" << i;

  // x = Cx+b
  constraints.applyAffineConstraints(x = x_init);

  std::vector x_translation{0., 10., 0., 5., 0., 0.};
for (std::size_t i = 0; i != x.size(); ++i)
    EXPECT_FLOAT_EQ(x[i], constraints.isConstrained( treePath(int(i)) ) ? x_linear[i] + x_translation[i] : x_init[i]) << "  i:=" << i;
}

TEST(TestConstraintsContainer, RandomConstraints) {

  using MultiIndex = Dune::TypeTree::HybridTreePath<std::size_t>;
  Dune::PDELab::Impl::AffineConstraintsContainerBase<double, MultiIndex> constraints;

  using Dune::TypeTree::treePath;
  std::size_t size = 10;

  std::default_random_engine re;
  re.seed(0);
  std::uniform_int_distribution<std::size_t> sunif(0, size-1);

  // generate a random set of constrained DOFs
  std::set<std::size_t> constrained_rows;
  for (std::size_t i = 0; i != size/2; ++i)
    constrained_rows.insert(sunif(re));

  std::uniform_int_distribution<std::size_t> bunif(0, 1);

  std::uniform_real_distribution<double> dunif(-100., 100.);

  std::vector<std::vector<double>> C(size);
  std::vector<double> b(size, 0.);
  std::vector<double> x_init(size, 0.);

  // start with unconstrained (C is the identity)
  for (std::size_t row = 0; row != size; ++row) {
    C[row].assign(size, 0.);
    C[row][row] = 1.;
  }

  // add constraints randomly and save them in C and b
  using Row = std::vector<std::pair<MultiIndex, double>>;
  std::vector<Row> C_rows(size);
  for (auto row : constrained_rows) {
    x_init[row] = dunif(re);
    C[row][row] = 0.;
    if (bunif(re))
      b[row] = dunif(re);
    if (bunif(re)) {
      std::set<std::size_t> constrained_cols;
      for (std::size_t i = 0; i != size; ++i)
        constrained_cols.insert(sunif(re));
      constrained_cols.erase(row);
      for (auto col : constrained_cols)
        C_rows[row].emplace_back(treePath(col), C[row][col] = dunif(re));
    }
    constraints.addAffineConstraint(treePath(row), b[row], C_rows[row]);
  }

  constraints.globalCompress();

  for(std::size_t row = 0; row != size; ++row)
    EXPECT_EQ(constraints.isConstrained( treePath(row) ), constrained_rows.contains(row) ) << "row:=" << row;

  auto x = x_init;

  // x = Cx
  constraints.applyLinearConstraints(x = x_init);

  auto x_linear = x_init;
  for(std::size_t row = 0; row != size; ++row) {
    x_linear[row] = 0.;
    for(std::size_t col = 0; col != size; ++col)
      x_linear[row] += C[row][col] * x_init[col];
  }

  for (std::size_t i = 0; i != x.size(); ++i)
    EXPECT_FLOAT_EQ(x[i], x_linear[i]) << "  i:=" << i << " constrained:=" << constrained_rows.contains(i);

  // x = C^Tx
  constraints.applyLinearTransposedConstraints(x = x_init);

  auto x_linear_t = x_init;
  for(std::size_t row = 0; row != size; ++row) {
    x_linear_t[row] = 0.;
    for(std::size_t col = 0; col != size; ++col)
      x_linear_t[row] += C[col][row] * x_init[col];
  }

  for (std::size_t i = 0; i != x.size(); ++i)
    EXPECT_FLOAT_EQ(x[i], x_linear_t[i]) << "  i:=" << i << " constrained:=" << constrained_rows.contains(i);

  // x = Cx+b
  constraints.applyAffineConstraints(x = x_init);

  auto x_affine = b;
  for(std::size_t row = 0; row != size; ++row)
    for(std::size_t col = 0; col != size; ++col)
      x_affine[row] += C[row][col] * x_init[col];

  for (std::size_t i = 0; i != x.size(); ++i) {
    EXPECT_FLOAT_EQ(x[i], x_affine[i]) << "  i:=" << i << " constrained:=" << constrained_rows.contains(i);
    if (not constrained_rows.contains(i)) {
      EXPECT_FLOAT_EQ(x[i], x_init[i]) << "  i:=" << i << " constrained:=" << constrained_rows.contains(i);
    }
  }
}

TEST_F(Basis2DQ2BlockedFixture, DirichletConstraints) {

  using EntitySet = typename Basis2DQ2BlockedFixture::Grid::LeafGridView;

  auto basis = makeFixtureBasis(_grid->leafGridView());
  using Basis = decltype(basis);
  using MultiIndex = typename Basis::LocalView::Tree::MultiIndex;
  using ConstraintsContainerNode = Dune::PDELab::AffineConstraintsContainer<double, MultiIndex, EntitySet>;
  auto container_node = std::make_shared<ConstraintsContainerNode>(basis.entitySet());

  Dune::PDELab::ConstraintsContainer container{container_node};

  auto dirichlet_function = Dune::Functions::makeAnalyticGridViewFunction(
    [](auto coordinate) -> std::optional<double> {
      // these values are always on the boundary of this grid (we rely on this for this test)
      if (Dune::FloatCmp::eq<double>(coordinate[0], 0.))
        return coordinate.two_norm2();
      else
        return std::nullopt;
    }, basis.entitySet());


  auto constraints_ops = Dune::TypeTree::makeTreeContainer(container.tree(), [&](auto& pre_basis_node){
    return Dune::PDELab::DirichletConstraints{dirichlet_function};
  });

  container.assembleConstraints(basis, constraints_ops);

  auto lbasis = basis.localView();
  auto lconstraints = container.localView(lbasis.tree());

  for (const auto& entity : elements(basis.entitySet())) {
    lbasis.bind(entity);
    lconstraints.bind(entity);

    const auto& geo = entity.geometry();
    const auto& refelem = referenceElement(geo);
    const auto& lkeys = lbasis.tree().finiteElement().localCoefficients();

    for (std::size_t dof = 0; dof != lbasis.size(); ++dof) {

      // the codim to which this dof is attached to
      unsigned int codim = lkeys.localKey(dof).codim();
      if (codim==0) continue;

      // find the reference sub_entity index for this degree of freedom
      int sub_entity = lkeys.localKey(dof).subEntity();

      for (int j=0; j != refelem.size(0,0,codim); ++j) {
        if (sub_entity == refelem.subEntity(0,0,j,codim)) {
          auto coord = geo.global(refelem.position(sub_entity, codim));
          // evaluate local function at the center of the sub-entity
          auto val = dirichlet_function(coord);
          EXPECT_EQ(lconstraints.tree().isConstrained(dof), val.has_value())
            << " dof:=" << lbasis.tree().index(dof)
            << " coord:=" << coord
            << " does not match local and global constraints";
          if (val.has_value()) {
            EXPECT_EQ(lconstraints.tree().translationValue(dof), val.value())
              << " dof:=" << lbasis.tree().index(dof)
              << " coord:=" << coord
              << " does not match local and global constraints";
            EXPECT_EQ(lconstraints.tree().linearCoefficients(dof).size(), 0)
              << " dof:=" << lbasis.tree().index(dof)
              << " coord:=" << coord
              << " does not match local and global constraints";
          }
        }
      }
    }

    lbasis.unbind();
    lconstraints.unbind();
  }
}

#if HAVE_DUNE_UGGRID

TEST(Basis2DQ1Flat, HangingNodeConstraints) {
  using Grid = Dune::UGGrid<2>;
  Dune::GridFactory<Grid> gf;

  BasicUnitCube<2>::insertVertices(gf);
  BasicUnitCube<2>::insertCubes(gf);

  auto _grid = std::unique_ptr<Grid>(gf.createGrid());
  _grid->setClosureType(Grid::ClosureType::NONE);

  auto refine_first_entity = [&]{
    auto entity = *(_grid->leafGridView().template begin<0>());
    _grid->mark(1,entity);
    _grid->preAdapt();
    _grid->adapt();
    _grid->postAdapt();
  };

  refine_first_entity();
  refine_first_entity();

  /*
    +--+--+-----+
    |  |  |     |
    +--+--o-----+
    |  |  |     |
    +--o--+-----+
    |     |     |
    |     |     |
    |     |     |
    +-----+-----+
  */

  using EntitySet = typename Grid::LeafGridView;
  EntitySet entity_set = _grid->leafGridView();

  // Dune::VTKWriter writer(entity_set, Dune::VTK::nonconforming);;
  // writer.write("nonconforming_grid");

  auto basis = makeFixtureBasis(makeQkBasis<1,false>(entity_set), entity_set);

  using Basis = decltype(basis);
  using MultiIndex = typename Basis::LocalView::Tree::MultiIndex;
  using ConstraintsContainerNode = Dune::PDELab::AffineConstraintsContainer<double, MultiIndex, EntitySet>;
  auto container_node = std::make_shared<ConstraintsContainerNode>(basis.entitySet());

  Dune::PDELab::ConstraintsContainer container{container_node};

  container.clear();
  // DoF 10 is the dof  non-conforming entity
  auto c0 = std::pair{MultiIndex(6), 0.45};
  auto c1 = std::pair{MultiIndex(7), 0.55};
  auto c2 = std::pair{MultiIndex(8), 0.45};
  container.tree().addLinearConstraint(MultiIndex(10), std::array{c1, c0});
  container.tree().addLinearConstraint(MultiIndex(11), std::array{c1, c2});
  container.compress(basis);
  // container.tree().debug();

  auto lbasis = basis.localView();
  auto lconstraints = container.localView(lbasis.tree());

  std::size_t count = 0;
  for (const auto& entity : elements(basis.entitySet())) {
    lbasis.bind(entity);
    lconstraints.bind(entity);

    // const auto& lkeys = lbasis.tree().finiteElement().localCoefficients();
    // std::cout << "center:=" << entity.geometry().center() << std::endl;

    for (std::size_t dof = 0; dof != lbasis.tree().size(); ++dof) {
      // auto pos = entity.template subEntity<2>(lkeys.localKey(dof).subEntity()).geometry().center();
      // std::cout << "dof:=" << dof << " index:=" << lbasis.tree().index(dof) << " pos:="  << pos << std::endl;
      if (lbasis.tree().index(dof) == MultiIndex(10)) {
        ++count;
        EXPECT_TRUE(lconstraints.tree().isConstrained(dof));
        EXPECT_DOUBLE_EQ(lconstraints.tree().translationValue(dof), 0.);
        auto linear_coefficients = lconstraints.tree().linearCoefficients(dof);
        // linear coefficeints are sorted lexicograpically by the multi-index
        EXPECT_EQ(linear_coefficients[0], c0);
        EXPECT_EQ(linear_coefficients[1], c1);
        if (linear_coefficients.size() == 3) {
          EXPECT_EQ(linear_coefficients[2], c2);
        }
      } else if (lbasis.tree().index(dof) == MultiIndex(11)) {
        ++count;
        EXPECT_TRUE(lconstraints.tree().isConstrained(dof));
        EXPECT_DOUBLE_EQ(lconstraints.tree().translationValue(dof), 0.);
        auto linear_coefficients = lconstraints.tree().linearCoefficients(dof);
        // linear coefficeints are sorted lexicograpically by the multi-index
        if (linear_coefficients.size() == 2) {
          EXPECT_EQ(linear_coefficients[0], c1);
          EXPECT_EQ(linear_coefficients[1], c2);
        } else if (linear_coefficients.size() == 3) {
          EXPECT_EQ(linear_coefficients[0], c0);
          EXPECT_EQ(linear_coefficients[1], c1);
          EXPECT_EQ(linear_coefficients[2], c2);
        }
      } else {
        EXPECT_FALSE(lconstraints.tree().isConstrained(dof));
      }
    }

    lbasis.unbind();
    lconstraints.unbind();
  }
  EXPECT_EQ(count, 4);
}
#endif // HAVE_DUNE_UGGRID

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  Dune::MPIHelper::instance(argc, argv);
  return RUN_ALL_TESTS();
}
