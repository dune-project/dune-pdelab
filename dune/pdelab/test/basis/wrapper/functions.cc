#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "../../fixture-grid.hh"
#include "../ordering.hh"

#include <dune/pdelab/basis/wrapper/functions.hh>
#include <dune/pdelab/basis/constraints/unconstrained.hh>
#include <dune/pdelab/basis/constraints/composite.hh>

#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>

TEST_F(StructuredGridFixture2D, TestFunctionsQ1Ordering) {

  using namespace Dune::Functions::BasisFactory;
  auto q1_basis = makeBasis(_grid->leafGridView(),lagrange<1>());
  Dune::PDELab::test_ordering( Dune::PDELab::Functions::Basis{q1_basis, Dune::PDELab::Unconstrained{}} );
}

TEST_F(StructuredGridFixture2D, TestFunctionsQ2Ordering) {
  using namespace Dune::Functions::BasisFactory;
  auto q2_basis = makeBasis(_grid->leafGridView(),lagrange<2>());
  Dune::PDELab::test_ordering( Dune::PDELab::Functions::Basis{q2_basis, Dune::PDELab::Unconstrained{}} );
}


TEST_F(StructuredGridFixture2D, TestFunctionsQ1x3Ordering) {
  using namespace Dune::Functions::BasisFactory;
  auto q1x3_basis = makeBasis(_grid->leafGridView(), power<3>(lagrange<1>(), flatInterleaved()));
  auto con = Dune::PDELab::makeCompositeConstraints(std::array<Dune::PDELab::Unconstrained,3>());
  Dune::PDELab::test_ordering( Dune::PDELab::Functions::Basis{q1x3_basis, con} );
}

TEST_F(StructuredGridFixture2D, TestFunctionsQ1x3BlockedOrdering) {
  using namespace Dune::Functions::BasisFactory;
  auto q1x3_basis = makeBasis(_grid->leafGridView(), power<3>(lagrange<1>(), blockedInterleaved()));
  auto con = Dune::PDELab::makeCompositeConstraints(std::array<Dune::PDELab::Unconstrained,3>());
  Dune::PDELab::test_ordering( Dune::PDELab::Functions::Basis{q1x3_basis, con} );
}


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  Dune::MPIHelper::instance(argc, argv);
  return RUN_ALL_TESTS();
}
