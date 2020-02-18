// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/uggrid.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/pdelab.hh>

#include "gridexamples.hh"

namespace Dune::PDELab {
template<typename S, typename B>
struct FunctionsBasisInfo
{
  struct Traits {
    using Basis = S;
    using Backend = B;
  };
};
}

template<typename Grid>
void test(const std::unique_ptr<Grid>& grid)
{
  using namespace Dune::Functions::BasisFactory;

  const auto gridView = grid->leafGridView();
  auto basis = makeBasis(gridView, lagrange<3>());

  using Backend = Dune::PDELab::ISTL::VectorBackend<>;
  using BasisInfo = Dune::PDELab::FunctionsBasisInfo<decltype(basis),Backend>;

  using Vec = Dune::PDELab::Backend::Vector<BasisInfo,double>;
  Vec x(basis);

}

int main(int argc, char** argv)
{
  try{

    Dune::MPIHelper::instance(argc,argv);

#if HAVE_UG
    test(UnitTriangleMaker<Dune::UGGrid<2>>::create());
    // test(TriangulatedUnitSquareMaker<Dune::UGGrid<2>>::create());
#endif // HAVE_UG

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
