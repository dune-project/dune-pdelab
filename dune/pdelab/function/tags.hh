//-*- tab-width: 4; c-basic-offset: 2; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_FUNCTION_TAGS_HH
#define DUNE_PDELAB_FUNCTION_TAGS_HH

namespace Dune {
namespace PDELab {

  struct DifferentiableFunctionBaseTag {};
  struct DifferentiableFunctionLocalViewBaseTag {};

  struct DifferentiableFunctionTag : public DifferentiableFunctionBaseTag {};
  struct DifferentiableFunctionLocalViewTag : public DifferentiableFunctionLocalViewBaseTag {};

  struct PowerDifferentiableFunctionTag : public DifferentiableFunctionBaseTag {};
  struct PowerDifferentiableFunctionLocalViewTag : public DifferentiableFunctionLocalViewBaseTag {};

  struct CompositeDifferentiableFunctionTag : public DifferentiableFunctionBaseTag {};
  struct CompositeDifferentiableFunctionLocalViewTag : public DifferentiableFunctionLocalViewBaseTag {};

}
}

#endif // DUNE_PDELAB_FUNCTION_TAGS_HH
