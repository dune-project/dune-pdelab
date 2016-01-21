// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_COMMON_TYPETRAITS_HH
#define DUNE_PDELAB_COMMON_TYPETRAITS_HH

#include <dune/common/typetraits.hh>
#include <dune/typetree/typetraits.hh>

namespace Dune {
  namespace PDELab {

    // Import AlwaysVoid from TypeTree library
    using TypeTree::AlwaysVoid;

    // forward decl of Tag defined in function.hh
    struct GridFunctionTag;
    struct PowerGridFunctionTag;
    struct CompositeGridFunctionTag;

    template<typename T, typename = void>
    struct IsGridFunction
    {
      static const bool value = false;
    };

    template<typename T>
    struct IsGridFunction<T, typename AlwaysVoid<typename T::ImplementationTag>::type >
    {
      typedef typename T::ImplementationTag A;
      static const bool value = std::is_same<A, GridFunctionTag>::value ||
        std::is_same<A, PowerGridFunctionTag>::value ||
        std::is_same<A, CompositeGridFunctionTag>::value;
    };

  } // end namespace PDELab
} // end namespace Dune

#endif // DUNE_PDELAB_COMMON_TYPETRAITS_HH
