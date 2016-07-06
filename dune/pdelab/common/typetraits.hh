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

#ifndef DOXYGEN

    namespace impl {

      template<typename T, typename = void>
      struct IsGridFunction
      {
        static const bool value = false;
      };

      template<typename T>
      struct IsGridFunction<T, typename AlwaysVoid<TypeTree::ImplementationTag<T>>::type >
      {
        using A = TypeTree::ImplementationTag<T>;
        static const bool value = std::is_same<A, GridFunctionTag>::value ||
          std::is_same<A, PowerGridFunctionTag>::value ||
          std::is_same<A, CompositeGridFunctionTag>::value;
      };

    } // namespace impl

#endif // DOXYGEN

    template<typename T>
    using IsGridFunction = std::integral_constant<bool,impl::IsGridFunction<std::decay_t<T>>::value>;

  } // end namespace PDELab
} // end namespace Dune

#endif // DUNE_PDELAB_COMMON_TYPETRAITS_HH
