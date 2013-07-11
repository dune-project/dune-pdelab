// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_COMMON_TYPETRAITS_HH
#define DUNE_PDELAB_COMMON_TYPETRAITS_HH

#include <dune/common/typetraits.hh>

#include "function.hh"

namespace Dune {
  namespace PDELab {

    template<typename>
    struct AlwaysVoid
    {
      typedef void type;
    };

    // // forward decl of Tag defined in function.hh
    // struct GridFunctionTag;

    template<typename T, typename = void>
    struct IsGridFunction
    {
      static const bool value = false;
    };

    template<typename T>
    struct IsGridFunction<T, typename AlwaysVoid<typename T::ImplementationTag>::type >
    {
      typedef typename T::ImplementationTag A;
      static const bool value = is_same<A, GridFunctionTag>::value ||
        is_same<A, PowerGridFunctionTag>::value ||
        is_same<A, CompositeGridFunctionTag>::value;
    };


#ifndef DOXYGEN

// Make sure we have decltype or a compatible fall back

#if HAVE_STD_DECLTYPE
#define DUNE_DECLTYPE decltype
#elif HAVE_GCC___TYPEOF__
#define DUNE_DECLTYPE __typeof__
#else
#error The TypeTree library (and by extension PDELab) require support for
#error C++11 decltype or a compatible fallback in your compiler.
#error Neither of those was found, aborting!!!!
#endif

#endif // DOXYGEN


    //! Helper function for generating a pointer to a value of type T in an unevaluated operand setting.
    template<typename T>
    T* declptr();

  } // end namespace PDELab
} // end namespace Dune

#endif // DUNE_PDELAB_COMMON_TYPETRAITS_HH
