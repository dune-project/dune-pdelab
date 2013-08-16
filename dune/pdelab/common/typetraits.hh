// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_COMMON_TYPETRAITS_HH
#define DUNE_PDELAB_COMMON_TYPETRAITS_HH

#include <dune/common/typetraits.hh>

#include "function.hh"

namespace Dune {

  // Provide some more C++11 TMP helpers.
  // These should be upstreamed to dune-common ASAP.

#if defined HAVE_TYPE_TRAITS

  // Tests whether the first template argument is a base class of the second one.
  using std::is_base_of;
  // C++11 equivalent of Dune::SelectType.
  using std::conditional;

#elif defined HAVE_TR1_TYPE_TRAITS

  // This is already in TR1
  using std::tr1::is_base_of;

  // We have to reimplement std::conditional, but that's trivial...
  template<bool B, typename T, typename F>
  struct conditional
  {
    typedef T type;
  };

  template<typename T, typename F>
  struct conditional<false,T,F>
  {
    typedef F type;
  };

#else

  // Every compiler we currently support should have TR1.
#error Your compiler supports neither C++11 nor TR1!
#error PDELab requires at least TR1 support in your compiler to work, bailing out...

#endif

  namespace PDELab {

    template<typename>
    struct AlwaysVoid
    {
      typedef void type;
    };

    // // forward decl of Tag defined in function.hh
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


    // Support for lazy evaluation of meta functions. This is required when doing
    // nested tag dispatch without C++11-style typedefs (based on using syntax).
    // The standard struct-based meta functions cause premature evaluation in a
    // context that is not SFINAE-compatible. We thus have to return the meta function
    // without evaluating it, placing that burden on the caller. On the other hand,
    // the lookup will often directly the target type, so here is some helper code
    // to automatically do the additional evaluation if necessary.
    // Too bad that the new syntax is GCC 4.6+...


    //! Marker tag declaring a meta function.
    /**
     * Just inherit from this type to cause lazy evaluation
     */
    struct meta_function {};

    //! Helper meta function to delay evaluation of F.
    template<typename F>
    struct lazy_evaluate
    {
      typedef typename F::type type;
    };

    //! Identity function.
    template<typename F>
    struct lazy_identity
    {
      typedef F type;
    };

    //! Meta function that evaluates its argument iff it inherits from meta_function.
    template<typename F>
    struct evaluate_if_meta_function
    {
      typedef typename conditional<
        is_base_of<meta_function,F>::value,
        lazy_evaluate<F>,
        lazy_identity<F>
        >::type::type type;
    };

  } // end namespace PDELab
} // end namespace Dune

#endif // DUNE_PDELAB_COMMON_TYPETRAITS_HH
