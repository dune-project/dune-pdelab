#ifndef DUNE_PDELAB_LOGGING_FMT_HH
#define DUNE_PDELAB_LOGGING_FMT_HH

#include <type_traits>

#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/time.h>
#include <fmt/ostream.h>
#include <fmt/ranges.h>
#include <fmt/chrono.h>

#include <dune/pdelab/common/checks.hh>

namespace Dune::PDELab {

  /**
   * \addtogroup fmt {fmt} Extensions
   * \brief Additional Dune-specific infrastructure for `{fmt}`-based string formatting
   * \ingroup logging
   * @{
   */

  //! A compile-time string that allows constructing constexpr string_views of the captured string.
  /**
   * This class is required for the compile-time validation of _fmt format strings and produced by
   * the _fmt literal when DUNE_PDELAB_CHECK_FORMAT_STRINGS is enabled.
   */
  template<typename Char, Char... chars>
  struct format_string
    : fmt::compile_string
  {

#ifndef DOXYGEN

    // required by the fmt library
    using char_type = Char;

    // static null-terminated string constructed from the literal
    constexpr static char_type string[] = { chars..., 0 };

#endif

    //! Constructs a constexpr fmt::basic_string_view used by fmt to perform compile-time format string validation.
    constexpr operator fmt::basic_string_view<char_type>() const
    {
      return { string, sizeof...(chars) };
    }

    //! Constructs a constexpr std::string_view containing the literal.
    constexpr operator std::string_view() const
    {
      return { string, sizeof...(chars) };
    }

  };

  //! A special string_view that lets us ensure that users always use the _fmt literal for format strings.
  /**
   * This class is produced by the _fmt literal when DUNE_PDELAB_CHECK_FORMAT_STRINGS is disabled.
   */
  struct format_string_view
    : std::string_view
  {

#ifndef DOXYGEN

    constexpr format_string_view() = default;

    constexpr format_string_view(const char* s, std::size_t size)
      : std::string_view(s,size)
    {}

#endif

  };

  //! Type trait for detecting whether a given type is a format_string.
  template<typename T>
  struct is_format_string
    : public std::false_type
  {};

#ifndef DOXYGEN

  template<typename Char, Char... chars>
  struct is_format_string<format_string<Char,chars...>>
    : public std::true_type
  {};

#endif

  //! boolean value indicating whether T is a format_string.
  template<typename T>
  constexpr bool is_format_string_v = is_format_string<T>::value;


  /**
   * \brief A fixed-type holder for a factory that produces its value only when mentioned in a fmt
   * format string.
   *
   * This is a useful construct for classes that optionally offer a number of expensive-to-compute
   * arguments to users specifying a format string, see e.g. PatternFormatSink.
   *
   * The {fmt} framework will detect format string arguments of this type, invoke the factory and
   * format the result of that function call.
   */
  template<typename Factory>
  struct LazyFormatArgument
  {

    LazyFormatArgument(Factory f)
      : factory(std::move(f))
    {}

    mutable Factory factory;
  };

} // end namespace Dune::PDELab


#ifndef DOXYGEN

namespace fmt {

  template <typename Factory>
  struct formatter<Dune::PDELab::LazyFormatArgument<Factory>>
    : public formatter<std::decay_t<decltype(std::declval<Factory>()())>>
  {

    using Arg  = Dune::PDELab::LazyFormatArgument<Factory>;
    using Base = formatter<std::decay_t<decltype(std::declval<Factory>()())>>;

    // we just inherit the parsing from the original formatter

    template <typename FormatContext>
    auto format(const Arg& arg, FormatContext& ctx) {
      return Base::format(arg.factory(),ctx);
    }

  };

}

#endif // DOXYGEN

namespace Dune::Literals {

#ifdef DOXYGEN

  //! String literal operator for marking {fmt} format strings.
  /**
   * The string literator `_fmt` must be used when writing format strings for passing
   * to fmt, in particular to the logging framework. To turn a normal string literal
   * into a format string literal, just append `_fmt` to the literal:
   *
   * ~~~
   * "This is a format string with {} and {:6.3e} some data"_fmt
   * ~~~
   *
   * In order to use the custom string literal, you need to import the namespace Dune::Literals
   * with
   *
   * ~~~
   * using namespace Dune::Literals;
   * ~~~
   *
   * This will only import the literal (and possibly other Dune-related literals) without polluting
   * your scope with other Dune-related names.
   *
   * \note You do not need to import this namespace if your code lives in namespace Dune or a nested
   *       namespace, as the literal is also imported into namespace Dune.
   *
   * \sa DUNE_PDELAB_CHECK_FORMAT_STRINGS
   */
  template<typename Char, Char... chars>
  constexpr Dune::PDELab::format_string<Char,chars...> operator""_fmt()
  {
    return {};
  }


#else

#if DUNE_PDELAB_CHECK_FORMAT_STRINGS

  template<typename Char, Char... chars>
  constexpr Dune::PDELab::format_string<Char,chars...> operator""_fmt()
  {
    return {};
  }

#else

  constexpr Dune::PDELab::format_string_view operator""_fmt(const char* string, std::size_t size)
  {
    return { string, size };
  }

#endif

#endif // DOXYGEN

} // end namespace Dune::Literals


namespace Dune {

  using Dune::Literals::operator""_fmt;

}

/**
 * @} fmt
 */

#endif // DUNE_PDELAB_LOGGING_FMT_HH
