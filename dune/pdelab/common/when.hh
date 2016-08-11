#ifndef DUNE_PDELAB_COMMON_WHEN_HH
#define DUNE_PDELAB_COMMON_WHEN_HH

#include <type_traits>

namespace Dune {
  namespace PDELab {

#ifndef DOXYGEN

    namespace impl {

      template<bool enabled>
      struct otherwise_holder;

      template<>
      struct otherwise_holder<true>
      {
        template<typename F>
        static void otherwise(F f)
        {
          f(0);
        }
      };

      template<>
      struct otherwise_holder<false>
      {
        template<typename F>
        static void otherwise(F f)
        {}
      };


      template<typename F>
      otherwise_holder<false> _when(std::true_type, F f)
      {
        f(0);
        return {};
      }

      template<typename F>
      otherwise_holder<true> _when(std::false_type, F f)
      {
        return {};
      }

    }

#endif // DOXYGEN

    template<bool condition, typename F>
    auto when(F f)
    {
      return impl::_when(std::integral_constant<bool,condition>{},f);
    }


  } // namespace PDELab
} //namespace Dune

#endif // DUNE_PDELAB_COMMON_WHEN_HH
