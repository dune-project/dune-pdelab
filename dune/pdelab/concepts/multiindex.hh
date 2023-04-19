#ifndef DUNE_PDELAB_CONCEPTS_MULTIINDEX_HH
#define DUNE_PDELAB_CONCEPTS_MULTIINDEX_HH

#include <concepts>
#include <type_traits>

namespace Dune::PDELab::inline Experimental::Concept {

  namespace Impl {

    template<class MI>
    concept MultiIndexEssentials = requires(MI mi) {
      { mi.size() }        -> std::convertible_to<std::size_t>;
      { MI::max_size() }   -> std::convertible_to<std::size_t>;
      requires std::copy_constructible<MI>;
      requires std::is_copy_assignable_v<MI>;
      requires std::move_constructible<MI>;
      requires std::is_move_assignable_v<MI>;
    };

    template<class MI>
    concept NonEmptyMultiIndex = requires(MI mi, std::integral_constant<std::size_t,0> _0)
    {
      requires MultiIndexEssentials<MI>;
      { back(mi) }                        -> std::convertible_to<std::size_t>;
      { front(mi) }                       -> std::convertible_to<std::size_t>;
      { pop_back(mi) }                    -> MultiIndexEssentials;
      { pop_front(mi) }                   -> MultiIndexEssentials;
      { back(accumulate_back(mi, _0)) }   -> std::convertible_to<std::size_t>;
      { back(accumulate_back(mi, 0)) }    -> std::convertible_to<std::size_t>;
      { front(accumulate_front(mi, _0)) } -> std::convertible_to<std::size_t>;
      { front(accumulate_front(mi, 0)) }  -> std::convertible_to<std::size_t>;
    };

  }

  /**
   * @brief Concept for container of multiple indices
   *
   * \f[ \alpha = (\alpha_1, \alpha_2,\ldots,\alpha_n)  \f]
   *
   */
  template<class MI>
  concept MultiIndex = requires(MI mi, std::integral_constant<std::size_t,0> _0)
  {
    requires Impl::MultiIndexEssentials<MI>;
    { push_front(mi, _0) }          -> Impl::NonEmptyMultiIndex;
    { push_front(mi, 0) }           -> Impl::NonEmptyMultiIndex;
    { push_back(mi, _0) }           -> Impl::NonEmptyMultiIndex;
    { push_back(mi, 0) }            -> Impl::NonEmptyMultiIndex;
    { join(push_back(mi, 0), mi) }  -> Impl::NonEmptyMultiIndex;
    { join(push_back(mi, _0), mi) } -> Impl::NonEmptyMultiIndex;
    { reverse(push_back(mi, _0)) }  -> Impl::NonEmptyMultiIndex;
    { reverse(push_back(mi, 0)) }   -> Impl::NonEmptyMultiIndex;
  };

  template<class MI>
  concept FixedSizeMultiIndex = requires(MI mi)
  {
    requires MultiIndex<MI>;
    { std::integral_constant<std::size_t, MI::size()>{} } -> std::convertible_to<std::size_t>;
  };


} // namespace Dune::PDELab::inline Experimental::Concept

#endif // DUNE_PDELAB_CONCEPTS_MULTIINDEX_HH
