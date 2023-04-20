#ifndef DUNE_PDELAB_CONCEPT_RANGE_HH
#define DUNE_PDELAB_CONCEPT_RANGE_HH

#include <dune/pdelab/concepts/indexable.hh>

#include <utility>
#include <ranges>
#include <concepts>

namespace Dune::PDELab::inline Experimental::Concept {

template<class T>
concept DenseDynamicRange = requires(T c) {
  std::begin(c); // Why not std::ranges::range ??
  std::end(c);
};

template<class T>
concept SparseDynamicRange = requires(T c) {
  requires DenseDynamicRange<T>;
  { std::ranges::begin(c).index() } -> std::integral;
};

template<class T>
concept DenseStaticRange = requires(T c) {
  requires (T::size() == 0) || requires {
    // this is the real check, but it fails on big sizes
    // []<std::size_t... I>(std::index_sequence<I...>)
    //   requires (BracketIndexable<T, std::integral_constant<std::size_t,I-1>> )
    // {}(std::make_index_sequence<T::size()>{});
    // ...instead, we just index some values and expect the rest to hold
    requires BracketIndexable<T, std::integral_constant<std::size_t,0>>;           // first
    requires BracketIndexable<T, std::integral_constant<std::size_t,T::size()/2>>; // middle
    requires BracketIndexable<T, std::integral_constant<std::size_t,T::size()-1>>; // last
  };
};

template<class T>
concept SparseStaticRange = requires(T c) {
    []<class I, I... i>(std::integer_sequence<I,i...>)
      requires (BracketIndexable<T, std::integral_constant<I,i-1>> and ...)
    {}(T::enumarate());
};

template<class T>
concept Range = DenseDynamicRange<T> || DenseDynamicRange<T> || DenseStaticRange<T> || SparseStaticRange<T>;



} // end namespace Dune::PDELab::inline Experimental::Concept


#endif // DUNE_PDELAB_CONCEPT_RANGE_HH
