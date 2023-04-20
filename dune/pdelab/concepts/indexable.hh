#ifndef DUNE_INDEXABLE_CONCEPT_INDEXABLE_HH
#define DUNE_INDEXABLE_CONCEPT_INDEXABLE_HH

#include <utility>
#include <ranges>
#include <concepts>

namespace Dune::PDELab::inline Experimental::Concept {

template<class C, class Index>
concept BracketIndexable = requires(C&& c, Index index)
{
  c[index];
};

} // end namespace Dune::PDELab::inline Experimental::Concept


#endif // DUNE_INDEXABLE_CONCEPT_INDEXABLE_HH
