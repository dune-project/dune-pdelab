#ifndef DUNE_PDELAB_CONCEPT_SIZE_PROVIDER_HH
#define DUNE_PDELAB_CONCEPT_SIZE_PROVIDER_HH

#include <dune/pdelab/concepts/multiindex.hh>

#include <concepts>

namespace Dune::PDELab::inline Experimental::Concept {

template<class SP>
concept SizeProvider = requires(const SP& size_provider, typename SP::SizePrefix prefix)
{
  { size_provider.size(prefix) } -> std::convertible_to<typename SP::size_type>;
  requires std::integral<typename SP::size_type>;
  requires MultiIndex<typename SP::SizePrefix>;
  requires std::copy_constructible<SP>;
};

} // end namespace Dune::PDELab::inline Experimental::Concept

#endif // DUNE_PDELAB_CONCEPT_SIZE_PROVIDER_HH
