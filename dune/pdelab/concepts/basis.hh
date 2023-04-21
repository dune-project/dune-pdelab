#ifndef DUNE_PDELAB_CONCEPTS_BASIS_HH
#define DUNE_PDELAB_CONCEPTS_BASIS_HH

#include <dune/pdelab/concepts/multiindex.hh>
#include <dune/pdelab/concepts/local_basis.hh>
#include <dune/pdelab/concepts/local_index_set.hh>
#include <dune/pdelab/concepts/local_constraints.hh>
#include <dune/pdelab/concepts/size_provider.hh>

#include <dune/grid/concepts/gridview.hh>

#include <concepts>

namespace Dune::PDELab::inline Experimental::Concept {

  /**
   * @brief Concept for discrete function spaces
   *
   * @tparam B
   */
  template<class B>
  concept Basis = requires(const B basis, std::size_t dim, std::size_t codim)
  {
    requires std::integral<typename B::size_type>;
    requires Dune::Concept::GridView<typename B::EntitySet>;
    requires SizeProvider<B>;
    requires LocalBasis<typename B::LocalView>;
    requires LocalIndexSet<typename B::LocalIndexSet>; // TODO loop over codim to bind each entity
    requires LocalConstraints<typename B::LocalConstraints>;
    requires std::equality_comparable<B>;
    { basis.entitySet()           } -> std::convertible_to<typename B::EntitySet>;
    { basis.localView()           } -> std::convertible_to<typename B::LocalView>;
    { basis.localIndexSet()       } -> std::convertible_to<typename B::LocalIndexSet>;
    { basis.localConstraints()    } -> std::convertible_to<typename B::LocalConstraints>;
    { basis.dimension()           } -> std::convertible_to<typename B::size_type>;
    { basis.size()                } -> std::convertible_to<typename B::size_type>;
    { basis.degree()              } -> std::convertible_to<typename B::size_type>;
    { basis.fixedSize(dim, codim) } -> std::convertible_to<bool>;
    { basis.contains(dim, codim)  } -> std::convertible_to<bool>;
    requires requires(B mutable_basis, typename B::EntitySet entity_set) {
      mutable_basis.update(entity_set);
    };
  };

} // namespace Dune::PDELab::inline Experimental::Concept

#endif // DUNE_PDELAB_CONCEPTS_BASIS_HH
