#ifndef DUNE_PDELAB_BASIS_CONSTRAINTS_UNCONSTRAINED_HH
#define DUNE_PDELAB_BASIS_CONSTRAINTS_UNCONSTRAINED_HH

#include <dune/pdelab/concepts/multiindex.hh>

#include <dune/pdelab/basis/constraints/container_empty.hh>

#include <dune/typetree/leafnode.hh>

#include <dune/grid/concepts/gridview.hh>

namespace Dune::PDELab::inline Experimental {

  struct Unconstrained : public TypeTree::LeafNode {

    // template<Concept::MultiIndex MultiIndex, Dune::Concept::GridView EntitySet> // gcc11.2 error: exceeds the maximum constraint complexity -_-
    template<class MultiIndex, class EntitySet>
    using Container = EmptyConstraintsContainer<MultiIndex, EntitySet>;

    Unconstrained() = default;

    static constexpr bool doConstrainBoundary() {return false;}
    static constexpr bool doConstrainSkeleton() {return false;}
    static constexpr bool doConstrainVolume()   {return false;}

    void constrainVolume(const Concept::LocalBasisTree auto& lbasis, auto& container) const
    {}

    void constrainSkeleton(
      const Dune::Concept::Intersection auto& intersection,
      const Concept::LocalBasisTree auto& lbasis_in,
      const Concept::LocalBasisTree auto& lbasis_out,
      auto& container) const
    {}

    void constrainBoundary(
      const Dune::Concept::Intersection auto& intersection,
      const Concept::LocalBasisTree auto& lbasis_in,
      auto& container) const
    {}
  };


} // namespace Dune::PDELab::inline Experimental

#endif // DUNE_PDELAB_BASIS_CONSTRAINTS_UNCONSTRAINED_HH
