#ifndef DUNE_ASSEMBLER_SPACE_CONSTRAINTS_DIRICHLET_HH
#define DUNE_ASSEMBLER_SPACE_CONSTRAINTS_DIRICHLET_HH

#include <dune/assembler/space/constraints/container_affine.hh>

#include <dune/functions/gridfunctions/gridfunction.hh>

#include <dune/typetree/leafnode.hh>

#include <type_traits>

namespace Dune::Assembler {

  // TODO check grid function concept
  template<class GridFunction>
  class DirichletConstraints : public TypeTree::LeafNode {

    // range may be a double (dimensionless) or the actual value type of the translation (unit)
    using Range = typename Functions::SignatureTraits<typename GridFunction::Signature>::Range;

    // allow value to be wrapped in std::optional
    template<class T> struct remove_optional                   : std::type_identity<T> {};
    template<class T> struct remove_optional<std::optional<T>> : std::type_identity<T> {};
    using Value = remove_optional<Range>::type;

    using LocalFunction = typename GridFunction::LocalFunction;

  public:

    static constexpr bool doConstrainBoundary() {return true;}
    static constexpr bool doConstrainSkeleton() {return false;}
    static constexpr bool doConstrainVolume()   {return false;}


    template<Concept::MultiIndex MultiIndex, Dune::Concept::GridView EntitySet>
    using Container = AffineConstraintsContainer<Value, MultiIndex, EntitySet>;

    DirichletConstraints(const GridFunction& grid_function)
      : _grid_function{grid_function}
      , _local_function{localFunction(_grid_function)}
    {}

    void constrainBoundary(
      const Dune::Concept::Intersection auto& intersection,
      const Concept::LocalSpaceLeaf auto& lspace_in,
      auto& container)
    {
      if (lspace_in.size() == 0) return;

      const auto& entity_in = lspace_in.element();
      // find dof indices that belong to the intersection
      const auto face = intersection.indexInInside();
      const auto& refelem = referenceElement(entity_in.geometry());
      _local_function.bind(entity_in);
      const auto& lkeys = lspace_in.finiteElement().localCoefficients();

      for (std::size_t dof=0; dof != lspace_in.size(); ++dof) {
        // the codim to which this dof is attached to
        unsigned int codim = lkeys.localKey(dof).codim();
        if (codim==0) continue;

        // find the reference sub_entity index for this degree of freedom
        int sub_entity = lkeys.localKey(dof).subEntity();
        for (int j=0; j != refelem.size(face,1,codim); ++j) {
          if (sub_entity == refelem.subEntity(face,1,j,codim)) {
            auto coord = refelem.position(sub_entity, codim);
            // evaluate local function at the center of the sub-entity
            std::optional<Value> dirichlet = _local_function(coord);
            if (dirichlet)
              container.addTranslationConstraint(lspace_in.index(dof), *dirichlet);
          }
        }
      }

      _local_function.unbind();
    }

    void constrainVolume(const Concept::LocalSpaceLeaf auto& lspace, auto& container) const
    {}

    void constrainSkeleton(
      const Dune::Concept::Intersection auto& intersection,
      const Concept::LocalSpaceLeaf auto& lspace_in,
      const Concept::LocalSpaceLeaf auto& lspace_out,
      auto& container) const
    {}

    GridFunction _grid_function;
    LocalFunction _local_function;
  };

} // namespace Dune::Assembler

#endif // DUNE_ASSEMBLER_SPACE_CONSTRAINTS_DIRICHLET_HH
