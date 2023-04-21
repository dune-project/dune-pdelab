#ifndef DUNE_PDELAB_BASIS_CONSTRAINTS_CONTAINER_EMPTY_HH
#define DUNE_PDELAB_BASIS_CONSTRAINTS_CONTAINER_EMPTY_HH

#include <dune/pdelab/concepts/basis.hh>
#include <dune/pdelab/concepts/tree.hh>

#include <dune/typetree/leafnode.hh>
#include <dune/typetree/dynamicpowernode.hh>

#include <utility>
#include <span>
#include <memory>

namespace Dune::PDELab::inline Experimental {

  template<Concept::MultiIndex ContainerIndex, Dune::Concept::GridView EntitySet>
  struct EmptyConstraintsContainer : public TypeTree::LeafNode {
    EmptyConstraintsContainer(const EntitySet&) {}

    [[nodiscard]] static constexpr std::false_type isConstrained(Concept::MultiIndex auto ci) { return {}; }

    void applyAffineConstraints(auto& x) const;

    void applyLinearConstraints(auto& x) const;

    void applyLinearTransposedConstraints(auto& x) const;

    void clear() {}

    static std::size_t size() { return 0; }

    void globalCompress(Concept::Basis auto space) {}
    void localCompress(Concept::LocalBasisTree auto& lbasis_tree) {}

    class LeafLocalView : public TypeTree::LeafNode {
    public:
      using size_type = std::size_t;

      LeafLocalView() {}

      [[nodiscard]] static std::false_type isConstrained(size_type dof) { return {}; }

      [[nodiscard]] static auto linearCoefficients(size_type dof) noexcept {
        return std::span<const std::pair<ContainerIndex,double>,0>{};
      }

      [[nodiscard]] static double translationValue(size_type dof) noexcept {
        return double{0.};
      }

      void bind(const Dune::Concept::Entity auto& entity) const {}

      void unbind() {}
    };

    class SubEntityLocalView : public TypeTree::DynamicPowerNode<LeafLocalView> {
      using Base = TypeTree::DynamicPowerNode<LeafLocalView>;
    public:

      SubEntityLocalView(typename Base::NodeStorage&& storage)
        : Base{ std::move(storage) }
      {}
    };


    template<Concept::Tree Tree>
    static auto makeLocalViewNode(const Tree& tree, const std::shared_ptr<const EmptyConstraintsContainer>& container) {
      if constexpr (Concept::LeafTreeNode<Tree>) {
        return std::make_shared<LeafLocalView>();
      } else { // skeleton finite element case
        static_assert(Concept::VectorTreeNode<Tree> and Concept::LeafTreeNode<typename Tree::ChildType>);
        typename SubEntityLocalView::NodeStorage storage(tree.degree());
        for (std::size_t i = 0; i != storage.size(); ++i)
          storage[i] = std::make_shared<LeafLocalView>();
        return std::make_shared<SubEntityLocalView>(std::move(storage));
      }
    }
  };


} // namespace Dune::PDELab::inline Experimental

#endif // DUNE_PDELAB_BASIS_CONSTRAINTS_CONTAINER_EMPTY_HH
