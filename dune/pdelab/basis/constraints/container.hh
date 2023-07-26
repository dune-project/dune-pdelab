#ifndef DUNE_PDELAB_BASIS_CONSTRAINTS_CONTAINER_HH
#define DUNE_PDELAB_BASIS_CONSTRAINTS_CONTAINER_HH

#include <dune/pdelab/concepts/treenode.hh>
#include <dune/pdelab/concepts/local_basis.hh>
#include <dune/pdelab/concepts/basis.hh>
#include <dune/pdelab/concepts/multiindex.hh>

#include <dune/pdelab/common/tree_traversal.hh>

#include <dune/pdelab/basis/prebasis/concept.hh>

#include <dune/typetree/dynamicpowernode.hh>
#include <dune/typetree/powernode.hh>
#include <dune/typetree/compositenode.hh>
#include <dune/grid/concepts/entity.hh>

#include <utility>
#include <span>

namespace Dune::PDELab::inline Experimental {

  template<Concept::TreeNode T>
  class VectorLocalConstraintsContainer : public TypeTree::DynamicPowerNode<T> {
    using TreeNode = TypeTree::DynamicPowerNode<T>;
  public:
    VectorLocalConstraintsContainer(const TreeNode::NodeStorage& storage)
      : TreeNode{storage}
    {}

    VectorLocalConstraintsContainer(const VectorLocalConstraintsContainer&) = delete;
    VectorLocalConstraintsContainer(VectorLocalConstraintsContainer&&) = delete;

    VectorLocalConstraintsContainer& operator=(const VectorLocalConstraintsContainer&) = delete;
    VectorLocalConstraintsContainer& operator=(VectorLocalConstraintsContainer&&) = default;
  };

  template<Concept::TreeNode T>
  class VectorConstraintsContainer : public TypeTree::DynamicPowerNode<T>
  {
    using TreeNode = TypeTree::DynamicPowerNode<T>;

  public:
    VectorConstraintsContainer(const typename TreeNode::NodeStorage& storage)
      : TreeNode{ storage }
    {}

    template<Concept::Tree Tree>
    static auto makeLocalViewNode(const Tree& tree, std::shared_ptr<const VectorConstraintsContainer> container) {
      static_assert(Concept::VectorTreeNode<Tree>);
      using ChildNode = std::decay_t<decltype(*container->child(0).makeLocalViewNode(tree.child(0), container->childStorage(0)))>;
      using LocalView = VectorLocalConstraintsContainer<ChildNode>;
      typename LocalView::NodeStorage storage(container->degree());
      for (std::size_t i = 0; i != storage.size(); ++i)
        storage[i] = container->child(i).makeLocalViewNode(tree.child(i), container->childStorage(i));
      return std::make_shared<LocalView>(std::move(storage));
    }
  };


  template<Concept::TreeNode T, std::size_t k>
  class ArrayLocalConstraintsContainer : public TypeTree::PowerNode<T,k> {
    using TreeNode = TypeTree::PowerNode<T,k>;
  public:
    ArrayLocalConstraintsContainer(const TreeNode::NodeStorage& storage)
      : TreeNode{storage}
    {}

    ArrayLocalConstraintsContainer(const ArrayLocalConstraintsContainer&) = delete;
    ArrayLocalConstraintsContainer(ArrayLocalConstraintsContainer&&) = delete;

    ArrayLocalConstraintsContainer& operator=(const ArrayLocalConstraintsContainer&) = delete;
    ArrayLocalConstraintsContainer& operator=(ArrayLocalConstraintsContainer&&) = default;
  };

  template<Concept::TreeNode T, std::size_t k>
  class ArrayConstraintsContainer : public TypeTree::PowerNode<T, k>
  {
    using TreeNode = TypeTree::PowerNode<T, k>;

  public:
    ArrayConstraintsContainer(const typename TreeNode::NodeStorage& storage)
      : TreeNode{ storage }
    {}

    template<Concept::Tree Tree>
    static auto makeLocalViewNode(const Tree& tree, std::shared_ptr<const ArrayConstraintsContainer> container) {
      static_assert(Concept::ArrayTreeNode<Tree>);
      using ChildNode = std::decay_t<decltype(*container->child(0).makeLocalViewNode(tree.child(0), container->childStorage(0)))>;
      using LocalView = ArrayLocalConstraintsContainer<ChildNode,k>;
      typename LocalView::NodeStorage storage;
      for (std::size_t i = 0; i != storage.size(); ++i)
        storage[i] = container->child(i).makeLocalViewNode(tree.child(i), container->childStorage(i));
      return std::make_shared<LocalView>(std::move(storage));
    }
  };


  template<Concept::TreeNode... T>
  class TupleLocalConstraintsContainer : public TypeTree::CompositeNode<T...> {
    using TreeNode = TypeTree::CompositeNode<T...>;
  public:
    TupleLocalConstraintsContainer(const TreeNode::NodeStorage& storage)
      : TreeNode{storage}
    {}

    TupleLocalConstraintsContainer(const TupleLocalConstraintsContainer&) = delete;
    TupleLocalConstraintsContainer(TupleLocalConstraintsContainer&&) = delete;

    TupleLocalConstraintsContainer& operator=(const TupleLocalConstraintsContainer&) = delete;
    TupleLocalConstraintsContainer& operator=(TupleLocalConstraintsContainer&&) = default;
  };

  template<Concept::TreeNode... T>
  class TupleConstraintsContainer : public TypeTree::CompositeNode<T...>
  {
    using TreeNode = TypeTree::CompositeNode<T...>;

  public:
    TupleConstraintsContainer(const typename TreeNode::NodeStorage& storage)
      : TreeNode{ storage }
    {}

    template<Concept::Tree Tree>
    static auto makeLocalViewNode(const Tree& tree, std::shared_ptr<const TupleConstraintsContainer> container) {
      static_assert(Concept::TupleTreeNode<Tree>);
      auto unfold_children = [&](auto... i) {
        using LocalView = TupleLocalConstraintsContainer<std::decay_t<decltype(*container->child(i).makeLocalViewNode(tree.child(i), container->childStorage(i)))>...>;
        typename LocalView::NodeStorage storage{ container->child(i).makeLocalViewNode(tree.child(i), container->childStorage(i))... };
        return std::make_shared<LocalView>(std::move(storage));
      };
      return unpackIntegerSequence(unfold_children, std::make_index_sequence<sizeof...(T)>{});
    }
  };


  template<Concept::Tree ConstraintsContainerTree>
  class ConstraintsContainer {
  public:

    // using Tree = ConstraintsContainerTree;

    ConstraintsContainer(std::shared_ptr<ConstraintsContainerTree> tree)
      : _tree{std::move(tree)}
    {}

    ConstraintsContainerTree& tree() {
      assert(not _assembled);
      return *_tree;
    }

    template<Concept::Tree SourceTree, Concept::MultiIndex SubBasisPath>
    class LocalView {

      static auto makeLocalTree(const Concept::Tree auto& source_tree, const std::shared_ptr<ConstraintsContainerTree>& gtree, SubBasisPath sub_basis_path) {
        // source tree is already on the sub space path, whereas the gtree is the root node of the tree
        // source tree is only to give the same tree structure to this tree, its contents (other than degrees) are ignored. Thus, once constructed, source_tree may be as well this->tree()
        if constexpr (SubBasisPath::size() == 0)
          return gtree->makeLocalViewNode(source_tree, gtree);
        else
          return TypeTree::childStorage(*gtree, sub_basis_path)->makeLocalViewNode(source_tree, TypeTree::childStorage(*gtree, sub_basis_path));
      }

      using TreeStorage = decltype(makeLocalTree(std::declval<SourceTree>(), std::shared_ptr<ConstraintsContainerTree>{}, SubBasisPath{}));
    public:
      using Tree = typename TreeStorage::element_type;
      static_assert(std::is_same_v<TreeStorage, decltype(makeLocalTree(std::declval<Tree>(), std::shared_ptr<ConstraintsContainerTree>{}, SubBasisPath{}))>);

      LocalView(const SourceTree& source_tree,
                std::shared_ptr<ConstraintsContainerTree> gtree,
                SubBasisPath sub_basis_path)
        : _sub_basis_path{ sub_basis_path }
        , _gtree{ std::move(gtree) }
        , _ltree{ makeLocalTree(source_tree, _gtree, _sub_basis_path) }
      {}

      LocalView(const LocalView& other) {
        (*this) = other;
      }

      LocalView(LocalView&&) = default;

      LocalView& operator=(const LocalView& other) {
        _ltree = makeLocalTree(other.tree(), _gtree = other._gtree, _sub_basis_path = other._sub_basis_path);
        return *this;
      }

      LocalView& operator=(LocalView&&) = default;

      const Tree& tree() const {
        return *_ltree;
      }

      void bind(const Dune::Concept::Entity auto& entity) {
        // _padlocks_view = _constraints_container.getEntityLockHandlerSpan(entity);
        forEachLeafNode(*_ltree, [&](auto& node){
          node.bind(entity);
        });
      }

      void unbind() {
        // _padlocks_view = std::span<LockHandle>{};
        forEachLeafNode(*_ltree, [&](auto& node){
          node.unbind();
        });
      }

//   // lock all the constrained dofs non-local to the currently bound entity
//   // notice that dofs local to the entity should be locked by the local space
//   void lock() {
//     while (not try_lock()) [[unlikely]] {
//       __builtin_ia32_pause();
//     }
//   }

//   [[nodiscard]] bool try_lock() {
//     for (std::size_t i = 0; i != _padlocks_view.size(); ++i) {
//       // try to lock every padlock in the vector
//       if (not _padlocks_view[i].try_lock()) [[unlikely]] {
//         // padlock was already locked, we have to roll back
//         for (std::size_t j = i; j != 0; --j)
//           _padlocks_view[j-1].unlock();
//         // ...and inform that we could not adquire the lock
//         return false;
//       }
//     }
//     // if all padlocks were locked, we succeded
//     return true;
//   }

//   void unlock() {
//     for (auto& lock : _padlocks_view)
//       lock.unlock();
//   }

    private:
      [[no_unique_address]] SubBasisPath _sub_basis_path;
      std::shared_ptr<ConstraintsContainerTree> _gtree;
      std::shared_ptr<Tree> _ltree;
      //   std::span<LockHandle> _padlocks_view;
    };

    template<Concept::Tree Tree, Concept::MultiIndex SubBasisPath>
    Concept::LocalConstraints auto localView(const Tree& tree, SubBasisPath sub_basis_path) const {
      assert(_assembled);
      return LocalView<Tree, SubBasisPath>{ tree, _tree, sub_basis_path };
    }

    template<Concept::Tree Tree>
    Concept::LocalConstraints auto localView(const Tree& tree) const {
      return localView(tree, TypeTree::treePath());
    }

    void assembleConstraints(const Concept::Basis auto& basis, auto constraints_ops) {
      _assembled = false;
      // this should only be done with root nodes
      auto lbasis_in = basis.localView();
      auto lbasis_out = basis.localView();

      bool constrained = false;
      bool intersection_constrained = false;
      this->clear();
      forEachLeafNode(this->tree(), [&](const auto& container_node, auto path){
        intersection_constrained |= constraints_ops[path].doConstrainSkeleton() | constraints_ops[path].doConstrainBoundary();
        constrained |= intersection_constrained | constraints_ops[path].doConstrainVolume();
      });

      _assembled = not constrained;
      if (_assembled) return;

      for (const auto& entity : elements(basis.entitySet())) {
        lbasis_in.bind(entity);
        forEachLeafNode(this->tree(), [&](auto& container_node, auto path){
          auto& constraints_node = constraints_ops[path];
          const auto& lbasis_in_node = PDELab::containerEntry(lbasis_in.tree(), path);
          if (constraints_node.doConstrainVolume())
            constraints_node.constrainVolume(lbasis_in_node, container_node);
        });
        if (intersection_constrained) {
          // notice that there is double visit, this can be optimized if necessary
          for (const auto& intersection : intersections(basis.entitySet(), entity)) {
            if (intersection.neighbor()) {
              lbasis_out.bind(intersection.outside());
              forEachLeafNode(this->tree(), [&](auto& container_node, auto path){
                auto& constraints_node = constraints_ops[path];
                const auto& lbasis_in_node = PDELab::containerEntry(lbasis_in.tree(), path);
                const auto& lbasis_out_node = PDELab::containerEntry(lbasis_out.tree(), path);
                if (constraints_node.doConstrainSkeleton())
                  constraints_node.constrainSkeleton(intersection, lbasis_in_node, lbasis_out_node, container_node);
              });
              lbasis_out.unbind();
            } else {
              forEachLeafNode(this->tree(), [&](auto& container_node, auto path){
                auto& constraints_node = constraints_ops[path];
                const auto& lbasis_in_node = PDELab::containerEntry(lbasis_in.tree(), path);
                if (constraints_node.doConstrainBoundary())
                  constraints_node.constrainBoundary(intersection, lbasis_in_node, container_node);
              });
            }
          }
        }
        lbasis_in.unbind();
      }

      this->compress(basis);
    }

    void clear() {
      forEachLeafNode(*_tree, [](auto& node){
        node.clear();
      });
      _assembled = false;
    }

    void compress(Concept::Basis auto basis) {
      // basis must be a root basis to map all entities and all sub basiss(i.e. no sub-basis)
      forEachLeafNode(*_tree, [&](auto& node){
        node.globalCompress(basis);
      });
      auto lbasis = basis.localView();
      for (const auto& entity : elements(basis.entitySet())) {
        lbasis.bind(entity);
        forEachLeafNode(*_tree, [&](auto& node, auto path){
          const auto& lbasis_node = PDELab::containerEntry(lbasis.tree(), path);
          node.localCompress(lbasis_node);
        });
        lbasis.unbind();
      }
      _assembled = true;
    }

    // std::vector<std::size_t> _padlocks_entity_offset;
    // std::vector<LockHandle> _padlocks_entity;

    bool _assembled = false;
    std::shared_ptr<ConstraintsContainerTree> _tree;
  };

template<Concept::Tree Tree, Concept::MultiIndex Path, class Callable>
auto makeConstraintsContainerNode(const Tree& tree, Path path, Callable callable)
{
  if constexpr (Concept::LeafTreeNode<Tree>) {
    static_assert(std::invocable<Callable, Tree, Path>);
    return callable(tree, path);
  } else if constexpr (Concept::ArrayTreeNode<Tree>) {
    using ChidlNode = std::decay_t<decltype(*makeConstraintsContainerNode(tree.child(0),push_back(path,0), callable))>;
    using Node = ArrayConstraintsContainer<ChidlNode, Tree::degree()>;
    typename Node::NodeStorage storage;
    for (std::size_t i = 0; i < tree.degree(); ++i)
      storage[i] = makeConstraintsContainerNode(tree.child(i), push_back(path,i), callable);
    return std::make_shared<Node>(storage);
  } else if constexpr (Concept::VectorTreeNode<Tree>) {
    using ChidlNode = std::decay_t<decltype(*makeConstraintsContainerNode(tree.child(0), push_back(path,0), callable))>;
    using Node = VectorConstraintsContainer<ChidlNode>;
    typename Node::NodeStorage storage(tree.degree());
    for (std::size_t i = 0; i < tree.degree(); ++i)
      storage[i] = makeConstraintsContainerNode(tree.child(i), push_back(path,i), callable);
    return std::make_shared<Node>(storage);
  } else {
    static_assert(Concept::TupleTreeNode<Tree>);
    return unpackIntegerSequence(
      [&](auto... i) {
        using ChildNodes = std::tuple<std::decay_t<decltype(*makeConstraintsContainerNode(tree.child(i), push_back(path,i), callable))>...>;
        using Node = TupleConstraintsContainer<std::tuple_element_t<i, ChildNodes>...>;
        typename Node::NodeStorage storage{ makeConstraintsContainerNode(tree.child(i), push_back(path,i), callable)... };
        return std::make_shared<Node>(storage);
      },
      std::make_index_sequence<Tree::degree()>{});
  }
}

// for a given tree and a callable generator `callable(leaf_node, path)`, build a constrains container tree
auto makeConstraintsContainer(Concept::Tree auto& tree, auto callable)
{
  auto tree_ptr = makeConstraintsContainerNode(tree, TypeTree::treePath(), callable);
  using ConstraintsContainerTree = decltype(tree_ptr)::element_type;
  static_assert(Concept::Tree<ConstraintsContainerTree>);
  return std::make_shared<ConstraintsContainer<ConstraintsContainerTree>>(tree_ptr);
}


} // namespace Dune::PDELab::inline Experimental

#endif // DUNE_PDELAB_BASIS_CONSTRAINTS_CONTAINER_HH
