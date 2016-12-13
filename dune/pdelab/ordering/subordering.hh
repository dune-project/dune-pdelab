// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_ORDERING_SUBORDERING_HH
#define DUNE_PDELAB_ORDERING_SUBORDERING_HH

#include <algorithm>
#include <iterator>
#include <memory>

#include <dune/common/array.hh>

#include <dune/typetree/treepath.hh>
#include <dune/typetree/proxynode.hh>
#include <dune/typetree/childextraction.hh>

#include <dune/pdelab/ordering/utility.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup Ordering
    //! \{

    //! A view on a subtree of a multi-component ordering.
    /**
     * SubOrdering presents a read-only view on the components of BaseOrdering contained
     * in the subtree pointed at by TreePath. The TreePath has to be an instantiation of
     * Dune::TypeTree::TreePath and is interpreted from left to right, i.e. a
     * TypeTree<3,1> means take the fourth child of BaseOrdering and continue to the second
     * child of that ordering.
     *
     * \note SubOrdering is not (yet) as feature-complete as a regular ordering. The most
     *       important missing feature is support for per-entity indices.
     *
     * \warning SubOrdering does *not* support nesting! Trying to construct a Subordering
     *          with a SubOrdering as BaseOrdering will not work and might either just fail
     *          to compile or fail in subtle ways at runtime, so don't try it. SubOrdering
     *          is not really an ordering, but more of a view into a subset of a real Ordering.
     *          Overall, this restriction should not be problematic, as users will usually
     *          construct GridFunctionSubSpaces, which take care of avoiding this recursion.
     *
     * \tparam BaseOrdering_  The type of the underlying base ordering.
     * \tparam TreePath       The path to the subtree for which to provide a view.
     *                        \warning This path will usually differ from the path into the
     *                                 GridFunctionSpace!
     */
    template<typename BaseOrdering_, typename TreePath>
    class SubOrdering
      : public TypeTree::ProxyNode<const TypeTree::ChildForTreePath<BaseOrdering_,TreePath>>
    {

      using NodeT = TypeTree::ProxyNode<const TypeTree::ChildForTreePath<BaseOrdering_,TreePath>>;

    public:

      //! The type of the BaseOrdering for which to represent a SubOrdering view.
      typedef BaseOrdering_ BaseOrdering;

      //! The target ordering that makes up the root of this SubOrdering view.
      using TargetOrdering = TypeTree::ChildForTreePath<BaseOrdering,TreePath>;

      //! Forwarded Ordering traits from BaseOrdering.
      typedef typename BaseOrdering::Traits Traits;

      //! Forwarded tag from BaseOrdering, required by PDELab internals.
      typedef typename BaseOrdering::ContainerAllocationTag ContainerAllocationTag;

      //! Forwarded tag from BaseOrdering, required by PDELab internals.
      typedef typename BaseOrdering::CacheTag CacheTag;

      //! Forwarded ordering property from TargetOrdering, required by PDELab internals.
      static const bool has_dynamic_ordering_children = TargetOrdering::has_dynamic_ordering_children;

      //! Forwarded ordering property from TargetOrdering, required by PDELab internals.
      static const bool consume_tree_index = TargetOrdering::consume_tree_index;


      //! Constructs a SubOrdering for base_ordering.
      /**
       * Constructor for a subordering view of base_ordering rooted in the subtree pointed to by
       * TreePath. The second parameter DOFIndexTreePath represents the subtree in the DOFIndex
       * space, i.e. the path to append to passed-in DOFIndices to obtain valid DOFIndices that
       * can be mapped using base_ordering (usually, this will be the TreePath of a GridFunctionSubSpace).
       * This information has to be passed in separately, as the two paths will generally not be
       * identical due to synthesized nodes in the ordering tree.
       *
       * \param base_ordering      Storage pointer to the BaseOrdering.
       * \tparam DOFIndexTreePath  Instantiation of TypeTree::TreePath, representing the path from
       *                           the root of the DOFIndex tree to the DOFIndices passed to this
       *                           ordering.
       */
      explicit SubOrdering(std::shared_ptr<const BaseOrdering> base_ordering)
        : NodeT(base_ordering->child(TreePath()))
        , _base_ordering(base_ordering)
      {
        update();
      }

      //! Updates this SubOrdering.
      void update()
      {}


      template<typename ItIn, typename ItOut>
      void map_lfs_indices(ItIn begin, const ItIn end, ItOut out) const
      {
        // Do the mapping up to the root ordering.
        // Avoid spelling out the type of ItIn here (it has to be a DOFIndexViewIterator),
        // so we don't need to include lfsindexcache.hh.
        map_lfs_indices_to_root_space(TreePath(),
                                      begin,
                                      end,
                                      out);
      }


    private:


      //! Performs the single-level index mapping in the ordering pointed at by TP.
      template<typename TP, typename ItIn, typename ItOut>
      void map_lfs_indices_in_ancestor(TP tp, ItIn& begin, ItIn& end, ItOut out) const
      {
        using Ordering = TypeTree::ChildForTreePath<BaseOrdering,TP>;

        // This logic needs to be replicated from the IndexCache visitor, as we bypass
        // the tree-visiting algorithm and work our way up the tree all by ourselves.
        // Don't consume the first entry in the tree path to the parent before it has
        // been used!
        if (!std::is_same<TreePath,TP>::value && Ordering::consume_tree_index)
          {
            begin.restore_back();
            end.restore_back();
          }

        // Call the single-level mapping step of our ancestor ordering.
        TypeTree::child(baseOrdering(),tp).map_lfs_indices(begin,end,out);
      }

      // Template recursion for walking up the TreePath to the BaseOrdering
      template<typename TP, typename ItIn, typename ItOut>
      typename std::enable_if<
        (TypeTree::TreePathSize<TP>::value > 0)
          >::type
      map_lfs_indices_to_root_space(TP, ItIn begin, ItIn end, ItOut out) const
      {
        map_lfs_indices_in_ancestor(TP(),begin,end,out);
        // recurse further up to the tree
        map_lfs_indices_to_root_space(typename TypeTree::TreePathPopBack<TP>::type(),begin,end,out);
      }

      // End of template recursion for walking up the TreePath to the BaseOrdering
      template<typename TP, typename ItIn, typename ItOut>
      typename std::enable_if<
        (TypeTree::TreePathSize<TP>::value == 0)
          >::type
      map_lfs_indices_to_root_space(TP, ItIn begin, ItIn end, ItOut out) const
      {
        map_lfs_indices_in_ancestor(TP(),begin,end,out);
      }

    public:

      //! Maps di from the DOFIndex subtree to the ContainerIndex in the BaseOrdering.
      typename Traits::ContainerIndex mapIndex(const typename Traits::DOFIndex& di) const
      {
        typename Traits::ContainerIndex ci;
        mapIndex(di,ci);
        return ci;
      }

      //! Maps di from the DOFIndex subtree to the ContainerIndex in the BaseOrdering - inplace version.
      void mapIndex(typename Traits::DOFIndexView di, typename Traits::ContainerIndex& ci) const
      {
        baseOrdering().mapIndex(di,ci);
      }

      //! Returns the size of the BaseOrdering.
      typename Traits::SizeType size() const
      {
        return baseOrdering().size();
      }

      //! Returns the block count of the BaseOrdering.
      typename Traits::SizeType blockCount() const
      {
        return baseOrdering().blockCount();
      }

      //! Returns the maximum per-entity size of the TargetOrdering.
      typename Traits::SizeType maxLocalSize() const
      {
        return targetOrdering().maxLocalSize();
      }

      //! Returns whether the TargetOrdering has DOFs attached to entities of codimension codim.
      bool contains(typename Traits::SizeType codim) const
      {
        return targetOrdering().contains(codim);
      }

      //! Returns whether the TargetOrdering is of fixed size for entities of codimension codim.
      bool fixedSize(typename Traits::SizeType codim) const
      {
        return targetOrdering().fixedSize(codim);
      }

      //! Returns the BaseOrdering.
      const BaseOrdering& baseOrdering() const
      {
        return *_base_ordering;
      }

      //! Returns the TargetOrdering.
      const TargetOrdering& targetOrdering() const
      {
        return this->proxiedNode();
      }

    private:

      std::shared_ptr<const BaseOrdering> _base_ordering;

    };

    //! \}

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ORDERING_SUBORDERING_HH
