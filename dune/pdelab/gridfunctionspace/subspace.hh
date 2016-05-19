// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_SUBSPACE_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_SUBSPACE_HH

/** \file
 *  \brief GridFunctionSubSpace implementation.
 */

#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/powergridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/compositegridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/subspacelocalfunctionspace.hh>
#include <dune/pdelab/ordering/subordering.hh>

namespace Dune {
  namespace PDELab {

    namespace gfs {

      // forward declaration for use in build_dof_index_type specialization and
      // in feature mixins.
      template<typename GFS, typename TreePath>
      class GridFunctionSubSpace;

    } // namespace gfs

#ifndef DOXYGEN

    // Specialization of DOFIndex type deduction TMP - the DOFIndex
    // of a subspace must be large enough to contain DOFIndex values
    // for the complete tree rooted in the base space.
    template<typename GFS, typename TP>
    struct build_dof_index_type<gfs::GridFunctionSubSpace<GFS,TP> >
    {
      typedef typename GFS::Ordering::Traits::DOFIndex type;
    };

#endif // DOXYGEN

    //! \addtogroup GridFunctionSpace grid function space
    //! \ingroup PDELab
    //! \{


    //! Implementation namespace for GridFunctionSpace-specific features.
    namespace gfs {

      namespace {



        // ********************************************************************************
        // Helper TMPs
        // ********************************************************************************

        //! TMP for deducing the TreePath to the ordering subtree corresponding to a GridFunctionSubSpace
        /**
         * This TMP will recursively walk down the GFS and the Ordering tree and insert additional
         * entries into the Ordering tree as required.
         */
        template<typename Ordering, typename GFS, typename GFSTP, typename OrderingTP = TypeTree::TreePath<> >
        struct find_ordering_treepath_for_sub_gfs
        {

          // Get the ordering at the current subtree position.
          using SubOrdering = TypeTree::ChildForTreePath<Ordering,OrderingTP>;

          // Only descend in the GFS tree if the current ordering child consumes a tree index entry.
          typedef typename std::conditional<
            SubOrdering::consume_tree_index,
            typename GFS::template Child<TypeTree::TreePathFront<GFSTP>::value>::type,
            GFS
            >::type SubGFS;

          // Insert either GFS child index or synthesized child index (always 0) in the ordering treepath.
          typedef typename TypeTree::TreePathPushBack<
            OrderingTP,
            (SubOrdering::consume_tree_index ? TypeTree::TreePathFront<GFSTP>::value : 0)
            >::type SubOrderingTP;

          // Keep (synthesized ordering node) or drop (ordering with associated GFS) first entry of GFS TreePath.
          typedef typename std::conditional<
            SubOrdering::consume_tree_index,
            typename TypeTree::TreePathPopFront<GFSTP>::type,
            GFSTP
            >::type SubGFSTP;

          // Recurse into child trees.
          typedef typename find_ordering_treepath_for_sub_gfs<
            Ordering,
            SubGFS,
            SubGFSTP,
            SubOrderingTP
            >::type type;

        };

        //! End of recursion for TreePath-deducing TMP.
        template<typename Ordering, typename GFS, typename OrderingTP>
        struct find_ordering_treepath_for_sub_gfs<Ordering,GFS,TypeTree::TreePath<>,OrderingTP>
        {

          // We have found the correct ordering TreePath, so let's return it.
          typedef OrderingTP type;

        };

      } // anonymous namespace




      // *****************************************************************************************
      // Feature provider mixins
      // *****************************************************************************************

      //! Default features used by every subspace implementation.
      template<typename GFS, typename TreePath, typename Tag>
      class DefaultSubSpaceFeatures
      {

        //! CRTP typedef
        typedef GridFunctionSubSpace<GFS,TreePath> SubSpace;

        //! CRTP access to actual class
        const SubSpace& subSpace() const
        {
          return static_cast<const SubSpace&>(*this);
        }

      public:

        //! \name Default Functionality for all GridFunctionSpaces
        //! \{

        //! The TreePath from the root of the space hierarchy to this subspace.
        typedef TreePath SubSpacePath;

        //! The base GridFunctionSpace that this GridFunctionSubSpace is based on.
        typedef GFS BaseGridFunctionSpace;

        //! The type of the original GridFunctionSpace that is the root of this GridFunctionSpace.
        using ChildGridFunctionSpace = TypeTree::ChildForTreePath<GFS,TreePath>;

        //! Re-exported Traits from the original GridFunctionSpace.
        typedef typename ChildGridFunctionSpace::Traits Traits;

        //! Re-exported OrderingTag from the original GridFunctionSpace.
        typedef typename ChildGridFunctionSpace::OrderingTag OrderingTag;


        //! Re-exported constraints container from the original GridFunctionSpace.
        template<typename E>
        using Constraintscontainer = typename GFS::template ConstraintsContainer<E>;

        //! The ordering used by this GridFunctionSubSpace.
        typedef SubOrdering<
          typename GFS::Ordering,
          typename find_ordering_treepath_for_sub_gfs<
            typename GFS::Ordering,
            GFS,
            TreePath
            >::type
          > Ordering;

        std::size_t subSpaceDepth() const
        {
          return TypeTree::TreePathSize<SubSpacePath>::value;
        }

        //! Returns the ordering associated with this GridFunctionSubSpace.
        const Ordering& ordering() const
        {
          return _ordering;
        }

        //! Returns the underlying EntitySet.
        const typename Traits::EntitySet& entitySet() const
        {
          return subSpace().childGridFunctionSpace().entitySet();
        }

        //! Returns the underlying GridView.
        const typename Traits::GridViewType& gridView() const
        {
          return subSpace().childGridFunctionSpace().gridView();
        }

        //! Returns the global size of the root space.
        typename Traits::SizeType globalSize() const
        {
          return _ordering.size();
        }

        //! Returns the global size of the root space.
        /**
         * \warning The semantics of this methods have changed with the introduction
         *          of Orderings: While this method used to return the size of the subspace
         *          only, it now behaves like globalSize() and returns the overall size of
         *          the root space! Calculating the size of the subspace might be a very
         *          expensive operation depending on the underlying orderings.
         */
        typename Traits::SizeType size() const
        {
          return _ordering.size();
        }

        //! Returns the maximum number of DOFs per cells in this subspace.
        typename Traits::SizeType maxLocalSize() const
        {
          return _ordering.maxLocalSize();
        }

        //! \}

      protected:

        DefaultSubSpaceFeatures(const GFS& gfs)
          : _ordering(gfs.orderingStorage())
        {}

      private:

        Ordering _ordering;

      };


      //! Additional features used by leaf subspaces.
      template<typename GFS, typename TreePath, typename Tag>
      class LeafSubSpaceFeatures
      {

        //! CRTP typedef
        typedef GridFunctionSubSpace<GFS,TreePath> SubSpace;

        //! CRTP access
        const SubSpace& subSpace() const
        {
          return static_cast<const SubSpace&>(*this);
        }

      public:

        //! The type of the original GridFunctionSpace that is the root of this GridFunctionSpace.
        using ChildGridFunctionSpace = TypeTree::ChildForTreePath<GFS,TreePath>;

        //! Re-exported Traits from the original GridFunctionSpace.
        typedef typename ChildGridFunctionSpace::Traits Traits;

        //! \name Additional Functionality for Leaf Spaces
        //! \{

        //! Returns the finite element map of this space.
        const typename Traits::FiniteElementMap& finiteElementMap() const
        {
          return subSpace().childGridFunctionSpace().finiteElementMap();
        }

        //! Returns the storage object for the finite element map of this space.
        std::shared_ptr<const typename Traits::FiniteElementMap> finiteElementMapStorage () const
        {
          return subSpace().childGridFunctionSpace().finiteElementMapStorage();
        }

        //! Returns the constraints engine of this space.
        const typename Traits::ConstraintsType& constraints() const
        {
          return subSpace().childGridFunctionSpace().constraints();
        }

        //! Returns the name of this space.
        const std::string& name() const
        {
          return subSpace().childGridFunctionSpace().name();
        }
        //! \}

      };


#ifdef DOXYGEN


      //! Feature provider for leaf spaces.
      //! \nosubgrouping
      template<typename GFS, typename TreePath, typename Tag>
      class SubSpaceFeatureProvider
        : public DefaultSubSpaceFeatures<GFS,TreePath,Tag>
        , public LeafSubSpaceFeatures<GFS,TreePath,Tag>
      {

      protected:

        SubSpaceFeatureProvider(const GFS& gfs)
          : DefaultSubSpaceFeatures<GFS,TreePath,Tag>(gfs)
        {}

      };

#else // DOXYGEN

      //! Default subspace feature provider.
      template<typename GFS, typename TreePath, typename Tag>
      class SubSpaceFeatureProvider
        : public DefaultSubSpaceFeatures<GFS,TreePath,Tag>
      {

      protected:

        SubSpaceFeatureProvider(const GFS& gfs)
          : DefaultSubSpaceFeatures<GFS,TreePath,Tag>(gfs)
        {}

      };

      //! Feature provider for leaf spaces.
      template<typename GFS, typename TreePath>
      class SubSpaceFeatureProvider<GFS,TreePath,LeafGridFunctionSpaceTag>
        : public DefaultSubSpaceFeatures<GFS,TreePath,LeafGridFunctionSpaceTag>
        , public LeafSubSpaceFeatures<GFS,TreePath,LeafGridFunctionSpaceTag>
      {

      protected:

        SubSpaceFeatureProvider(const GFS& gfs)
          : DefaultSubSpaceFeatures<GFS,TreePath,LeafGridFunctionSpaceTag>(gfs)
        {}

      };

#endif // DOXYGEN

      //! Mixin class which inherits from GridFunctionOutputParameters iff T inherits from GridFunctionOutputParameters
      template<class T, bool Enable = std::is_base_of<GridFunctionOutputParameters, T>::value>
      class GridFunctionSubSpaceOutputParameters
      {
      public:
        void inheritDataSetType(const T & t) {}
      };

#ifndef DOXYGEN
      template<class T>
      class GridFunctionSubSpaceOutputParameters<T,true> :
        public GridFunctionOutputParameters
      {
      public:
        void inheritDataSetType(const T & t)
        {
          setDataSetType(t.dataSetType());
        }
      };
#endif

      // ********************************************************************************
      // GridFunctionSubSpace implementation
      // ********************************************************************************

      //! Non-nesting implementation of GridFunctionSubSpace.
      /**
       * This is the actual implementation of GridFunctionSubSpace. It is based around the idea
       * of performing the mapping from the subspace to the rootspace in a single step using a
       * SubOrdering. As SubOrderings cannot be nested, this class needs some helper functionality
       * that constructs a new GridFunctionSubSpace directly on top of the underlying root GridFunctionSpace
       * if the user attempts to nest GridFunctionSubSpaces. On the other hand, this implementation should
       * render such usage mostly unnecessary, as it is now possible to directly construct a subspace for
       * a given leaf space.
       * If the compiler has support for template aliases, the de-nesting infrastructure will be able to
       * completely remove all traces of nesting, letting the type obtained by nesting two GridFunctionSubSpaces
       * look exactly like the type obtained by constructing a single GridFunctionSubSpace.
       * Alternatively, without template aliases, it is only possible to have the actual implementation be
       * non-nested. In this case, the Dune::PDELab::GridFunctionSubSpace classes will still nest.
       *
       * \note This class should always be used as Dune::PDELab::GridFunctionSubSpace. Never attempt
       *       to directly use the class Dune::PDELab::gfs::GridFunctionSubSpace!
       *
       * \tparam GFS       The root GridFunctionSpace.
       * \tparam TreePath  Path from the root GridFunctionSpace to the represented subspace.
       *
       * \nosubgrouping
       */
      template<typename GFS, typename TreePath>
      class GridFunctionSubSpace
        : public TypeTree::ProxyNode<const TypeTree::ChildForTreePath<GFS,TreePath>>
        , public SubSpaceFeatureProvider<GFS,TreePath,TypeTree::ImplementationTag<TypeTree::ChildForTreePath<GFS,TreePath>>>
        , public GridFunctionSubSpaceOutputParameters<TypeTree::ChildForTreePath<GFS,TreePath>>
      {

        using NodeT = TypeTree::ProxyNode<const TypeTree::ChildForTreePath<GFS,TreePath>>;

        using FeatureT = SubSpaceFeatureProvider<
          GFS,
          TreePath,
          TypeTree::ImplementationTag<TypeTree::ChildForTreePath<GFS,TreePath>>
          >;

      public:

        //! Construct a GridFunctionSubSpace from the storage object of a root space.
        explicit GridFunctionSubSpace(std::shared_ptr<const GFS> gfs_storage)
          : NodeT(TypeTree::childStorage(*gfs_storage,TreePath()))
          , FeatureT(*gfs_storage)
          , _base_gfs(gfs_storage)
        {
          this->inheritDataSetType(childGridFunctionSpace());
        }

        // We can mask out the following constructors if we don't have template aliases,
        // as we perform the necessary reference <-> shared_ptr conversions in the derived
        // interface class.

        //! Construct a GridFunctionSubSpace from a root space.
        explicit GridFunctionSubSpace(const GFS& gfs)
          : NodeT(TypeTree::childStorage(gfs,TreePath()))
          , FeatureT(gfs)
          , _base_gfs(stackobject_to_shared_ptr(gfs))
        {
          this->inheritDataSetType(childGridFunctionSpace());
        }

        //! Construct a GridFunctionSubSpace from the storage of another GridFunctionSubSpace.
        /**
         * This constructor is used to implement the non-nesting behavior by extracting the
         * original root space from the GridFunctionSubSpace and using that space for initialization.
         * In order to work correctly, this relies on a support wrapper that correctly sets up the
         * TreePath for the new space.
         *
         * \note The second parameter is a little SFINAE helper that removes the constructor from
         *       the overload set if the two spaces are identical to avoid masking the standard
         *       copy constructor.
         */
        template<typename TP>
        explicit GridFunctionSubSpace(std::shared_ptr<const GridFunctionSubSpace<GFS,TP> > gfs_storage, typename std::enable_if<!std::is_same<TP,TreePath>::value,void*>::type = nullptr)
          : NodeT(TypeTree::childStorage(gfs_storage->baseGridFunctionSpace(),TreePath()))
          , FeatureT(gfs_storage->baseGridFunctionSpace())
          , _base_gfs(gfs_storage->baseGridFunctionSpaceStorage())
        {
          setDataSetType(childGridFunctionSpace().dataSetType());
        }

        //! Construct a GridFunctionSubSpace from another GridFunctionSubSpace.
        /**
         * This constructor is used to implement the non-nesting behavior by extracting the
         * original root space from the GridFunctionSubSpace and using that space for initialization.
         * In order to work correctly, this relies on a support wrapper that correctly sets up the
         * TreePath for the new space.
         *
         * \note The second parameter is a little SFINAE helper that removes the constructor from
         *       the overload set if the two spaces are identical to avoid masking the standard
         *       copy constructor.
         */
        template<typename TP>
        explicit GridFunctionSubSpace(const GridFunctionSubSpace<GFS,TP>& gfs, typename std::enable_if<!std::is_same<TP,TreePath>::value,void*>::type = nullptr)
          : NodeT(TypeTree::childStorage(gfs.baseGridFunctionSpace(),TreePath()))
          , FeatureT(gfs.baseGridFunctionSpace())
          , _base_gfs(gfs.baseGridFunctionSpaceStorage())
        {
          this->inheritDataSetType(childGridFunctionSpace());
        }

      public:

        //! The base GridFunctionSpace that this GridFunctionSubSpace is based on.
        typedef GFS BaseGridFunctionSpace;

        //! The type of the original GridFunctionSpace that is the root of this GridFunctionSpace.
        using ChildGridFunctionSpace = TypeTree::ChildForTreePath<GFS,TreePath>;

        //! Re-exported Traits from the original GridFunctionSpace.
        typedef typename ChildGridFunctionSpace::Traits Traits;

        //! Our ImplementationTag is derived from the tag of the original GridFunctionSpace.
        typedef GridFunctionSubSpaceTag<
          TypeTree::ImplementationTag<ChildGridFunctionSpace>
          > ImplementationTag;

        //! Returns the root GridFunctionSpace that this subspace view is based on.
        const BaseGridFunctionSpace& baseGridFunctionSpace() const
        {
          return *_base_gfs;
        }

        //! Returns the storage object of the root GridFunctionSpace that this subspace view is based on.
        std::shared_ptr<const BaseGridFunctionSpace> baseGridFunctionSpaceStorage() const
        {
          return _base_gfs;
        }

        //! Returns the original GridFunctionSpace that we provide a view for.
        /**
         * \warning Users should think *at least* twice before using this object in their
         *          code, as it will usually not do what they want! Due to the way GridFunctionSpaces
         *          are constructed, it is not aware of the overall structure of the space!
         */
        const ChildGridFunctionSpace& childGridFunctionSpace() const
        {
          return this->proxiedNode();
        }

        //! Returns the storage object of the original GridFunctionSpace that we provide a view for.
        /**
         * \warning Users should think *at least* twice before using this object in their
         *          code, as it will usually not do what they want! Due to the way GridFunctionSpaces
         *          are constructed, it is not aware of the overall structure of the space!
         */
        std::shared_ptr<const ChildGridFunctionSpace> childGridFunctionSpaceStorage() const
        {
          return this->proxiedNodeStorage();
        }

        std::string name() const
        {
          return childGridFunctionSpace().name();
        }

        void name(const std::string& name)
        {
          childGridFunctionSpace().name(name);
        }

      private:

        std::shared_ptr<const GFS> _base_gfs;

      };

#ifndef DOXYGEN


      //! Helper TMP to construct non-nested GridFunctionSubSpaces - default case.
      template<typename GFS, typename TreePath>
      struct construct_sub_space
      {
        typedef GridFunctionSubSpace<
          GFS,
          TreePath
          > type;
      };

      //! Helper TMP to construct non-nested GridFunctionSubSpaces - nested case.
      template<typename BaseGFS, typename SubGFSTreePath, typename TreePath>
      struct construct_sub_space<Dune::PDELab::gfs::GridFunctionSubSpace<
                                   BaseGFS,
                                   SubGFSTreePath
                                   >,
                                 TreePath
                                 >
      {
        typedef GridFunctionSubSpace<
          BaseGFS,
          typename TypeTree::TreePathConcat<
            SubGFSTreePath,
            TreePath
            >::type
          > type;
      };

#endif // DOXYGEN

    } // namespace gfs


#if DOXYGEN

    //! \copydoc Dune::PDELab::gfs::GridFunctionSubSpace
    template<typename GFS, typename TreePath>
    using GridFunctionSubSpace = gfs::GridFunctionSubSpace<GFS,TreePath>;

#else // DOXYGEN

    //! \copydoc Dune::PDELab::gfs::GridFunctionSubSpace
    template<typename GFS, typename TreePath>
    using GridFunctionSubSpace = typename gfs::construct_sub_space<GFS,TreePath>::type;

#endif // DOXYGEN

    //! \}

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_SUBSPACE_HH
