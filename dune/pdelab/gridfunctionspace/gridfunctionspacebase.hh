// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_GRIDFUNCTIONSPACEBASE_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_GRIDFUNCTIONSPACEBASE_HH

#include <optional>

#include <dune/typetree/visitor.hh>
#include <dune/typetree/traversal.hh>

#include <dune/pdelab/common/exceptions.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup GridFunctionSpace grid function space
    //! \ingroup PDELab
    //! \{

#ifndef DOXYGEN

    // forward declaration for friend declaration
    template<typename GFS, typename GFSTraits>
    class GridFunctionSpaceBase;

    namespace impl {

      struct reset_root_space_flag
        : public TypeTree::DirectChildrenVisitor
        , public TypeTree::DynamicTraversal
      {

        template<typename GFS, typename Child, typename TreePath, typename ChildIndex>
        void afterChild(const GFS& gfs, Child& child, TreePath, ChildIndex) const
        {
          if (child._initialized && child._is_root_space)
            {
              DUNE_THROW(GridFunctionSpaceHierarchyError,"initialized space cannot become part of larger GridFunctionSpace tree");
            }
          child._is_root_space = false;
        }

      };

      template<typename size_type>
      struct update_ordering_data;

      // helper class with minimal dependencies. Orderings keep a pointer to this structure and populate it
      // during their update procedure.

      template<typename size_type>
      class GridFunctionSpaceOrderingData
      {

        template<typename,typename>
        friend class ::Dune::PDELab::GridFunctionSpaceBase;

        template<typename>
        friend struct update_ordering_data;

        GridFunctionSpaceOrderingData()
          : _size(0)
          , _block_count(0)
          , _global_size(0)
          , _max_local_size(0)
          , _is_root_space(true)
          , _initialized(false)
          , _size_available(true)
        {}

        size_type _size;
        size_type _block_count;
        size_type _global_size;
        size_type _max_local_size;
        bool _is_root_space;
        bool _initialized;
        bool _size_available;

      };

      template<typename size_type>
      struct update_ordering_data
        : public TypeTree::TreeVisitor
        , public TypeTree::DynamicTraversal
      {

        typedef GridFunctionSpaceOrderingData<size_type> Data;

        template<typename Ordering>
        void update(const Ordering& ordering, bool is_root)
        {
          if (ordering._gfs_data)
            {
              Data& data = *ordering._gfs_data;
              // if (data._initialized && data._is_root_space && !is_root)
              //   {
              //     DUNE_THROW(GridFunctionSpaceHierarchyError,"former root space is now part of a larger tree");
              //   }
              data._initialized = true;
              data._global_size = _global_size;
              data._max_local_size = _max_local_size;
              data._size_available = ordering.update_gfs_data_size(data._size,data._block_count);
            }
        }

        template<typename Ordering, typename TreePath>
        void leaf(const Ordering& ordering, TreePath tp)
        {
          update(ordering,tp.size() == 0);
        }

        template<typename Ordering, typename TreePath>
        void post(const Ordering& ordering, TreePath tp)
        {
          update(ordering,tp.size() == 0);
        }

        template<typename Ordering>
        explicit update_ordering_data(const Ordering& ordering)
          : _global_size(ordering.size())
          , _max_local_size(ordering.maxLocalSize())
        {}

        const size_type _global_size;
        const size_type _max_local_size;

      };


      //! Checks that every leaf node has the same entity set
      template<class EntitySet>
      struct common_entity_set
        : public TypeTree::TreeVisitor
        , public TypeTree::DynamicTraversal
      {
        template<typename T, typename TreePath>
        void leaf(T&& t, TreePath treePath) {
          if (not _entity_set)
            _entity_set = t.entitySet();
          else if (*_entity_set != t.entitySet())
            DUNE_THROW(GridFunctionSpaceHierarchyError, "Entity sets should match!");
        }

        std::optional<EntitySet> _entity_set;
      };

      /**
       * @brief  Updates every entity set in leaf nodes
       * @details Potentially, every leaf node has a different entity set.
       *  We only update the first of common entity sets
       */
      template<class EntitySet>
      struct update_leaf_entity_set
        : public TypeTree::TreeVisitor
        , public TypeTree::DynamicTraversal
      {
        update_leaf_entity_set(const std::optional<EntitySet>& entity_set, bool force_update)
          : _entity_set{entity_set}
          , _force_update{force_update}
        {}

        template<typename GFSNode, typename TreePath>
        void leaf(GFSNode&& gfs_node, TreePath treePath) {
          if (not _entity_set)
            _entity_set = gfs_node.entitySet();
          if (*_entity_set != gfs_node.entitySet()) {
            gfs_node.entitySet().update(_force_update);
            _entity_set = gfs_node.entitySet();
          }
        }

        bool _force_update;
        std::optional<EntitySet> _entity_set;
      };

    } // namespace impl

#endif // DOXYGEN


    template<typename GFS, typename GFSTraits>
    class GridFunctionSpaceBase
      : public impl::GridFunctionSpaceOrderingData<typename GFSTraits::SizeType>
    {

      friend struct impl::reset_root_space_flag;

    public:

      typedef GFSTraits Traits;

      template<typename Backend_, typename OrderingTag_>
      GridFunctionSpaceBase(Backend_&& backend, OrderingTag_&& ordering_tag)
        : _backend(std::forward<Backend_>(backend))
        , _ordering_tag(std::forward<OrderingTag_>(ordering_tag))
      {
        TypeTree::applyToTree(gfs(),impl::reset_root_space_flag());
      }

      typename Traits::SizeType size() const
      {
        if (!_initialized)
          {
            DUNE_THROW(UninitializedGridFunctionSpaceError,"space is not initialized");
          }
        if (!_size_available)
          {
            DUNE_THROW(GridFunctionSpaceHierarchyError,
                       "Size cannot be calculated at this point in the GFS tree.");
          }
        return _size;
      }

      typename Traits::SizeType blockCount() const
      {
        if (!_initialized)
          {
            DUNE_THROW(UninitializedGridFunctionSpaceError,"space is not initialized");
          }
        if (!_size_available)
          {
            DUNE_THROW(GridFunctionSpaceHierarchyError,
                       "Block count cannot be calculated at this point in the GFS tree.");
          }
        return _block_count;
      }

      typename Traits::SizeType globalSize() const
      {
        if (!_initialized)
          {
            DUNE_THROW(UninitializedGridFunctionSpaceError,"space is not initialized");
          }
        return _global_size;
      }

      //! get max dimension of shape function space
      typename Traits::SizeType maxLocalSize () const
      {
        if (!_initialized)
          {
            DUNE_THROW(UninitializedGridFunctionSpaceError,"space is not initialized");
          }
        return _max_local_size;
      }

      //! Update the indexing information of the GridFunctionSpace.
      /**
       *
       * \ param force   Set to true if the underlying grid has changed (e.g. due to adaptivity)
       *                 to force an update of the embedded EntitySet.
       */
      void update(bool force = false)
      {
        _entity_set->update(force);
        auto update_leaf_es = impl::update_leaf_entity_set{_entity_set, force};
        TypeTree::applyToTree(*this, update_leaf_es);
        // We bypass the normal access using ordering() here to avoid a double
        // update if the Ordering has not been created yet.
        if (!gfs()._ordering)
          gfs().create_ordering();
        update(*gfs()._ordering);
      }

      const std::string& name() const
      {
        return _name;
      }

      void name(const std::string& name)
      {
        _name = name;
      }

      typename Traits::Backend& backend()
      {
        return _backend;
      }

      const typename Traits::Backend& backend() const
      {
        return _backend;
      }

      //! get grid view
      const typename Traits::GridView& gridView () const
      {
        return entitySet().gridView();
      }

      //! get EntitySet
      const typename Traits::EntitySet& entitySet () const
      {
        assert(_entity_set && "No entity set has been assigned to this node");
        return *_entity_set;
      }

      typename Traits::EntitySet& entitySet ()
      {
        assert(_entity_set && "No entity set has been assigned to this node");
        return *_entity_set;
      }


      //! get EntitySet
      bool hasEntitySet () const
      {
        return _entity_set.has_value();
      }

      void setEntitySet(const typename Traits::EntitySet& entity_set)
      {
        _entity_set.emplace(entity_set);
      }

      typename Traits::OrderingTag& orderingTag()
      {
        return _ordering_tag;
      }

      const typename Traits::OrderingTag& orderingTag() const
      {
        return _ordering_tag;
      }

      bool isRootSpace() const
      {
        return _is_root_space;
      }

    protected:

      template<typename Ordering>
      void update(Ordering& ordering) const
      {
        if (!_is_root_space)
          {
            DUNE_THROW(GridFunctionSpaceHierarchyError,"update() may only be called on the root of the function space hierarchy");
          }
        ordering.update();
        TypeTree::applyToTree(ordering,impl::update_ordering_data<typename Traits::SizeType>(ordering));
      }

      mutable std::optional<typename Traits::EntitySet> _entity_set;

    private:

      typedef impl::GridFunctionSpaceOrderingData<typename GFSTraits::SizeType> BaseT;

      GFS& gfs()
      {
        return static_cast<GFS&>(*this);
      }

      const GFS& gfs() const
      {
        return static_cast<const GFS&>(*this);
      }

      std::string _name;
      typename Traits::Backend _backend;
      typename Traits::OrderingTag _ordering_tag;

      using BaseT::_size;
      using BaseT::_block_count;
      using BaseT::_global_size;
      using BaseT::_max_local_size;
      using BaseT::_is_root_space;
      using BaseT::_initialized;
      using BaseT::_size_available;

    };


  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_GRIDFUNCTIONSPACEBASE_HH
