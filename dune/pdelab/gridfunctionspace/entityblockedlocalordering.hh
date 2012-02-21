// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_ENTITYBLOCKEDLOCALORDERING_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_ENTITYBLOCKEDLOCALORDERING_HH

#include <cstddef>
#include <ostream>
#include <string>

#include <dune/common/classname.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/stdstreams.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/deprecated.hh>

#include <dune/pdelab/common/typetree/compositenodemacros.hh>
#include <dune/pdelab/common/typetree/powernode.hh>
#include <dune/pdelab/common/typetree/traversal.hh>
#include <dune/pdelab/common/typetree/visitor.hh>

#include <dune/pdelab/gridfunctionspace/gridviewordering.hh>

namespace Dune {
  namespace PDELab {

    //! Interface for merging index spaces
    template<typename ChildOrdering, std::size_t k>
    class PowerEntityBlockedLocalOrdering
      : public TypeTree::PowerNode<ChildOrdering,k>
      , public LocalOrderingBase<typename ChildOrdering::Traits::GridView,
                                 typename ChildOrdering::Traits::DOFIndex,
                                 typename ChildOrdering::Traits::ContainerIndex>
    {

      typedef TypeTree::PowerNode<ChildOrdering,k> NodeT;
      typedef LocalOrderingBase<typename ChildOrdering::Traits::GridView,
                                typename ChildOrdering::Traits::DOFIndex,
                                typename ChildOrdering::Traits::ContainerIndex> BaseT;

    public:

      static const bool consume_tree_index = true;

      typedef typename BaseT::Traits Traits;

      PowerEntityBlockedLocalOrdering(const typename NodeT::NodeStorage& child_storage, bool container_blocked)
        : NodeT(child_storage)
        , BaseT(*this,container_blocked)
      {}

      const typename Traits::GridView& gridView() const
      {
        return this->child(0).gridView();
      }

      //! update internal data structures
      /**
       * This method must be called after initialization and every time the
       * structure of the dof-vector of one of gfs's children changes.  All
       * the children must have been set up properly before the call to
       * update().
       */
      void update() {
        /*
        Dune::dinfo << asImp().name() << ":" << std::endl;

        if(!NonLeafOrderingBase<SizeType, Imp>::fixedSize())
          DUNE_THROW(InvalidStateException, className<Imp>() << " works "
                     "only with children that have a uniform size for all "
                     "entities of a given geometry type/all intersections");

        asImp().sizeCheck();

        printInfo(dinfo);
        */
      }

#if 0
      //! \brief offset of the block of dofs attached to a given entity (of
      //!        arbitrary codimension)
      /**
       * \note This method should be available for all types of entities
       *       required by the grid.  It is mostly required for
       *       communcation, so if it is known that this method is not
       *       actually called for a given entity type the implementation
       *       may throw NotImplemented.
       * \note If the grid does not support a given entity type, it may
       *       still be possible to get this information using
       *       entityOffset(const Element &e, std::size_t codim, std::size_t
       *       subentity).
       *
       * \throw NotImplemented        If this EntityType is not supported by
       *                              the ordering.
       * \throw InvalidStateException If blocked()==false.
       */
      template<class Entity>
      SizeType entityOffset(const Entity &e) const {
        return asImp().subMap(0,
                              asImp().template child<0>().entityOffset(e));
      }
      //! \brief offset of the blocks of dofs attached to a given subentity
      //!        of an element
      /**
       * This method determines the starting offset of the block of dofs
       * attached to a subentity of the given codim 0 entity.  If the grid
       * (and the ordering) directly support entities of the given
       * codimension, this is equivalent to calling
       * entityOffset(*e.subEntity<codim>(subentity)).
       */
      template<class Element>
      SizeType entityOffset(const Element &e, std::size_t codim,
                            std::size_t subentity) const {
        return asImp().subMap
          (0, asImp().template child<0>().entityOffset(e, codim, subentity));
      }
      //! offset of the block of dofs attached to a given intersection
      template<class Intersection>
      SizeType intersectionOffset(const Intersection &i) const {
        return asImp().subMap
          (0, asImp().template child<0>().intersectionOffset(i));
      }

#endif // 0

    private:

      using BaseT::_container_blocked;

    };

    template<typename GFS, typename Transformation, typename OrderingTag>
    struct power_gfs_to_local_ordering_descriptor;


    template<typename GFS, typename Transformation>
    struct power_gfs_to_local_ordering_descriptor<GFS,Transformation,EntityBlockedOrderingTag>
    {

      static const bool recursive = true;

      template<typename TC>
      struct result
      {
        typedef PowerEntityBlockedLocalOrdering<TC,GFS::CHILDREN> type;
        typedef shared_ptr<type> storage_type;
      };

      template<typename TC>
      static typename result<TC>::type transform(const GFS& gfs, const Transformation& t, const array<shared_ptr<TC>,GFS::CHILDREN>& children)
      {
        return typename result<TC>::type(children,gfs.backend().blocked());
      }

      template<typename TC>
      static typename result<TC>::storage_type transform_storage(shared_ptr<const GFS> gfs, const Transformation& t, const array<shared_ptr<TC>,GFS::CHILDREN>& children)
      {
        return make_shared<typename result<TC>::type>(children,gfs->backend().blocked());
      }

    };

    template<typename GFS, typename Params>
    power_gfs_to_local_ordering_descriptor<
      GFS,
      gfs_to_local_ordering<Params>,
      typename GFS::OrderingTag
      >
    lookupNodeTransformation(GFS* gfs, gfs_to_local_ordering<Params>* t, PowerGridFunctionSpaceTag tag);


    template<typename GFS, typename Transformation>
    struct power_gfs_to_ordering_descriptor<GFS,Transformation,EntityBlockedOrderingTag>
    {

      static const bool recursive = false;

      typedef TypeTree::TransformTree<GFS,gfs_to_local_ordering<Transformation> > LocalOrderingTransformation;
      typedef typename LocalOrderingTransformation::Type LocalOrdering;

      typedef GridViewOrdering<LocalOrdering> transformed_type;

      typedef shared_ptr<transformed_type> transformed_storage_type;

      static transformed_type transform(const GFS& gfs, const Transformation& t)
      {
        transformed_type r(transformed_type(make_tuple(make_shared<LocalOrdering>(LocalOrderingTransformation::transform(gfs,gfs_to_local_ordering<Transformation>()))),gfs.backend().blocked()));
        // r.template child<0>()._container_blocked = false;
        return std::move(r);
      }

      static transformed_storage_type transform_storage(shared_ptr<const GFS> gfs, const Transformation& t)
      {
        transformed_storage_type r(make_shared<transformed_type>(make_tuple(LocalOrderingTransformation::transform_storage(gfs,gfs_to_local_ordering<Transformation>())),gfs->backend().blocked()));
        // r->template child<0>()._container_blocked = false;
        return std::move(r);
      }

    };

    template<typename GFS, typename Params>
    power_gfs_to_ordering_descriptor<
      GFS,
      gfs_to_ordering<Params>,
      typename GFS::OrderingTag
      >
    lookupNodeTransformation(GFS* gfs, gfs_to_ordering<Params>* t, PowerGridFunctionSpaceTag tag);










    //! Interface for merging index spaces
    template<DUNE_TYPETREE_COMPOSITENODE_TEMPLATE_CHILDREN>
    class CompositeEntityBlockedLocalOrdering
      : public DUNE_TYPETREE_COMPOSITENODE_BASETYPE
      , public LocalOrderingBase<typename DUNE_TYPETREE_COMPOSITENODE_FIRST_CHILD::Traits::GridView,
                                 typename DUNE_TYPETREE_COMPOSITENODE_FIRST_CHILD::Traits::DOFIndex,
                                 typename DUNE_TYPETREE_COMPOSITENODE_FIRST_CHILD::Traits::ContainerIndex>
    {

      typedef DUNE_TYPETREE_COMPOSITENODE_BASETYPE Node;
      typedef LocalOrderingBase<typename DUNE_TYPETREE_COMPOSITENODE_FIRST_CHILD::Traits::GridView,
                                typename DUNE_TYPETREE_COMPOSITENODE_FIRST_CHILD::Traits::DOFIndex,
                                typename DUNE_TYPETREE_COMPOSITENODE_FIRST_CHILD::Traits::ContainerIndex> Base;

    public:

      typedef typename Base::Traits Traits;

      static const bool consume_tree_index = true;

      CompositeEntityBlockedLocalOrdering(bool container_blocked, DUNE_TYPETREE_COMPOSITENODE_STORAGE_CONSTRUCTOR_SIGNATURE)
        : Node(DUNE_TYPETREE_COMPOSITENODE_CHILDVARIABLES)
        , Base(*this,container_blocked)
      {}

      const typename Traits::GridView& gridView() const
      {
        return this->template child<0>().gridView();
      }

      //! update internal data structures
      /**
       * This method must be called after initialization and every time the
       * structure of the dof-vector of one of gfs's children changes.  All
       * the children must have been set up properly before the call to
       * update().
       */
      void update() {
        /*
        Dune::dinfo << asImp().name() << ":" << std::endl;

        if(!NonLeafOrderingBase<SizeType, Imp>::fixedSize())
          DUNE_THROW(InvalidStateException, className<Imp>() << " works "
                     "only with children that have a uniform size for all "
                     "entities of a given geometry type/all intersections");

        asImp().sizeCheck();

        printInfo(dinfo);
        */
      }

#if 0
      //! \brief offset of the block of dofs attached to a given entity (of
      //!        arbitrary codimension)
      /**
       * \note This method should be available for all types of entities
       *       required by the grid.  It is mostly required for
       *       communcation, so if it is known that this method is not
       *       actually called for a given entity type the implementation
       *       may throw NotImplemented.
       * \note If the grid does not support a given entity type, it may
       *       still be possible to get this information using
       *       entityOffset(const Element &e, std::size_t codim, std::size_t
       *       subentity).
       *
       * \throw NotImplemented        If this EntityType is not supported by
       *                              the ordering.
       * \throw InvalidStateException If blocked()==false.
       */
      template<class Entity>
      SizeType entityOffset(const Entity &e) const {
        return asImp().subMap(0,
                              asImp().template child<0>().entityOffset(e));
      }
      //! \brief offset of the blocks of dofs attached to a given subentity
      //!        of an element
      /**
       * This method determines the starting offset of the block of dofs
       * attached to a subentity of the given codim 0 entity.  If the grid
       * (and the ordering) directly support entities of the given
       * codimension, this is equivalent to calling
       * entityOffset(*e.subEntity<codim>(subentity)).
       */
      template<class Element>
      SizeType entityOffset(const Element &e, std::size_t codim,
                            std::size_t subentity) const {
        return asImp().subMap
          (0, asImp().template child<0>().entityOffset(e, codim, subentity));
      }
      //! offset of the block of dofs attached to a given intersection
      template<class Intersection>
      SizeType intersectionOffset(const Intersection &i) const {
        return asImp().subMap
          (0, asImp().template child<0>().intersectionOffset(i));
      }

#endif // 0

    };

    template<typename GFS, typename Transformation, typename OrderingTag>
    struct composite_gfs_to_local_ordering_descriptor;


    template<typename GFS, typename Transformation>
    struct composite_gfs_to_local_ordering_descriptor<GFS,Transformation,EntityBlockedOrderingTag>
    {

      static const bool recursive = true;

      template<typename... TC>
      struct result
      {
        typedef CompositeEntityBlockedLocalOrdering<TC...> type;
        typedef shared_ptr<type> storage_type;
      };

      template<typename... TC>
      static typename result<TC...>::type transform(const GFS& gfs, const Transformation& t, shared_ptr<TC>... children)
      {
        return typename result<TC...>::type(gfs.backend().blocked(),children...);
      }

      template<typename... TC>
      static typename result<TC...>::storage_type transform_storage(shared_ptr<const GFS> gfs, const Transformation& t, shared_ptr<TC>... children)
      {
        return make_shared<typename result<TC...>::type>(gfs.backend().blocked(),children...);
      }

    };

    template<typename GFS, typename Params>
    composite_gfs_to_local_ordering_descriptor<
      GFS,
      gfs_to_local_ordering<Params>,
      typename GFS::OrderingTag
      >
    lookupNodeTransformation(GFS* gfs, gfs_to_local_ordering<Params>* t, CompositeGridFunctionSpaceTag tag);


    template<typename GFS, typename Transformation>
    struct composite_gfs_to_ordering_descriptor<GFS,Transformation,EntityBlockedOrderingTag>
    {

      static const bool recursive = false;

      typedef TypeTree::TransformTree<GFS,gfs_to_local_ordering<Transformation> > LocalOrderingTransformation;
      typedef typename LocalOrderingTransformation::Type LocalOrdering;

      typedef GridViewOrdering<LocalOrdering> transformed_type;

      typedef shared_ptr<transformed_type> transformed_storage_type;

      static transformed_type transform(const GFS& gfs, const Transformation& t)
      {
        transformed_type r(make_tuple(make_shared<LocalOrdering>(LocalOrderingTransformation::transform(gfs,gfs_to_local_ordering<Transformation>()))),gfs.backend().blocked());
        // r.template child<0>()._container_blocked = false;
        return std::move(r);
      }

      static transformed_storage_type transform_storage(shared_ptr<const GFS> gfs, const Transformation& t)
      {
        transformed_storage_type r(make_shared<transformed_type>(make_tuple(LocalOrderingTransformation::transform_storage(gfs,gfs_to_local_ordering<Transformation>())),gfs->backend().blocked()));
        // r->template child<0>()._container_blocked = false;
        return std::move(r);
      }

    };

    template<typename GFS, typename Params>
    composite_gfs_to_ordering_descriptor<
      GFS,
      gfs_to_ordering<Params>,
      typename GFS::OrderingTag
      >
    lookupNodeTransformation(GFS* gfs, gfs_to_ordering<Params>* t, CompositeGridFunctionSpaceTag tag);


   //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_ENTITYBLOCKEDORDERING_HH
