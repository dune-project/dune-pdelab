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

    struct EntityBlockedOrderingTag {};

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

      template<typename, typename, typename>
      friend struct power_gfs_to_ordering_descriptor;

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
        r.template child<0>()._container_blocked = false;
        return std::move(r);
      }

      static transformed_storage_type transform_storage(shared_ptr<const GFS> gfs, const Transformation& t)
      {
        transformed_storage_type r(make_shared<transformed_type>(make_tuple(LocalOrderingTransformation::transform_storage(gfs,gfs_to_local_ordering<Transformation>())),gfs->backend().blocked()));
        r->template child<0>()._container_blocked = false;
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
        r.template child<0>()._container_blocked = false;
        return std::move(r);
      }

      static transformed_storage_type transform_storage(shared_ptr<const GFS> gfs, const Transformation& t)
      {
        transformed_storage_type r(make_shared<transformed_type>(make_tuple(LocalOrderingTransformation::transform_storage(gfs,gfs_to_local_ordering<Transformation>())),gfs->backend().blocked()));
        r->template child<0>()._container_blocked = false;
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



#if 0
    //! Interface for merging index spaces
    template<class SizeType, class Tag, class Child, std::size_t k>
    class PowerBlockwiseOrdering :
      public TypeTree::PowerNode<Child, k>,
      public BlockwiseOrderingImp::Base<
        SizeType, Tag, PowerBlockwiseOrdering<SizeType, Tag, Child, k>
      >
    {
      typedef TypeTree::PowerNode<Child, k> Node;

      friend class BlockwiseOrderingImp::Base<
        SizeType, Tag, PowerBlockwiseOrdering
        >;

      // make sure the ordering tag wasn't given more than one argument
      dune_static_assert
      ((IsBaseOf<ComponentBlockwiseOrderingTag<Tag::s0>, Tag>::value),
       "At most one blocksize may be given to a ComponentBlockwiseOrderingTag "
       "to be used in a PowerNode.");

      // check for compatible sizes
      void sizeCheck() {
        static const std::size_t &blockSize =
          BlockwiseOrderingImp::Size<typename Tag::NormalizedTag>::value[0];

        SizeType maxLocalSize = this->child(0).maxLocalSize();
        if(maxLocalSize%blockSize!=0)
          DUNE_THROW(InvalidStateException,
                     className<PowerBlockwiseOrdering>() << ": Number of DOFs "
                     "(" << maxLocalSize << ") per component must be a "
                     "multiple of the BlockSize (" << blockSize << ")");
        for(std::size_t childIndex = 1; childIndex < Node::CHILDREN;
            ++childIndex)
          if(this->child(childIndex).maxLocalSize() != maxLocalSize)
            DUNE_THROW(InvalidStateException,
                       className<PowerBlockwiseOrdering>() << ": Components "
                       "must be of equal size");
      }

    public:
      //! Construct ordering object
      /**
       * In general, an ordering object is not properly setup after
       * construction.  This must be done by a seperate call to update() after
       * all the children have been properly set up.
       *
       * \note This constructor must be present for ordering objects not at
       *       the leaf of the tree.
       */
      template<class GFS>
      PowerBlockwiseOrdering(const GFS &gfs,
                             const typename Node::NodeStorage &children) :
        Node(children)
      { }

      //! map a global dof index from a child
      /**
       * Given the index of a dof in the global dof-vector of one of the
       * children, compute the index of the same dof in the global dof-vector
       * of this ordering.
       *
       * \note update() must have been called before this method may be used.
       */
      SizeType subMap(SizeType child, SizeType indexInChild) const {
        static const SizeType &size = BlockwiseOrderingImp::Size<typename Tag::NormalizedTag>::value[0];

        return indexInChild % size
          + size*child
          + (indexInChild/size) * size*Node::CHILDREN;
      }

      std::string name() const { return "PowerBlockwiseOrdering"; }
    };

    template<std::size_t s0, std::size_t s1, std::size_t s2, std::size_t s3,
             std::size_t s4, std::size_t s5, std::size_t s6, std::size_t s7,
             std::size_t s8, std::size_t s9>
    struct TransformPowerGFSToOrdering<
      ComponentBlockwiseOrderingTag<s0, s1, s2, s3, s4, s5, s6, s7, s8, s9>
      >
    {
      template<class GFSTraits, class TransformedChild, std::size_t k>
      struct result {
        typedef PowerBlockwiseOrdering<
          typename GFSTraits::SizeType,
          ComponentBlockwiseOrderingTag<s0, s1, s2, s3, s4, s5, s6, s7, s8,
                                        s9>,
          TransformedChild,
          k> type;
      };
    };

    template<std::size_t s0, std::size_t s1, std::size_t s2, std::size_t s3,
             std::size_t s4, std::size_t s5, std::size_t s6, std::size_t s7,
             std::size_t s8, std::size_t s9>
    struct TransformPowerGFSToOrdering<
      GridFunctionSpaceComponentBlockwiseMapper<s0, s1, s2, s3, s4, s5, s6, s7, s8, s9>
      > : public TransformPowerGFSToOrdering<
      ComponentBlockwiseOrderingTag<s0, s1, s2, s3, s4, s5, s6, s7, s8, s9>
      >
    {};


    //! Interface for merging index spaces
    template<class SizeType, class Tag,
             DUNE_TYPETREE_COMPOSITENODE_TEMPLATE_CHILDREN>
    class CompositeBlockwiseOrdering :
      public DUNE_TYPETREE_COMPOSITENODE_BASETYPE,
      public BlockwiseOrderingImp::Base<
        SizeType, Tag,
        CompositeBlockwiseOrdering<
          SizeType, Tag, DUNE_TYPETREE_COMPOSITENODE_CHILDTYPES
          >
      >
    {
      typedef DUNE_TYPETREE_COMPOSITENODE_BASETYPE Node;

      friend class BlockwiseOrderingImp::Base<
        SizeType, Tag, CompositeBlockwiseOrdering
        >;

      //////////////////////////////////////////////////////////////////////
      class SizeCheckVisitor :
        public TypeTree::DirectChildrenVisitor,
        public TypeTree::DynamicTraversal
      {
        std::size_t sizeRatio;

      public:
        template<class T, class Child, class TreePath, class ChildIndex>
        void beforeChild(const T &t, const Child& child, TreePath,
                         ChildIndex childIndex)
        {
          static const std::size_t *blockSize =
            BlockwiseOrderingImp::Size<typename Tag::NormalizedTag>::value;

          if(child.maxLocalSize()%blockSize[childIndex]!=0)
            DUNE_THROW(InvalidStateException,
                       className<CompositeBlockwiseOrdering>() << ": Number "
                       "of DOFs (" << child.maxLocalSize() << ") per "
                       "component must be a multiple of the BlockSize "
                       "(" << blockSize[childIndex] << ")");
          if(childIndex == 0)
            sizeRatio = child.maxLocalSize()/blockSize[0];
          else
            if(child.maxLocalSize()/blockSize[childIndex] != sizeRatio)
              DUNE_THROW(InvalidStateException,
                         className<CompositeBlockwiseOrdering>() << ": "
                         "Components must be of equal size");
        }
      };
      //////////////////////////////////////////////////////////////////////

      // check for compatible sizes
      void sizeCheck() {
        SizeCheckVisitor visitor;
        TypeTree::applyToTree(*this, visitor);
      }

    public:
      //! Construct ordering object
      /**
       * In general, an ordering object is not properly setup after
       * construction.  This must be done by a seperate call to update() after
       * all the children have been properly set up.
       *
       * \note This constructor must be present for ordering objects not at
       *       the leaf of the tree.
       */
      template<class GFS>
      CompositeBlockwiseOrdering
      ( const GFS &gfs,
        DUNE_TYPETREE_COMPOSITENODE_STORAGE_CONSTRUCTOR_SIGNATURE) :
        Node(DUNE_TYPETREE_COMPOSITENODE_CHILDVARIABLES)
      { }

      //! map a global dof index from a child
      /**
       * Given the index of a dof in the global dof-vector of one of the
       * children, compute the index of the same dof in the global dof-vector
       * of this ordering.
       *
       * \note update() must have been called before this method may be used.
       */
      SizeType subMap(SizeType child, SizeType indexInChild) const {
        static const SizeType *size = BlockwiseOrderingImp::Size<typename Tag::NormalizedTag>::value;
        static const SizeType *offset =
          BlockwiseOrderingImp::Offset<typename Tag::NormalizedTag>::value;

        return indexInChild % size[child]
          + offset[child]
          + (indexInChild/size[child]) * offset[Node::CHILDREN];
      }

      std::string name() const { return "CompositeBlockwiseOrdering"; }
    };


#if HAVE_VARIADIC_TEMPLATES

    //! Node transformation descriptor for CompositeGridFunctionSpace -> BlockwiseOrdering (with variadic templates).
    template<typename GFSNode, typename Transformation>
    struct VariadicCompositeGFSToBlockwiseOrderingTransformation
    {

      static const bool recursive = true;

      template<typename... TC>
      struct result
      {
        typedef CompositeBlockwiseOrdering<typename Transformation::GridFunctionSpace::Traits::SizeType,
                                           typename Transformation::Ordering,
                                           TC...> type;
        typedef shared_ptr<type> storage_type;
      };

      template<typename... TC>
      static typename result<TC...>::type transform(const GFSNode& s, const Transformation& t, shared_ptr<TC>... children)
      {
        return typename result<TC...>::type(t.asGridFunctionSpace(s),children...);
      }

      template<typename... TC>
      static typename result<TC...>::storage_type transform_storage(shared_ptr<const GFSNode> s, const Transformation& t, shared_ptr<TC>... children)
      {
        return make_shared<typename result<TC...>::type>(t.asGridFunctionSpace(*s),children...);
      }

    };

    // Register transformation.
    template<typename GFSNode, typename GFS,
             std::size_t s0, std::size_t s1, std::size_t s2, std::size_t s3,
             std::size_t s4, std::size_t s5, std::size_t s6, std::size_t s7,
             std::size_t s8, std::size_t s9>
    VariadicCompositeGFSToBlockwiseOrderingTransformation<GFSNode,gfs_to_ordering<GFS,ComponentBlockwiseOrderingTag<s0, s1, s2, s3, s4, s5, s6, s7,s8, s9> > >
    lookupNodeTransformation(GFSNode*, gfs_to_ordering<GFS,ComponentBlockwiseOrderingTag<s0, s1, s2, s3, s4, s5, s6, s7,s8, s9> >*, CompositeGridFunctionSpaceBaseTag);

    // register old GridFunctionSpaceComponentBlockwiseMapper separately - unfortunately, this is necessary to make the specialization match correctly
    template<typename GFSNode, typename GFS,
             std::size_t s0, std::size_t s1, std::size_t s2, std::size_t s3,
             std::size_t s4, std::size_t s5, std::size_t s6, std::size_t s7,
             std::size_t s8, std::size_t s9>
    VariadicCompositeGFSToBlockwiseOrderingTransformation<GFSNode,gfs_to_ordering<GFS,GridFunctionSpaceComponentBlockwiseMapper<s0, s1, s2, s3, s4, s5, s6, s7,s8, s9> > >
    lookupNodeTransformation(GFSNode*, gfs_to_ordering<GFS,GridFunctionSpaceComponentBlockwiseMapper<s0, s1, s2, s3, s4, s5, s6, s7,s8, s9> >*, CompositeGridFunctionSpaceBaseTag);

#else // HAVE_VARIADIC_TEMPLATES

    //! Node transformation descriptor for CompositeGridFunctionSpace -> BlockwiseOrdering (without variadic templates).
    template<typename GFSNode, typename Transformation>
    struct CompositeGFSToBlockwiseOrderingTransformation
    {

      static const bool recursive = true;

      template<typename TC0,
               typename TC1,
               typename TC2,
               typename TC3,
               typename TC4,
               typename TC5,
               typename TC6,
               typename TC7,
               typename TC8,
               typename TC9>
      struct result
      {
        typedef CompositeBlockwiseOrdering<typename Transformation::GridFunctionSpace::Traits::SizeType,
                                           typename Transformation::Ordering,
                                           TC0,TC1,TC2,TC3,TC4,TC5,TC6,TC7,TC8,TC9> type;
        typedef shared_ptr<type> storage_type;
      };

      template<typename TC0,
               typename TC1,
               typename TC2,
               typename TC3,
               typename TC4,
               typename TC5,
               typename TC6,
               typename TC7,
               typename TC8,
               typename TC9>
      static typename result<TC0,TC1,TC2,TC3,TC4,TC5,TC6,TC7,TC8,TC9>::type
      transform(const GFSNode& s,
                const Transformation& t,
                shared_ptr<TC0> c0,
                shared_ptr<TC1> c1,
                shared_ptr<TC2> c2,
                shared_ptr<TC3> c3,
                shared_ptr<TC4> c4,
                shared_ptr<TC5> c5,
                shared_ptr<TC6> c6,
                shared_ptr<TC7> c7,
                shared_ptr<TC8> c8,
                shared_ptr<TC9> c9)
      {
        return typename result<TC0,TC1,TC2,TC3,TC4,TC5,TC6,TC7,TC8,TC9>::type(t.asGridFunctionSpace(s),c0,c1,c2,c3,c4,c5,c6,c7,c8,c9);
      }

      template<typename TC0,
               typename TC1,
               typename TC2,
               typename TC3,
               typename TC4,
               typename TC5,
               typename TC6,
               typename TC7,
               typename TC8,
               typename TC9>
      static typename result<TC0,TC1,TC2,TC3,TC4,TC5,TC6,TC7,TC8,TC9>::storage_type
      transform_storage(shared_ptr<const GFSNode> s,
                        const Transformation& t,
                        shared_ptr<TC0> c0,
                        shared_ptr<TC1> c1,
                        shared_ptr<TC2> c2,
                        shared_ptr<TC3> c3,
                        shared_ptr<TC4> c4,
                        shared_ptr<TC5> c5,
                        shared_ptr<TC6> c6,
                        shared_ptr<TC7> c7,
                        shared_ptr<TC8> c8,
                        shared_ptr<TC9> c9)
      {
        return make_shared<typename result<TC0,TC1,TC2,TC3,TC4,TC5,TC6,TC7,TC8,TC9>::type>(t.asGridFunctionSpace(s),c0,c1,c2,c3,c4,c5,c6,c7,c8,c9);
      }

    };


    // Register transformation.
    template<typename GFSNode, typename GFS,
             std::size_t s0, std::size_t s1, std::size_t s2, std::size_t s3,
             std::size_t s4, std::size_t s5, std::size_t s6, std::size_t s7,
             std::size_t s8, std::size_t s9>
    CompositeGFSToBlockwiseOrderingTransformation<GFSNode,gfs_to_ordering<GFS,GridFunctionSpaceComponentBlockwiseMapper<s0, s1, s2, s3, s4, s5, s6, s7,s8, s9> > >
    lookupNodeTransformation(GFSNode*, gfs_to_ordering<GFS,ComponentBlockwiseOrderingTag<s0, s1, s2, s3, s4, s5, s6, s7,s8, s9> >*, CompositeGridFunctionSpaceBaseTag);

    // register old GridFunctionSpaceComponentBlockwiseMapper separately - unfortunately, this is necessary to make the specialization match correctly
    template<typename GFSNode, typename GFS,
             std::size_t s0, std::size_t s1, std::size_t s2, std::size_t s3,
             std::size_t s4, std::size_t s5, std::size_t s6, std::size_t s7,
             std::size_t s8, std::size_t s9>
    CompositeGFSToBlockwiseOrderingTransformation<GFSNode,gfs_to_ordering<GFS,GridFunctionSpaceComponentBlockwiseMapper<s0, s1, s2, s3, s4, s5, s6, s7,s8, s9> > >
    lookupNodeTransformation(GFSNode*, gfs_to_ordering<GFS,ComponentBlockwiseOrderingTag<s0, s1, s2, s3, s4, s5, s6, s7,s8, s9> >*, CompositeGridFunctionSpaceBaseTag);

#endif // HAVE_VARIADIC_TEMPLATES

#endif // 0


   //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_ENTITYBLOCKEDORDERING_HH
