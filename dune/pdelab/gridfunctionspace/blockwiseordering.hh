// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_BLOCKWISEORDERING_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_BLOCKWISEORDERING_HH

#include <cstddef>
#include <ostream>
#include <string>

#include <dune/common/classname.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/stdstreams.hh>

#include <dune/pdelab/common/typetree/compositenodemacros.hh>
#include <dune/pdelab/common/typetree/powernode.hh>
#include <dune/pdelab/common/typetree/traversal.hh>
#include <dune/pdelab/common/typetree/visitor.hh>
#include <dune/pdelab/gridfunctionspace/nonleaforderingbase.hh>
#include <dune/pdelab/gridfunctionspace/orderingbase.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup GridFunctionSpace
    //! \ingroup PDELab
    //! \{

    //! \brief Indicate blockwise ordering of the unknowns of non-leaf grid
    //!        function spaces.
    /**
     * This class instructs the non-leaf GridFunctionSpaces to order the dofs
     * of the child-GridFunctionSpaces in a blockwise manner in the
     * combined dof-vector, i.e. { [s0 dofs of child 0] ... [s(k-1) dofs of
     * child (k-1)] } repeated for each grid element.
     *
     * \note This works only as long as there are no dofs on entities with
     *       codim > 0.
     * \note It is assumed that all grid elements have the same number of
     *       dofs, otherwise use DynamicBlockwiseOrdering.
     */
    template<std::size_t s0 = 1, std::size_t s1 = s0, std::size_t s2 = s1,
             std::size_t s3 = s2, std::size_t s4 = s3, std::size_t s5 = s4,
             std::size_t s6 = s5, std::size_t s7 = s6, std::size_t s8 = s7,
             std::size_t s9 = s8>
    struct ComponentBlockwiseOrderingTag { };
    template<int s0 = 1, int s1 = s0, int s2 = s1, int s3 = s2, int s4 = s3,
             int s5 = s4, int s6 = s5, int s7 = s6, int s8 = s7, int s9 = s8>
    struct GridFunctionSpaceComponentBlockwiseMapper :
      public ComponentBlockwiseOrderingTag<s0, s1, s2, s3, s4, s5, s6, s7, s8,
                                           s9>
    { };

    typedef ComponentBlockwiseOrderingTag<> BlockwiseOrderingTag;
    typedef BlockwiseOrderingTag GridFunctionSpaceBlockwiseMapper;

    namespace BlockwiseOrderingImp {

      //////////////////////////////////////////////////////////////////////
      //
      //  Size
      //

      template<class ComponentBlockwiseOrderingTag>
      struct Size;

      template<std::size_t s0, std::size_t s1, std::size_t s2, std::size_t s3,
               std::size_t s4, std::size_t s5, std::size_t s6, std::size_t s7,
               std::size_t s8, std::size_t s9>
      struct Size<ComponentBlockwiseOrderingTag<s0, s1, s2, s3, s4, s5, s6, s7,
                                                s8, s9> >
      {
        static const std::size_t value[10];
      };
      template<std::size_t s0, std::size_t s1, std::size_t s2, std::size_t s3,
               std::size_t s4, std::size_t s5, std::size_t s6, std::size_t s7,
               std::size_t s8, std::size_t s9>
      const std::size_t Size<
        ComponentBlockwiseOrderingTag<s0, s1, s2, s3, s4, s5, s6, s7, s8, s9>
        >::value[10] = { s0, s1, s2, s3, s4, s5, s6, s7, s8, s9 };

      //////////////////////////////////////////////////////////////////////
      //
      //  Offset
      //


      template<class ComponentBlockwiseOrderingTag>
      struct Offset;

      template<std::size_t s0, std::size_t s1, std::size_t s2, std::size_t s3,
               std::size_t s4, std::size_t s5, std::size_t s6, std::size_t s7,
               std::size_t s8, std::size_t s9>
      struct Offset<ComponentBlockwiseOrderingTag<s0, s1, s2, s3, s4, s5, s6,
                                                  s7, s8, s9> >
      {
        static const std::size_t value[11];
      };
      template<std::size_t s0, std::size_t s1, std::size_t s2, std::size_t s3,
               std::size_t s4, std::size_t s5, std::size_t s6, std::size_t s7,
               std::size_t s8, std::size_t s9>
      const std::size_t Offset<
        ComponentBlockwiseOrderingTag<s0, s1, s2, s3, s4, s5, s6, s7, s8, s9>
        >::value[11] =
          { 0, s0, s0+s1, s0+s1+s2, s0+s1+s2+s3, s0+s1+s2+s3+s4,
            s0+s1+s2+s3+s4+s5, s0+s1+s2+s3+s4+s5+s6, s0+s1+s2+s3+s4+s5+s6+s7,
            s0+s1+s2+s3+s4+s5+s6+s7+s8, s0+s1+s2+s3+s4+s5+s6+s7+s8+s9 };

      //////////////////////////////////////////////////////////////////////

      template<class OrderingImp, class Tag>
      class SizeCheckVisitor :
        public TypeTree::DirectChildrenVisitor,
        public TypeTree::DynamicTraversal
      {
        std::size_t sizeRatio;

      public:
        template<class T, class Child, class TreePath, class ChildIndex>
        void beforeChild(const T &t, const Child& child, TreePath, ChildIndex childIndex) {
          static const std::size_t *blockSize = Size<Tag>::value;

          if(child.maxLocalSize()%blockSize[childIndex]!=0)
            DUNE_THROW(InvalidStateException, className<OrderingImp>() << ": "
                       "Number of DOFs (" << child.maxLocalSize() << ") per "
                       "component must be a multiple of the BlockSize "
                       "(" << blockSize[childIndex] << ")");
          if(childIndex == 0)
            sizeRatio = child.maxLocalSize()/blockSize[0];
          else
            if(child.maxLocalSize()/blockSize[childIndex] != sizeRatio)
              DUNE_THROW(InvalidStateException,
                         className<OrderingImp>() << ": Components must be of "
                         "equal size");
        }
      };

      //////////////////////////////////////////////////////////////////////

      //! Interface for merging index spaces
      template<class SizeType, class Tag, class Imp>
      class Base :
        public NonLeafOrderingBase<SizeType, Imp>
      {
      protected:
        using NonLeafOrderingBase<SizeType, Imp>::asImp;
        using NonLeafOrderingBase<SizeType, Imp>::printInfo;

      public:
        Base() { asImp().update(); }

        //! update internal data structures
        /**
         * This method must be called after initialization and every time the
         * structure of the dof-vector of one of gfs's children changes.  All
         * the children must have been set up properly before the call to
         * update().
         */
        void update() {
          Dune::dinfo << asImp().name() << ":" << std::endl;

          if(!NonLeafOrderingBase<SizeType, Imp>::blocked())
            DUNE_THROW(InvalidStateException, className<Imp>() << " works "
                       "only with blocking children");

          if(!NonLeafOrderingBase<SizeType, Imp>::fixedSize())
            DUNE_THROW(InvalidStateException, className<Imp>() << " works "
                       "only with children that have a uniform size for all "
                       "entities of a given geometry type/all intersections");

          // check for compatible sizes
          {
            SizeCheckVisitor<Imp, Tag> visitor;
            TypeTree::applyToTree(asImp(), visitor);
          }

          printInfo(dinfo);
        }

        //! map a global dof index from a child
        /**
         * Given the index of a dof in the global dof-vector of one of the
         * children, compute the index of the same dof in the global
         * dof-vector of this ordering.
         *
         * \note update() must have been called before this method may be
         *       used.
         */
        SizeType subMap(SizeType child, SizeType indexInChild) const {
          static const SizeType *size = Size<Tag>::value;
          static const SizeType *offset = Offset<Tag>::value;

          return indexInChild % size[child]
            + offset[child]
            + (indexInChild/size[child]) * offset[Imp::CHILDREN];
        }

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
      };
    }

    //! Interface for merging index spaces
    template<class SizeType, class Tag, class Child, std::size_t k>
    class PowerBlockwiseOrdering :
      public TypeTree::PowerNode<Child, k>,
      public BlockwiseOrderingImp::Base<
        SizeType, Tag, PowerBlockwiseOrdering<SizeType, Tag, Child, k>
      >
    {
      typedef TypeTree::PowerNode<Child, k> Node;

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

      std::string name() const { return "CompositeBlockwiseOrdering"; }
    };

    template<std::size_t s0, std::size_t s1, std::size_t s2, std::size_t s3,
             std::size_t s4, std::size_t s5, std::size_t s6, std::size_t s7,
             std::size_t s8, std::size_t s9>
    struct TransformCompositeGFSToOrdering<
      ComponentBlockwiseOrderingTag<s0, s1, s2, s3, s4, s5, s6, s7, s8, s9>
      >
    {
      template<class GFSTraits, DUNE_TYPETREE_COMPOSITENODE_TEMPLATE_CHILDREN>
      struct result {
        typedef CompositeBlockwiseOrdering<
          typename GFSTraits::SizeType,
          ComponentBlockwiseOrderingTag<s0, s1, s2, s3, s4, s5, s6, s7, s8,
                                        s9>,
          DUNE_TYPETREE_COMPOSITENODE_CHILDTYPES
          > type;
      };
    };

   //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_BLOCKWISEORDERING_HH
