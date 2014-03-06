// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_ORDERING_ORDERINGINTERFACE_HH
#define DUNE_PDELAB_ORDERING_ORDERINGINTERFACE_HH

#include <cstddef>

#include <dune/common/documentation.hh>

#include <dune/typetree/nodeinterface.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup Ordering
    //! \{

    //===============================================================
    // Utilities for the power and composite gfs
    // ===============================================================

    //! Interface for merging index spaces
    struct GridFunctionSpaceOrderingInterface :
      public TypeTree::NodeInterface
    {
      typedef ImplementationDefined SizeType;

      //! Construct ordering object
      /**
       * In general, an ordering object is not properly setup after
       * construction.  This must be done by a seperate call to update() after
       * all the children have been properly set up.
       *
       * \note This constructor must be present for ordering objects not at
       *       the leaf of the tree.
       */
      GridFunctionSpaceOrderingInterface
      (const TypeTree::NodeInterface::NodeStorage &children);

      //! Construct ordering object
      /**
       * In general, an ordering object is not properly setup after
       * construction.  This must be done by a seperate call to update() after
       * all the children have been properly set up.
       *
       * \note This constructor must be present for ordering objects at the
       *       leaf of the tree.
       */
      GridFunctionSpaceOrderingInterface();

      //! update internal data structures
      /**
       * This method must be called after initialization and every time the
       * structure of the dof-vector of one of gfs's children changes.  All
       * the children must have been set up properly before the call to
       * update().
       */
      void update();

      //! whether dofs are blocked per entity/intersection
      bool blocked() const;

      //! \brief whether all entites of the same geometry type/all
      //!        intersections have the same number of dofs
      /**
       * \note Even if fixedSize()==true the number of dofs may still vary
       *       between entities of different geometry type or between entities
       *       and intersections.
       */
      bool fixedSize() const;

      //! map a global dof index from a child
      /**
       * Given the index of a dof in the global dof-vector of one of the
       * children, compute the index of the same dof in the global dof-vector
       * of this ordering.
       *
       * \note update() must have been called before this method may be used.
       * \note This is only required on non-leafs.
       */
      SizeType subMap(SizeType child, SizeType indexInChild) const;

      //! number of indices in this ordering
      SizeType size() const;

      //! \brief maximum number of dofs attached to any given element and all
      //!        of its subentities and intersections
      /**
       * This is generally not an exact maximum and may be bigger than the
       * actual maximum.  There is however one special case: it is guaranteed
       * to be the exact maximum for fixedSize()==true.
       */
      SizeType maxLocalSize() const;

      //! \brief number of indices attached to a given entity (of arbitrary
      //!        codimension)
      /**
       * \note This method should be available for all types of entities
       *       required by the grid.  It is mostly required for communication,
       *       so if it is known that this method is not actually called for a
       *       given entity type the implementation may throw NotImplemented.
       * \note If the grid does not support a given entity type, it may still
       *       be possible to get this information using entitySize(const
       *       Element &e, std::size_t codim, std::size_t subentity).
       *
       * \throw NotImplemented If this EntityType is not supported by the
       *                       ordering.
       */
      template<class Entity>
      SizeType entitySize(const Entity &e) const;
      //! number of indices attached to a given subentity of an element
      /**
       * This method determines the number of indices attached to a subentity
       * of the given codim 0 entity.  If the grid (and the ordering) directly
       * supports entities of the given codimension, this is equivalent to
       * calling entitySize((*e.subEntity<codim>(subentity)).
       */
      template<class Element>
      SizeType entitySize(const Element &e, std::size_t codim,
                          std::size_t subentity) const;
      //! number of indices attached to a given intersection
      template<class Intersection>
      SizeType intersectionSize(const Intersection &i) const;

      //! \brief offset of the block of dofs attached to a given entity (of
      //!        arbitrary codimension)
      /**
       * \note This method should be available for all types of entities
       *       required by the grid.  It is mostly required for communcation,
       *       so if it is known that this method is not actually called for a
       *       given entity type the implementation may throw NotImplemented.
       * \note If the grid does not support a given entity type, it may still
       *       be possible to get this information using entityOffset(const
       *       Element &e, std::size_t codim, std::size_t subentity).
       *
       * \throw NotImplemented        If this EntityType is not supported by
       *                              the ordering.
       * \throw InvalidStateException If blocked()==false.
       */
      template<class Entity>
      SizeType entityOffset(const Entity &e) const;
      //! \brief offset of the blocks of dofs attached to a given subentity of
      //!        an element
      /**
       * This method determines the starting offset of the block of dofs
       * attached to a subentity of the given codim 0 entity.  If the grid
       * (and the ordering) directly support entities of the given
       * codimension, this is equivalent to calling
       * entityOffset(*e.subEntity<codim>(subentity)).
       */
      template<class Element>
      SizeType entityOffset(const Element &e, std::size_t codim,
                            std::size_t subentity) const;
      //! offset of the block of dofs attached to a given intersection
      template<class Intersection>
      SizeType intersectionOffset(const Intersection &i) const;
    };

   //! \} group Ordering
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ORDERING_ORDERINGINTERFACE_HH
