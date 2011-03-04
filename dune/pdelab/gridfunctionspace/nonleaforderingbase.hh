// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_NONLEAFORDERINGBASE_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_NONLEAFORDERINGBASE_HH

#include <cstddef>
#include <ostream>

#include <dune/pdelab/common/typetree/traversal.hh>
#include <dune/pdelab/common/typetree/visitor.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup GridFunctionSpace
    //! \ingroup PDELab
    //! \{

    //===============================================================
    // Utilities for the power and composite gfs
    // ===============================================================

    namespace NonLeafOrderingImp {
      template<class Stream>
      class PrintSizesVisitor :
        public TypeTree::DirectChildrenVisitor,
        public TypeTree::DynamicTraversal
      {
        Stream &stream;

      public:
        PrintSizesVisitor(Stream &stream_) : stream(stream_) {}

        template<class T, class TreePath>
        void pre(const T &t, const TreePath &) const
        { stream << t.maxLocalSize() << " "; }
      };

      class BlockedVisitor :
        public TypeTree::DirectChildrenVisitor,
        public TypeTree::DynamicTraversal
      {
        bool &blocked;

      public:
        BlockedVisitor(bool &blocked_) : blocked(blocked_) { }

        template<class T, class TreePath>
        void pre(const T &t, const TreePath &) const
        { blocked = blocked && t.blocked(); }
      };

      class FixedSizeVisitor :
        public TypeTree::DirectChildrenVisitor,
        public TypeTree::DynamicTraversal
      {
        bool &fixedSize;

      public:
        FixedSizeVisitor(bool &fixedSize_) : fixedSize(fixedSize_) { }

        template<class T, class TreePath>
        void pre(const T &t, const TreePath &) const
        { fixedSize = fixedSize && t.fixedSize(); }
      };

      template<class SizeType>
      class SizeVisitor :
        public TypeTree::DirectChildrenVisitor,
        public TypeTree::DynamicTraversal
      {
        SizeType &size;

      public:
        SizeVisitor(SizeType &size_) : size(size_) { }

        template<class T, class TreePath>
        void pre(const T &t, const TreePath &) const
        { size += t.size(); }
      };

      template<class SizeType>
      class MaxLocalSizeVisitor :
        public TypeTree::DirectChildrenVisitor,
        public TypeTree::DynamicTraversal
      {
        SizeType &size;

      public:
        MaxLocalSizeVisitor(SizeType &size_) : size(size_) { }

        template<class T, class TreePath>
        void pre(const T &t, const TreePath &) const
        { size += t.maxLocalSize(); }
      };

      template<class SizeType, class Entity>
      class EntitySizeVisitor :
        public TypeTree::DirectChildrenVisitor,
        public TypeTree::DynamicTraversal
      {
        SizeType &size;
        const Entity &e;

      public:
        EntitySizeVisitor(SizeType &size_, const Entity &e_) :
          size(size_), e(e_)
        { }

        template<class T, class TreePath>
        void pre(const T &t, const TreePath &) const
        { size += t.entitySize(e); }
      };

      template<class SizeType, class Element>
      class SubentitySizeVisitor :
        public TypeTree::DirectChildrenVisitor,
        public TypeTree::DynamicTraversal
      {
        SizeType &size;
        const Element &e;
        std::size_t codim;
        std::size_t subentity;

      public:
        SubentitySizeVisitor(SizeType &size_, const Element &e_,
                             std::size_t codim_, std::size_t subentity_) :
          size(size_), e(e_), codim(codim_), subentity(subentity_)
        { }

        template<class T, class TreePath>
        void pre(const T &t, const TreePath &) const
        { size += t.entitySize(e, codim, subentity); }
      };

      template<class SizeType, class Intersection>
      class IntersectionSizeVisitor :
        public TypeTree::DirectChildrenVisitor,
        public TypeTree::DynamicTraversal
      {
        SizeType &size;
        const Intersection &i;

      public:
        IntersectionSizeVisitor(SizeType &size_, const Intersection &i_) :
          size(size_), i(i_)
        { }

        template<class T, class TreePath>
        void pre(const T &t, const TreePath &) const
        { size += t.intersectionSize(i); }
      };
    }

    //! Interface for merging index spaces
    template<class SizeType_, class Imp>
    class NonLeafOrderingBase {
    protected:
      const Imp &asImp() const { return *static_cast<const Imp*>(this); }
      Imp &asImp() { return *static_cast<Imp*>(this); }

      template<class Stream>
      void printInfo(Stream &stream) {
        stream << "( ";
        TypeTree::applyToTree
          (asImp(), NonLeafOrderingImp::PrintSizesVisitor<Stream>(stream));
        stream << ") total size = " << asImp().size()
               << " max local size = " << asImp().maxLocalSize()
               << std::endl;
      }

    public:
      typedef SizeType_ SizeType;

      //! whether dofs are blocked per entity/intersection
      bool blocked() const {
        bool result = true;
        TypeTree::applyToTree(asImp(),
                              NonLeafOrderingImp::BlockedVisitor(result));
        return result;
      }

      //! \brief whether all entites of the same geometry type/all
      //!        intersections have the same number of dofs
      /**
       * This default implementation just checks whether all children have
       * fixedSize()==true.
       *
       * \note Even if fixedSize()==true the number of dofs may still vary
       *       between entities od different geometry type or between entities
       *       and intersections.
       */
      bool fixedSize() const {
        bool result = true;
        TypeTree::applyToTree(asImp(),
                              NonLeafOrderingImp::FixedSizeVisitor(result));
        return result;
      }

      //! number of indices in this ordering
      /**
       * This default implementation just adds the sizes of all children up
       * (which should be the right thing almost always).
       */
      SizeType size() const {
        SizeType result = 0;
        TypeTree::applyToTree
          (asImp(), NonLeafOrderingImp::SizeVisitor<SizeType>(result));
        return result;
      }

      //! \brief maximum number of dofs attached to any given element and all
      //!        of its subentities and intersections
      /**
       * This is generally not an exact maximum and may be bigger than the
       * actual maximum.  There is however one special case: it is guaranteed
       * to be the exact maximum for fixedSize()==true.
       */
      SizeType maxLocalSize() const {
        SizeType result = 0;
        TypeTree::applyToTree
          (asImp(), NonLeafOrderingImp::MaxLocalSizeVisitor<SizeType>(result));
        return result;
      }

      //! \brief number of indices attached to a given entity (of arbitrary
      //!        codimension)
      /**
       * This default implementation just adds the sizes of all children up
       * (which should be the right thing almost always).
       *
       * \note This method should be available for all types of entities
       *       required by the grid.  It is mostly required for communcation,
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
      SizeType entitySize(const Entity &e) const {
        SizeType result = 0;
        TypeTree::applyToTree
          (asImp(),
           NonLeafOrderingImp::EntitySizeVisitor<SizeType, Entity>(result, e));
        return result;
      }
      //! number of indices attached to a given subentity of an element
      /**
       * This default implementation just adds the sizes of all children up
       * (which should be the right thing almost always).
       *
       * This method determines the number of indices attached to a subentity
       * of the given codim 0 entity.  If the grid (and the ordering) directly
       * supports entities of the given codimension, this is equivalent to
       * calling entitySize((*e.subEntity<codim>(subentity)).
       */
      template<class Element>
      SizeType entitySize(const Element &e, std::size_t codim,
                          std::size_t subentity) const {
        SizeType result = 0;
        TypeTree::applyToTree
          (asImp(),
           NonLeafOrderingImp::SubentitySizeVisitor<SizeType, Element>
             (result, e, codim, subentity));
        return result;
      }
      //! number of indices attached to a given intersection
      /**
       * This default implementation just adds the sizes of all children up
       * (which should be the right thing almost always).
       */
      template<class Intersection>
      SizeType intersectionSize(const Intersection &i) const {
        SizeType result = 0;
        TypeTree::applyToTree
          (asImp(),
           NonLeafOrderingImp::IntersectionSizeVisitor<SizeType, Intersection>
             (result, i));
        return result;
      }
    };

   //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_NONLEAFORDERINGBASE_HH
