// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_DYNAMICBLOCKWISEORDERING_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_DYNAMICBLOCKWISEORDERING_HH

#include <cstddef>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include <dune/common/stdstreams.hh>

#include <dune/grid/common/mcmgmapper.hh>

#include <dune/pdelab/common/typetree/compositenodemacros.hh>
#include <dune/pdelab/common/typetree/powernode.hh>
#include <dune/pdelab/common/typetree/traversal.hh>
#include <dune/pdelab/common/typetree/visitor.hh>
#include <dune/pdelab/gridfunctionspace/nonleaforderingbase.hh>
#include <dune/pdelab/gridfunctionspace/orderingbase.hh>
#include <dune/pdelab/gridfunctionspace/compositeorderingutilities.hh>

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
    struct DynamicBlockwiseOrderingTag { };
    typedef DynamicBlockwiseOrderingTag
      GridFunctionSpaceDynamicBlockwiseMapper;

    namespace DynamicBlockwiseOrderingImp {

      template<class SizeType>
      class CollectSizesVisitor :
        public TypeTree::DirectChildrenVisitor,
        public TypeTree::DynamicTraversal
      {
        typedef std::vector<std::vector<std::pair<SizeType, SizeType> > >
          IndexRanges;
        IndexRanges &indexRanges;

      public:
        CollectSizesVisitor(IndexRanges &indexRanges_) :
          indexRanges(indexRanges_)
        { }

        // template<class T, class TreePath>
        // void pre(const T &t, const TreePath &tp) const
        // { indexRanges[tp.back()].back().first = t.size(); }

        template<class T, class Child, class TreePath, class ChildIndex>
        void beforeChild(const T &t, const Child& child, TreePath, ChildIndex childIndex)
        { indexRanges[childIndex].back().first = child.size(); }
      };

      template<class SizeType, class Element>
      class CollectOffsetsVisitor :
        public TypeTree::DirectChildrenVisitor,
        public TypeTree::DynamicTraversal
      {
        typedef std::vector<std::vector<std::pair<SizeType, SizeType> > >
          IndexRanges;
        IndexRanges &indexRanges;
        std::size_t eindex;
        const Element &e;

      public:
        CollectOffsetsVisitor(IndexRanges &indexRanges_, std::size_t eindex_,
                              const Element &e_) :
          indexRanges(indexRanges_), eindex(eindex_), e(e_)
        { }

        // template<class T, class TreePath>
        // void pre(const T &t, const TreePath &tp) const
        // { indexRanges[tp.back()][eindex].first = t.entityOffset(e); }

        template<class T, class Child, class TreePath, class ChildIndex>
        void beforeChild(const T &t, const Child& child, TreePath, ChildIndex childIndex)
        { indexRanges[childIndex][eindex].first = child.entityOffset(e); }
      };

      struct PrintMaxLocalSizesVisitor :
        public TypeTree::DirectChildrenVisitor,
        public TypeTree::DynamicTraversal
      {
        // template<class T, class TreePath>
        // void pre(const T &t, const TreePath &) const
        // { dinfo << t.maxLocalSize() << " "; }

        template<class T, class Child, class TreePath, class ChildIndex>
        void beforeChild(const T &t, const Child& child, TreePath, ChildIndex childIndex)
        { dinfo << child.maxLocalSize() << " "; }
      };

      //! Interface for merging index spaces
      template<class SizeType, class GV, class Imp>
      class Base :
        public NonLeafOrderingBase<SizeType, Imp>
      {
      protected:
        using NonLeafOrderingBase<SizeType, Imp>::asImp;
        using NonLeafOrderingBase<SizeType, Imp>::printInfo;

      private:
        const GV &gv;
        MultipleCodimMultipleGeomTypeMapper<GV, MCMGElementLayout> mapper;

        // first: childIndex, second: parentIndex
        typedef std::pair<SizeType, SizeType> IndexPair;
        typedef std::vector<IndexPair> IndexRangeList;
        std::vector<IndexRangeList> indexRanges;
        typedef typename IndexRangeList::const_iterator IndexRangeIterator;
        mutable std::vector<IndexRangeIterator> indexRangeIterators;

      public:

        Base(const Base &b) : gv(b.gv), mapper(b.mapper), indexRanges(b.indexRanges)
        {
          indexRangeIterators.clear();
          for(std::size_t child=0; child<Imp::CHILDREN; ++child)
            indexRangeIterators.push_back(indexRanges[child].begin());
        }

        Base(const GV &gv_) : gv(gv_), mapper(gv) { asImp().update(); }

        //! update internal data structures
        /**
         * This method must be called after initialization and every time the
         * structure of the dof-vector of one of gfs's children changes.  All
         * the children must have been set up properly before the call to
         * update().
         */
        void update() {
          dinfo << asImp().name() << ":" << std::endl;

          mapper.update();

          indexRanges.assign(Imp::CHILDREN, IndexRangeList());
          for(std::size_t i=0; i<Imp::CHILDREN; ++i)
            indexRanges[i].resize(mapper.size()+1);

          // collect all the per-element offsets of the children
          typedef typename GV::template Codim<0>::Iterator Iterator;
          typedef typename GV::template Codim<0>::Entity Element;
          const Iterator &eend = gv.template end<0>();
          for(Iterator eit = gv.template begin<0>(); eit != eend; ++eit)
            TypeTree::applyToTree(asImp(),
                                  CollectOffsetsVisitor<SizeType, Element>
                                    (indexRanges, mapper.map(*eit), *eit));
          // collect the sizes of the children in the last elements
          TypeTree::applyToTree(asImp(),
                                CollectSizesVisitor<SizeType>(indexRanges));

          // now assign indices in the parent to the child-blocks
          SizeType parentIndex = 0;
          for(std::size_t block = 0; block < std::size_t(mapper.size()); ++block)
            for(std::size_t child = 0; child < Imp::CHILDREN; ++child) {
              indexRanges[child][block].second = parentIndex;
              parentIndex += indexRanges[child][block+1].first -
                indexRanges[child][block].first;
            }
          // assign parent size to final entries
          for(std::size_t child = 0; child < Imp::CHILDREN; ++child)
            indexRanges[child].back().second = parentIndex;

          indexRangeIterators.clear();
          for(std::size_t child=0; child<Imp::CHILDREN; ++child)
            indexRangeIterators.push_back(indexRanges[child].begin());

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
          IndexRangeIterator &it = indexRangeIterators[child];

          while(it->first <= indexInChild){
            ++it;
            assert(it != indexRanges[child].end());
          }
          while(it->first > indexInChild){ 
            assert(it != indexRanges[child].begin());
            --it;
          }
          return it->second + (indexInChild - it->first);
        }

        //! number of indices in this ordering
        SizeType size() const { return indexRanges[0].back().second; }

        //! \brief number of indices attached to a given entity (of arbitrary
        //!        codimension)
        /**
         * \note This method should be available for all types of entities
         *       required by the grid.  It is mostly required for
         *       communication, so if it is known that this method is not
         *       actually called for a given entity type the implementation
         *       may throw NotImplemented.
         * \note If the grid does not support a given entity type, it may
         *       still be possible to get this information using
         *       entitySize(const Element &e, std::size_t codim, std::size_t
         *       subentity).
         *
         * \throw NotImplemented If this EntityType is not supported by the
         *                       ordering.
         */
        template<class Entity>
        SizeType entitySize(const Entity &e) const {
          if(Entity::codimension == 0) {
            std::size_t eindex = mapper.map(e);
            return indexRanges[0][eindex+1].second -
              indexRanges[0][eindex].second;
          }
          else
            return 0;
        }
        //! number of indices attached to a given subentity of an element
        /**
         * This method determines the number of indices attached to a
         * subentity of the given codim 0 entity.  If the grid (and the
         * ordering) directly supports entities of the given codimension, this
         * is equivalent to calling
         * entitySize((*e.subEntity<codim>(subentity)).
         */
        template<class Element>
        SizeType entitySize(const Element &e, std::size_t codim,
                            std::size_t subentity) const
        {
          if(codim == 0) {
            std::size_t eindex = mapper.map(e);
            return indexRanges[0][eindex+1].second -
              indexRanges[0][eindex].second;
          }
          else
            return 0;
        }
        //! number of indices attached to a given intersection
        template<class Intersection>
        SizeType intersectionSize(const Intersection &i) const { return 0; }

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
          if(Entity::codimension == 0)
            return indexRanges[0][mapper.map(e)].second;
          else
            return size();
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
                              std::size_t subentity) const
        {
          if(codim == 0)
            return indexRanges[0][mapper.map(e)].second;
          else
            return size();
        }
        //! offset of the block of dofs attached to a given intersection
        template<class Intersection>
        SizeType intersectionOffset(const Intersection &i) const
        { return size(); }
      };
    }

    //! Interface for merging index spaces
    template<class SizeType, class GV, class Child, std::size_t k>
    class PowerDynamicBlockwiseOrdering :
      public TypeTree::PowerNode<Child, k>,
      public DynamicBlockwiseOrderingImp::Base<
        SizeType, GV, PowerDynamicBlockwiseOrdering<SizeType, GV, Child, k>
      >
    {
      typedef TypeTree::PowerNode<Child, k> Node;
      typedef DynamicBlockwiseOrderingImp::Base<
        SizeType, GV, PowerDynamicBlockwiseOrdering
        > Base;

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
      PowerDynamicBlockwiseOrdering
      ( const GFS &gfs,
        const typename Node::NodeStorage &children) :
        Node(children), Base(gfs.gridview())
      { }

      std::string name() const { return "PowerDynamicBlockwiseOrdering"; }

    };

    template<>
    struct TransformPowerGFSToOrdering<DynamicBlockwiseOrderingTag> {
      template<class GFSTraits, class TransformedChild, std::size_t k>
      struct result {
        typedef PowerDynamicBlockwiseOrdering<typename GFSTraits::SizeType,
                                              typename GFSTraits::GridViewType,
                                              TransformedChild,
                                              k> type;
      };
    };

    //! Interface for merging index spaces
    template<class SizeType, class GV,
             DUNE_TYPETREE_COMPOSITENODE_TEMPLATE_CHILDREN>
    class CompositeDynamicBlockwiseOrdering :
      public DUNE_TYPETREE_COMPOSITENODE_BASETYPE,
      public DynamicBlockwiseOrderingImp::Base<
        SizeType, GV,
        CompositeDynamicBlockwiseOrdering<
          SizeType, GV, DUNE_TYPETREE_COMPOSITENODE_CHILDTYPES
          >
      >
    {
      typedef DUNE_TYPETREE_COMPOSITENODE_BASETYPE Node;
      typedef DynamicBlockwiseOrderingImp::Base<
        SizeType, GV, CompositeDynamicBlockwiseOrdering
        > Base;

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
      CompositeDynamicBlockwiseOrdering
      ( const GFS &gfs,
        DUNE_TYPETREE_COMPOSITENODE_STORAGE_CONSTRUCTOR_SIGNATURE) :
        Node(DUNE_TYPETREE_COMPOSITENODE_CHILDVARIABLES), Base(gfs.gridview())
      { }

      std::string name() const { return "CompositeDynamicBlockwiseOrdering"; }

    };

#if HAVE_VARIADIC_TEMPLATES

    //! Node transformation descriptor for CompositeGridFunctionSpace -> DynamicBlockwiseOrdering (with variadic templates).
    template<typename GFSNode, typename Transformation>
    struct VariadicCompositeGFSToDynamicBlockwiseOrderingTransformation
    {

      static const bool recursive = true;

      template<typename... TC>
      struct result
      {
        typedef CompositeDynamicBlockwiseOrdering<typename Transformation::GridFunctionSpace::Traits::SizeType,
                                                  typename Transformation::GridFunctionSpace::Traits::GridViewType,
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
    template<typename GFSNode, typename GFS>
    VariadicCompositeGFSToDynamicBlockwiseOrderingTransformation<GFSNode,gfs_to_ordering<GFS,DynamicBlockwiseOrderingTag> >
    lookupNodeTransformation(GFSNode*, gfs_to_ordering<GFS,DynamicBlockwiseOrderingTag>*, CompositeGridFunctionSpaceBaseTag);

#else // HAVE_VARIADIC_TEMPLATES

    //! Node transformation descriptor for CompositeGridFunctionSpace -> DynamicBlockwiseOrdering (without variadic templates).
    template<typename GFSNode, typename Transformation>
    struct CompositeGFSToDynamicBlockwiseOrderingTransformation
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
        typedef CompositeDynamicBlockwiseOrdering<typename Transformation::GridFunctionSpace::Traits::SizeType,
                                                  typename Transformation::GridFunctionSpace::Traits::GridViewType,
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

    // Register transformation
    template<typename GFSNode, typename GFS>
    CompositeGFSToDynamicBlockwiseOrderingTransformation<GFSNode,gfs_to_ordering<GFS,DynamicBlockwiseOrderingTag> >
    lookupNodeTransformation(GFSNode*, gfs_to_ordering<GFS,DynamicBlockwiseOrderingTag>*, CompositeGridFunctionSpaceBaseTag);

#endif // HAVE_VARIADIC_TEMPLATES


   //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_DYNAMICBLOCKWISEORDERING_HH
