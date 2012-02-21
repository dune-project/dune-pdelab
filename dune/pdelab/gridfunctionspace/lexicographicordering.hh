// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_LEXICOGRAPHICORDERING_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_LEXICOGRAPHICORDERING_HH

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
//#include <dune/pdelab/gridfunctionspace/nonleaforderingbase.hh>
//#include <dune/pdelab/gridfunctionspace/orderingbase.hh>
//#include <dune/pdelab/gridfunctionspace/compositeorderingutilities.hh>

#include <dune/pdelab/gridfunctionspace/tags.hh>
#include <dune/pdelab/gridfunctionspace/orderingutility.hh>
#include <dune/pdelab/gridfunctionspace/orderingdynamicbase.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup GridFunctionSpace
    //! \ingroup PDELab
    //! \{

    //! \brief Indicate lexicographic ordering of the unknowns of non-leaf
    //!        grid function spaces.
    /**
     * This class instructs the non-leaf GridFunctionSpaces to order the dofs
     * of the child-GridFunctionSpaces in a lexicographic manner in the
     * combined dof-vector, i.e. first all dofs of child 0, then all dofs of
     * child 1 and so on.
     */
    struct LexicographicOrderingTag { };

    namespace lexicographic_ordering {

      //! Interface for merging index spaces
      template<typename DI, typename CI, typename Backend, typename Node>
      class Base
        : public OrderingBase<DI,CI>
      {

      public:

        typedef LexicographicOrderingTag OrderingTag;

        //! Construct ordering object
        /**
         * In general, an ordering object is not properly setup after
         * construction.  This must be done by a seperate call to update()
         * after all the children have been properly set up.
         */
        Base(Node& node)
        : OrderingBase<DI,CI>(node,false,nullptr)
        {
        }

        /*        Base(shared_ptr<Backend> backend)
          : OrderingBase<MI,CI>(static_cast<Node&>(*this),backend->blocked())
          , _backend(backend)
        {
        }*/

        //! update internal data structures
        /**
         * This method must be called after initialization and every time the
         * structure of the dof-vector of one of gfs's children changes.  All
         * the children must have been set up properly before the call to
         * update().
         */
        void updatae() {
          /*
          Dune::dinfo << asImp().name() << ":" << std::endl;

          TypeTree::applyToTree(asImp(),update_children());
          static_cast<OrderingBase<DI,CI>*>(this)->update();

          printInfo(dinfo);
          */
        }

        //! whether dofs are blocked per entity/intersection (they are not)
        bool blocked() const { return false; }

        //! map a global dof index from a child
        /**
         * Given the index of a dof in the global dof-vector of one of the
         * children, compute the index of the same dof in the global
         * dof-vector of this ordering.
         *
         * \note update() must have been called before this method may be
         *       used.
         */

        //! number of indices in this ordering
#if 0
        //! \brief offset of the block of dofs attached to a given entity (of
        //!        arbitrary codimension)
        /**
         * This implementation just throws NotImplemented since there are no
         * per-entity blocks for lexicographic ordering.
         *
         * \throw NotImplemented If this EntityType is not supported by the
         *                       ordering.
         */
        template<class Entity>
        SizeType entityOffset(const Entity &e) const {
          DUNE_THROW(NotImplemented, className<Imp>() << "::entityOffset() "
                     "does not make sense since the ordering is non-blocking");
        }
        //! \brief offset of the blocks of dofs attached to a given subentity
        //!        of an element
        /**
         * This implementation just throws NotImplemented since there are no
         * per-entity blocks for lexicographic ordering.
         *
         * \throw NotImplemented If this EntityType is not supported by the
         *                       ordering.
         */
        template<class Element>
        SizeType entityOffset(const Element &e, std::size_t codim,
                              std::size_t subentity) const {
          DUNE_THROW(NotImplemented, className<Imp>() << "::entityOffset() "
                     "does not make sense since the ordering is non-blocking");
        }
        //! offset of the block of dofs attached to a given intersection
        /**
         * This implementation just throws NotImplemented since there are no
         * per-intersection blocks for lexicographic ordering.
         *
         * \throw NotImplemented If this EntityType is not supported by the
         *                       ordering.
         */
        template<class Intersection>
        SizeType intersectionOffset(const Intersection &i) const {
          DUNE_THROW(NotImplemented,
                     className<Imp>() << "::intersectionOffset() does not "
                     "make sense since the ordering is non-blocking");
        }

#endif // 0

      };
    }

    //! Interface for merging index spaces
    template<typename DI, typename CI, typename Backend, typename Child, std::size_t k>
    class PowerLexicographicOrdering
      : public TypeTree::PowerNode<Child, k>
      , public lexicographic_ordering::Base<DI,
                                            CI,
                                            Backend,
                                            PowerLexicographicOrdering<DI,CI,Backend,Child,k>
                                            >
    {
      typedef TypeTree::PowerNode<Child, k> Node;

      typedef lexicographic_ordering::Base<DI,
                                           CI,
                                           Backend,
                                           PowerLexicographicOrdering<DI,CI,Backend,Child,k>
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
      PowerLexicographicOrdering(const typename Node::NodeStorage& children)
        : Node(children)
        , Base(*this)
      { }

      void update()
      {
        std::cout << "********************" << std::endl;
        for (std::size_t i = 0; i < Node::CHILDREN; ++i)
          {
            this->child(i).update();
            for (auto it = this->child(i)._child_offsets.begin(); it != this->child(i)._child_offsets.end(); ++it)
              std::cout << *it << " ";
            std::cout << std::endl;
          }
        Base::update();
        for (auto it = this->_child_offsets.begin(); it != this->_child_offsets.end(); ++it)
          std::cout << *it << " ";
        std::cout << std::endl;
        std::cout << "********************" << std::endl;
      }

      std::string name() const { return "PowerLexicographicOrdering"; }
    };


    template<typename GFS, typename Transformation, typename OrderingTag>
    struct power_gfs_to_ordering_descriptor;

    template<typename GFS, typename Transformation>
    struct power_gfs_to_ordering_descriptor<GFS,Transformation,LexicographicOrderingTag>
    {

      static const bool recursive = true;

      template<typename TC>
      struct result
      {

        typedef PowerLexicographicOrdering<
          typename Transformation::DOFIndex,
          typename Transformation::ContainerIndex,
          typename GFS::Traits::Backend,
          TC,
          GFS::CHILDREN
          > type;

        typedef shared_ptr<type> storage_type;

      };

      template<typename TC>
      static typename result<TC>::type transform(const GFS& gfs, const Transformation& t, const array<shared_ptr<TC>,GFS::CHILDREN>& children)
      {
        return typename result<TC>::type(children);
      }

      template<typename TC>
      static typename result<TC>::storage_type transform_storage(shared_ptr<const GFS> gfs, const Transformation& t, const array<shared_ptr<TC>,GFS::CHILDREN>& children)
      {
        return make_shared<typename result<TC>::type>(children);
      }

    };


    template<typename GridFunctionSpace, typename Params>
    power_gfs_to_ordering_descriptor<
      GridFunctionSpace,
      gfs_to_ordering<Params>,
      typename GridFunctionSpace::OrderingTag
      >
    lookupNodeTransformation(GridFunctionSpace* gfs, gfs_to_ordering<Params>* t, PowerGridFunctionSpaceTag tag);


    struct update_direct_children
      : public TypeTree::DirectChildrenVisitor
      , public TypeTree::DynamicTraversal
    {

      template<typename GFS, typename Child, typename TreePath, typename ChildIndex>
      void afterChild(const GFS& gfs, Child& child, TreePath tp, ChildIndex childIndex) const
      {
        child.update();
      }

    };


    //! Interface for merging index spaces
    template<typename DI, typename CI, typename Backend, DUNE_TYPETREE_COMPOSITENODE_TEMPLATE_CHILDREN>
    class CompositeLexicographicOrdering :
      public DUNE_TYPETREE_COMPOSITENODE_BASETYPE,
      public lexicographic_ordering::Base<DI,
                                          CI,
                                          Backend,
                                          CompositeLexicographicOrdering<
                                            DI,
                                            CI,
                                            Backend,
                                            DUNE_TYPETREE_COMPOSITENODE_CHILDTYPES
                                            >
                                          >
    {
      typedef DUNE_TYPETREE_COMPOSITENODE_BASETYPE Node;

      typedef lexicographic_ordering::Base<
        DI,
        CI,
        Backend,
        CompositeLexicographicOrdering<
          DI,
          CI,
          Backend,
          DUNE_TYPETREE_COMPOSITENODE_CHILDTYPES
          >
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
      CompositeLexicographicOrdering(DUNE_TYPETREE_COMPOSITENODE_STORAGE_CONSTRUCTOR_SIGNATURE)
        : Node(DUNE_TYPETREE_COMPOSITENODE_CHILDVARIABLES)
        , Base(*this)
      { }

      std::string name() const { return "CompositeLexicographicOrdering"; }

      void update()
      {
        std::cout << "********************" << std::endl;
        TypeTree::applyToTree(*this,update_direct_children());
        for (std::size_t i = 0; i < Node::CHILDREN; ++i)
          {
            for (auto it = this->dynamic_child(i)._child_offsets.begin(); it != this->dynamic_child(i)._child_offsets.end(); ++it)
              std::cout << *it << " ";
            std::cout << std::endl;
          }
        Base::update();
        for (auto it = this->_child_offsets.begin(); it != this->_child_offsets.end(); ++it)
          std::cout << *it << " ";
        std::cout << std::endl;
        std::cout << "********************" << std::endl;
      }
    };

#if HAVE_VARIADIC_TEMPLATES

    template<typename GFS, typename Transformation, typename OrderingTag>
    struct composite_gfs_to_ordering_descriptor;

    template<typename GFS, typename Transformation>
    struct composite_gfs_to_ordering_descriptor<GFS,Transformation,LexicographicOrderingTag>
    {

      static const bool recursive = true;

      template<typename... TC>
      struct result
      {

        typedef CompositeLexicographicOrdering<
          typename Transformation::DOFIndex,
          typename Transformation::ContainerIndex,
          typename GFS::Traits::Backend,
          TC...
          > type;

        typedef shared_ptr<type> storage_type;

      };

      template<typename... TC>
      static typename result<TC...>::type transform(const GFS& gfs, const Transformation& t, shared_ptr<TC>... children)
      {
        return typename result<TC...>::type(children...);
      }

      template<typename... TC>
      static typename result<TC...>::storage_type transform_storage(shared_ptr<const GFS> gfs, const Transformation& t, shared_ptr<TC>... children)
      {
        return make_shared<typename result<TC...>::type>(children...);
      }

    };


    template<typename GridFunctionSpace, typename Params>
    composite_gfs_to_ordering_descriptor<
      GridFunctionSpace,
      gfs_to_ordering<Params>,
      typename GridFunctionSpace::OrderingTag
      >
    lookupNodeTransformation(GridFunctionSpace* gfs, gfs_to_ordering<Params>* t, CompositeGridFunctionSpaceTag tag);


#else // HAVE_VARIADIC_TEMPLATES

    //! Node transformation descriptor for CompositeGridFunctionSpace -> LexicographicOrdering (without variadic templates).
    template<typename GFSNode, typename Transformation>
    struct CompositeGFSToLexicographicOrderingTransformation
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
        typedef CompositeLexicographicOrdering<typename Transformation::GridFunctionSpace::Traits::SizeType,
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
    CompositeGFSToLexicographicOrderingTransformation<GFSNode,gfs_to_ordering<GFS,LexicographicOrderingTag> >
    lookupNodeTransformation(GFSNode*, gfs_to_ordering<GFS,LexicographicOrderingTag>*, CompositeGridFunctionSpaceBaseTag);

#endif // HAVE_VARIADIC_TEMPLATES

   //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_LEXICOGRAPHICORDERING_HH
