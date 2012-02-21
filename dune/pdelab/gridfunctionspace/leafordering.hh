// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_LEAFORDERING_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_LEAFORDERING_HH

#include <cstddef>

#include <dune/pdelab/common/typetree/leafnode.hh>
#include <dune/pdelab/gridfunctionspace/gridviewordering.hh>
#include <dune/pdelab/gridfunctionspace/leaflocalordering.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup GridFunctionSpace
    //! \ingroup PDELab
    //! \{

    //! Dummy ordering for leaf gridfunctionspaces
    template<typename GVOrdering>
    class LeafOrdering
      : public TypeTree::VariadicCompositeNode<GVOrdering>
      , public OrderingBase<typename GVOrdering::Traits::MI,typename GVOrdering::Traits::CI>
    {

      typedef TypeTree::VariadicCompositeNode<GVOrdering> NodeT;

      typedef OrderingBase<typename GVOrdering::Traits::MI,typename GVOrdering::Traits::CI> BaseT;

    public:

      typedef ... Traits;

      //! Construct ordering object
      /**
       * In general, an ordering object is not properly setup after
       * construction.  This must be done by a seperate call to update().
       * This particular ordering however can be used right away.
       */
      LeafOrdering(const NodeT::NodeStorage& gv_ordering_storage)
        : NodeT(gv_ordering_storage)
        , BaseT(*this,this->template childStorage<0>.get())
      { }

      //! update internal data structures
      /**
       * In general this method must be called after initialization and every
       * time the structure of the dof-vector of one of gfs's children
       * changes.  For this particular ordering however this method does
       * nothing.
       */
      void update()
      {
        gridViewOrdering().update();
      }

      //! dofs are blocked per entity/intersection on the leafs
      bool blocked() const { return true; }

      //! \brief whether all entites of the same geometry type/all
      //!        intersections have the same number of dofs
      /**
       * On the leaf this is realized by iterating over the grid during update
       * an checking.
       *
       * \note Even if fixedSize()==true the number of dofs may still vary
       *       between entities od different geometry type or between entities
       *       and intersections.
       */
      bool fixedSize() const { return gridViewOrdering().fixedSize(); }

      //! number of indices in this ordering
      SizeType size() const { return gridViewOrdering().size(); }

      //! \brief maximum number of dofs attached to any given element and all
      //!        of its subentities and intersections
      /**
       * This is generally not an exact maximum and may be bigger than the
       * actual maximum.  There is however one special case: it is guaranteed
       * to be the exact maximum for fixedSize()==true.
       */
      SizeType maxLocalSize() const { return gridViewOrdering().maxLocalSize(); }

      //! \brief number of indices attached to a given entity (of arbitrary
      //!        codimension)
      /**
       * \note If the grid does not support a given entity type, it may still
       *       be possible to get this information using entitySize(const
       *       Element &e, std::size_t codim, std::size_t subentity).
       */
      template<class Entity>
      SizeType entitySize(const Entity &e) const { return gfs.entitySize(e); }
      //! number of indices attached to a given subentity of an element
      /**
       * This method determines the number of indices attached to a subentity
       * of the given codim 0 entity.  If the grid (and the ordering) directly
       * supports entities of the given codimension, this is equivalent to
       * calling entitySize((*e.subEntity<codim>(subentity)).
       */
      template<class Element>
      SizeType entitySize(const Element &e, std::size_t codim,
                          std::size_t subentity) const
      { return gfs.entitySize(e, codim, subentity); }
      //! number of indices attached to a given intersection
      template<class Intersection>
      SizeType intersectionSize(const Intersection &i) const
      { return gfs.intersectionSize(i); }

      //! \brief offset of the block of dofs attached to a given entity (of
      //!        arbitrary codimension)
      /**
       * \note If the grid does not support a given entity type, it may still
       *       be possible to get this information using entityOffset(const
       *       Element &e, std::size_t codim, std::size_t subentity).
       *
       * \throw NotImplemented        If this EntityType is not supported by
       *                              the ordering.
       * \throw InvalidStateException If blocked()==false.
       */
      template<class Entity>
      SizeType entityOffset(const Entity &e) const
      { return gfs.entityOffset(e); }
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
                            std::size_t subentity) const
      { return gfs.entityOffset(e, codim, subentity); }
      //! offset of the block of dofs attached to a given intersection
      template<class Intersection>
      SizeType intersectionOffset(const Intersection &i) const
      { return gfs.intersectionOffset(i); }
    };


    template<typename LeafGFS, typename Transformation>
    struct leaf_gfs_to_ordering_descriptor
    {

      static const bool recursive = false;

      typedef LeafLocalOrdering<LeafGFS,
                                typename Transformation::MultiIndex,
                                typename Transformation::ContainerIndex
                                > LocalOrdering;

      typedef LeafGridViewOrdering<typename LeafGFS::Traits::GridView,LocalOrdering> GridViewOrdering;

      typedef LeafOrdering<GridViewOrdering> transformed_type;
      typedef shared_ptr<transformed_type> transformed_storage_type;

      static transformed_type transform(const LeafGFS& s, const Transformation& t)
      {
        return transformed_type(make_tuple(make_shared<GridViewOrdering>(gfs.gridview(),make_tuple(make_shared<LocalOrdering>(gfs)))));
      }

      static transformed_storage_type transform_storage(shared_ptr<const LeafGFS> gfs, const Transformation& t)
      {
        return make_shared<transformed_type>(make_tuple(make_shared<GridViewOrdering>(gfs.gridview(),make_tuple(make_shared<LocalOrdering>(gfs)))));
      }

    };

    template<typename LeafGFS, typename Params>
    leaf_gfs_to_ordering_descriptor<LeafGFS,gfs_to_ordering<Params> >
    lookupNodeTransformation(LeafGFS* gfs, gfs_to_ordering<Params>* t, LeafGridFunctionSpaceTag tag);

   //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_LEAFORDERING_HH
