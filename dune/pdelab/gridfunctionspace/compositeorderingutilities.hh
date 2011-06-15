// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_COMPOSITEORDERINGUTILITIES_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_COMPOSITEORDERINGUTILITIES_HH

#include <dune/common/typetraits.hh>

#include <dune/pdelab/common/typetree/compositenodemacros.hh>
#include <dune/pdelab/common/typetree/transformation.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup GridFunctionSpace
    //! \ingroup PDELab
    //! \{

    //! Tag for the intermediate base class of the CompositeGridFunctionSpace.
    struct CompositeGridFunctionSpaceBaseTag {};


    //! Intermediate base class for the CompositeGridFunctionSpace.
    /**
     * The only purpose of this base class is to serve as a source tree for
     * the generation of the Ordering object corresponding to the
     * CompositeGridFunctionSpace.
     */
    template<DUNE_TYPETREE_COMPOSITENODE_TEMPLATE_CHILDREN>
    class CompositeGridFunctionSpaceBase
      : public DUNE_TYPETREE_COMPOSITENODE_BASETYPE
    {
      typedef DUNE_TYPETREE_COMPOSITENODE_BASETYPE NodeT;

    public:
      typedef CompositeGridFunctionSpaceBaseTag ImplementationTag;

    protected:
      CompositeGridFunctionSpaceBase(DUNE_TYPETREE_COMPOSITENODE_CONSTRUCTOR_SIGNATURE)
        : NodeT(DUNE_TYPETREE_COMPOSITENODE_CHILDVARIABLES)
      {}
    };


    //! The transformation from GFS to the ordering described by OrderingTag.
    template<typename GFS, typename OrderingTag>
    struct gfs_to_ordering
    {
      typedef GFS GridFunctionSpace;
      typedef OrderingTag Ordering;

      const GFS& asGridFunctionSpace(const typename GFS::NodeT& gfsnode) const
      {
        return static_cast<const GFS&>(gfsnode);
      }
    };

    //! Existing GridFunctionSpaces are converted to an Ordering by simply returning
    //! the ordering they already contain and not recursing deeper down the tree.
    template<typename ChildGFS>
    struct OrderingExtractingTransformation
    {

      // The ordering hierarchy for this subtree has already been constructed,
      // so there is no need for more recursion.
      static const bool recursive = false;

      typedef typename ChildGFS::Ordering transformed_type;
      typedef shared_ptr<transformed_type> transformed_storage_type;

      template<typename Transformation>
      static transformed_storage_type transform_storage(shared_ptr<const ChildGFS> cgfsp, const Transformation& t)
      {
        return cgfsp->orderingPtr();
      }

    };

    //! Register the default transformation descriptor for gfs_to_ordering, which will just extract the existing ordering.
    template<typename ChildGFS, typename GFS, typename OrderingTag, typename ImplementationTag>
    OrderingExtractingTransformation<ChildGFS>
    lookupNodeTransformation(ChildGFS*, gfs_to_ordering<GFS,OrderingTag>*, ImplementationTag);

   //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_COMPOSITEORDERINGUTILITIES_HH
