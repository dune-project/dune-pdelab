// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_COMPOSITEGRIDFUNCTIONSPACE_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_COMPOSITEGRIDFUNCTIONSPACE_HH

#include <dune/common/shared_ptr.hh>

#include <dune/pdelab/common/typetree/compositenodemacros.hh>
#include <dune/pdelab/common/typetree/utility.hh>
#include <dune/pdelab/gridfunctionspace/lexicographicordering.hh>
#include <dune/pdelab/gridfunctionspace/powercompositegridfunctionspacebase.hh>
#include <dune/pdelab/gridfunctionspace/tags.hh>

namespace Dune {
  namespace PDELab {

    //=======================================
    // composite grid function space
    //=======================================

    //! \addtogroup GridFunctionSpace grid function space
    //! \ingroup PDELab
    //! \{

    /** \brief base class for tuples of grid function spaces
        base class that holds implementation of the methods
        this is the default version with lexicographic ordering
        \tparam Mapper is the ordering parameter. Use e.g.
        \link GridFunctionSpaceLexicographicMapper GridFunctionSpaceLexicographicMapper \endlink
        or \link  GridFunctionSpaceComponentBlockwiseMapper  GridFunctionSpaceComponentBlockwiseMapper \endlink
        or \link  GridFunctionSpaceBlockwiseMapper  GridFunctionSpaceBlockwiseMapper \endlink
        or \link  GridFunctionSpaceDynamicBlockwiseMapper  GridFunctionSpaceDynamicBlockwiseMapper \endlink
        \tparam Ti are all grid function spaces
    */
    template<typename OrderingTag,
             DUNE_TYPETREE_COMPOSITENODE_TEMPLATE_CHILDREN>
    class CompositeGridFunctionSpace
      : public DUNE_TYPETREE_COMPOSITENODE_BASETYPE
      , public PowerCompositeGridFunctionSpaceBase<
          CompositeGridFunctionSpace<
            OrderingTag,
            DUNE_TYPETREE_COMPOSITENODE_CHILDTYPES>,
          typename DUNE_TYPETREE_COMPOSITENODE_BASETYPE::template Child<0>::
            Type::Traits::GridViewType,
          typename DUNE_TYPETREE_COMPOSITENODE_BASETYPE::template Child<0>::
            Type::Traits::BackendType,
          OrderingTag,
          DUNE_TYPETREE_COMPOSITENODE_BASETYPE::CHILDREN
        >
    {
      typedef DUNE_TYPETREE_COMPOSITENODE_BASETYPE BaseT;

      typedef PowerCompositeGridFunctionSpaceBase<
        CompositeGridFunctionSpace,
        typename BaseT::template Child<0>::Type::Traits::GridViewType,
        typename BaseT::template Child<0>::Type::Traits::BackendType,
        OrderingTag,
        BaseT::CHILDREN> ImplementationBase;

      friend class PowerCompositeGridFunctionSpaceBase<
        CompositeGridFunctionSpace,
        typename BaseT::template Child<0>::Type::Traits::GridViewType,
        typename BaseT::template Child<0>::Type::Traits::BackendType,
        OrderingTag,
        BaseT::CHILDREN>;

    public:
      typedef CompositeGridFunctionSpaceTag ImplementationTag;

      typedef typename TransformCompositeGFSToOrdering<OrderingTag>::
        template result<
          typename ImplementationBase::Traits,
          DUNE_TYPETREE_COMPOSITENODE_CHILDTYPES_NESTED_TYPE(Ordering)
        >::type Ordering;

      typedef typename ImplementationBase::Traits Traits;

      CompositeGridFunctionSpace(DUNE_TYPETREE_COMPOSITENODE_CONSTRUCTOR_SIGNATURE)
        : BaseT(DUNE_TYPETREE_COMPOSITENODE_CHILDVARIABLES_THROUGH_FUNCTION(TypeTree::assertGridViewType<typename BaseT::template Child<0>::Type>))
        , orderingp
          (make_shared<Ordering>
           (*this,
            DUNE_TYPETREE_COMPOSITENODE_CHILDVARIABLES_MEMBER(orderingPtr())))
      { }

      const Ordering &ordering() const { return *orderingp; }
      const shared_ptr<Ordering> &orderingPtr() { return orderingp; }

    private:
      shared_ptr<Ordering> orderingp;

    };

    //! \}

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_COMPOSITEGRIDFUNCTIONSPACE_HH
