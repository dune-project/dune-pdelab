// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_COMPOSITEGRIDFUNCTIONSPACE_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_COMPOSITEGRIDFUNCTIONSPACE_HH

#include <dune/common/shared_ptr.hh>

#include <dune/pdelab/common/typetree/compositenodemacros.hh>
#include <dune/pdelab/common/typetree/utility.hh>
#include <dune/pdelab/gridfunctionspace/powercompositegridfunctionspacebase.hh>
#include <dune/pdelab/gridfunctionspace/datahandleprovider.hh>
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
    template<typename Backend,
             typename OrderingTag,
             DUNE_TYPETREE_COMPOSITENODE_TEMPLATE_CHILDREN>
    class CompositeGridFunctionSpace
      : public DUNE_TYPETREE_COMPOSITENODE_BASETYPE
      , public PowerCompositeGridFunctionSpaceBase<
          CompositeGridFunctionSpace<
            Backend,
            OrderingTag,
            DUNE_TYPETREE_COMPOSITENODE_CHILDTYPES>,
          typename DUNE_TYPETREE_COMPOSITENODE_BASETYPE::template Child<0>::
            Type::Traits::GridViewType,
          Backend,
          OrderingTag,
          DUNE_TYPETREE_COMPOSITENODE_BASETYPE::CHILDREN
        >
      , public DataHandleProvider<CompositeGridFunctionSpace<Backend,OrderingTag,DUNE_TYPETREE_COMPOSITENODE_CHILDTYPES> >
    {
      typedef DUNE_TYPETREE_COMPOSITENODE_BASETYPE NodeT;

      typedef PowerCompositeGridFunctionSpaceBase<
        CompositeGridFunctionSpace,
        typename NodeT::template Child<0>::Type::Traits::GridViewType,
        Backend,
        OrderingTag,
        NodeT::CHILDREN> ImplementationBase;

      friend class PowerCompositeGridFunctionSpaceBase<
        CompositeGridFunctionSpace,
        typename NodeT::template Child<0>::Type::Traits::GridViewType,
        Backend,
        OrderingTag,
        NodeT::CHILDREN>;

      typedef TypeTree::TransformTree<CompositeGridFunctionSpace,
                                      gfs_to_ordering<CompositeGridFunctionSpace>
                                      > ordering_transformation;

    public:
      typedef CompositeGridFunctionSpaceTag ImplementationTag;

      typedef typename ordering_transformation::Type Ordering;

      typedef typename ImplementationBase::Traits Traits;

      CompositeGridFunctionSpace(const Backend& backend, DUNE_TYPETREE_COMPOSITENODE_CONSTRUCTOR_SIGNATURE)
        : NodeT(DUNE_TYPETREE_COMPOSITENODE_CHILDVARIABLES_THROUGH_FUNCTION(TypeTree::assertGridViewType<typename NodeT::template Child<0>::Type>))
        , ImplementationBase(backend)
        , _ordering(nullptr)
      { }

      CompositeGridFunctionSpace(DUNE_TYPETREE_COMPOSITENODE_CONSTRUCTOR_SIGNATURE)
        : NodeT(DUNE_TYPETREE_COMPOSITENODE_CHILDVARIABLES_THROUGH_FUNCTION(TypeTree::assertGridViewType<typename NodeT::template Child<0>::Type>))
        , ImplementationBase(Backend())
        , _ordering(nullptr)
      { }

      //! Direct access to the DOF ordering.
      const Ordering &ordering() const
      {
        return *orderingStorage();
      }

      //! Direct access to the DOF ordering.
      Ordering &ordering()
      {
        return *orderingStorage();
      }

      //! Direct access to the storage of the DOF ordering.
      shared_ptr<const Ordering> orderingStorage() const
      {
        if (!_ordering)
          {
            _ordering = make_shared<Ordering>(ordering_transformation::transform(*this));
            _ordering->update();
          }
        return _ordering;
      }

      //! Direct access to the storage of the DOF ordering.
      shared_ptr<Ordering> orderingStorage()
      {
        if (!_ordering)
          {
            _ordering = make_shared<Ordering>(ordering_transformation::transform(*this));
            _ordering->update();
          }
        return _ordering;
      }

    private:
      mutable shared_ptr<Ordering> _ordering;

    };

    //! \}

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_COMPOSITEGRIDFUNCTIONSPACE_HH
