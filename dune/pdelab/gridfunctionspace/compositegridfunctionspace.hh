// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_COMPOSITEGRIDFUNCTIONSPACE_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_COMPOSITEGRIDFUNCTIONSPACE_HH

#include <memory>

#include <dune/typetree/compositenode.hh>
#include <dune/typetree/utility.hh>

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
             typename... Children>
    class CompositeGridFunctionSpace
      : public TypeTree::CompositeNode<Children...>
      , public PowerCompositeGridFunctionSpaceBase<
          CompositeGridFunctionSpace<
            Backend,
            OrderingTag,
            Children...>,
          typename TypeTree::Child<TypeTree::CompositeNode<Children...>,0>::Traits::EntitySet,
          Backend,
          OrderingTag,
          sizeof...(Children)
        >
      , public DataHandleProvider<CompositeGridFunctionSpace<Backend,OrderingTag,Children...> >
    {
      typedef TypeTree::CompositeNode<Children...> NodeT;

      typedef PowerCompositeGridFunctionSpaceBase<
        CompositeGridFunctionSpace,
        typename TypeTree::Child<NodeT,0>::Traits::EntitySet,
        Backend,
        OrderingTag,
        sizeof...(Children)> ImplementationBase;

      friend class PowerCompositeGridFunctionSpaceBase<
        CompositeGridFunctionSpace,
        typename TypeTree::Child<NodeT,0>::Traits::EntitySet,
        Backend,
        OrderingTag,
        sizeof...(Children)>;

      typedef TypeTree::TransformTree<CompositeGridFunctionSpace,
                                      gfs_to_ordering<CompositeGridFunctionSpace>
                                      > ordering_transformation;

      template<typename,typename>
      friend class GridFunctionSpaceBase;

    public:
      typedef CompositeGridFunctionSpaceTag ImplementationTag;

      typedef typename ordering_transformation::Type Ordering;

      typedef typename ImplementationBase::Traits Traits;

      // ********************************************************************************
      // constructors for stack-constructed children passed in by reference
      // ********************************************************************************

      CompositeGridFunctionSpace(const Backend& backend, Children&... children)
        : NodeT(TypeTree::assertGridViewType<typename NodeT::template Child<0>::Type>(children)...)
        , ImplementationBase(backend,OrderingTag())
      { }

      CompositeGridFunctionSpace(const OrderingTag& ordering_tag, Children&... children)
        : NodeT(TypeTree::assertGridViewType<typename NodeT::template Child<0>::Type>(children)...)
        , ImplementationBase(Backend(),ordering_tag)
      { }

      CompositeGridFunctionSpace(const Backend& backend, const OrderingTag& ordering_tag, Children&... children)
        : NodeT(TypeTree::assertGridViewType<typename NodeT::template Child<0>::Type>(children)...)
        , ImplementationBase(backend,ordering_tag)
      { }

      CompositeGridFunctionSpace(Children&... children)
        : NodeT(TypeTree::assertGridViewType<typename NodeT::template Child<0>::Type>(children)...)
        , ImplementationBase(Backend(),OrderingTag())
      { }

      // ********************************************************************************
      // constructors for heap-constructed children passed in as shared_ptrs
      // ********************************************************************************

      CompositeGridFunctionSpace(const Backend& backend, std::shared_ptr<Children>... children)
        : NodeT(children...)
        , ImplementationBase(backend,OrderingTag())
      { }

      CompositeGridFunctionSpace(const OrderingTag& ordering_tag, std::shared_ptr<Children>... children)
        : NodeT(children...)
        , ImplementationBase(Backend(),ordering_tag)
      { }

      CompositeGridFunctionSpace(const Backend& backend, const OrderingTag& ordering_tag, std::shared_ptr<Children>... children)
        : NodeT(children...)
        , ImplementationBase(backend,ordering_tag)
      { }

      CompositeGridFunctionSpace(std::shared_ptr<Children>... children)
        : NodeT(children...)
        , ImplementationBase(Backend(),OrderingTag())
      { }


      //! Direct access to the DOF ordering.
      const Ordering &ordering() const
      {
        if (!this->isRootSpace())
          {
            DUNE_THROW(GridFunctionSpaceHierarchyError,
                       "Ordering can only be obtained for root space in GridFunctionSpace tree.");
          }
        if (!_ordering)
          {
            create_ordering();
            this->update(*_ordering);
          }
        return *_ordering;
      }

      //! Direct access to the DOF ordering.
      Ordering &ordering()
      {
        if (!this->isRootSpace())
          {
            DUNE_THROW(GridFunctionSpaceHierarchyError,
                       "Ordering can only be obtained for root space in GridFunctionSpace tree.");
          }
        if (!_ordering)
          {
            create_ordering();
            this->update(*_ordering);
          }
        return *_ordering;
      }

      //! Direct access to the storage of the DOF ordering.
      std::shared_ptr<const Ordering> orderingStorage() const
      {
        if (!this->isRootSpace())
          {
            DUNE_THROW(GridFunctionSpaceHierarchyError,
                       "Ordering can only be obtained for root space in GridFunctionSpace tree.");
          }
        if (!_ordering)
          {
            create_ordering();
            this->update(*_ordering);
          }
        return _ordering;
      }

      //! Direct access to the storage of the DOF ordering.
      std::shared_ptr<Ordering> orderingStorage()
      {
        if (!this->isRootSpace())
          {
            DUNE_THROW(GridFunctionSpaceHierarchyError,
                       "Ordering can only be obtained for root space in GridFunctionSpace tree.");
          }
        if (!_ordering)
          {
            create_ordering();
            this->update(*_ordering);
          }
        return _ordering;
      }


    private:

      // This method here is to avoid a double update of the Ordering when the user calls
      // GFS::update() before GFS::ordering().
      void create_ordering() const
      {
        _ordering = std::make_shared<Ordering>(ordering_transformation::transform(*this));
      }

      mutable std::shared_ptr<Ordering> _ordering;

    };

    //! \}

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_COMPOSITEGRIDFUNCTIONSPACE_HH
