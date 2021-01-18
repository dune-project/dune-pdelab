// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_DYNAMICPOWERGRIDFUNCTIONSPACE_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_DYNAMICPOWERGRIDFUNCTIONSPACE_HH

#include <cstddef>
#include <memory>

#include <dune/typetree/powernode.hh>

#include <dune/pdelab/gridfunctionspace/powercompositegridfunctionspacebase.hh>
#include <dune/pdelab/gridfunctionspace/datahandleprovider.hh>
#include <dune/pdelab/gridfunctionspace/tags.hh>

namespace Dune {
  namespace PDELab {

    //! Trait class for the multi component grid function spaces
    template<typename G, typename B, typename O>
    struct DynamicTreeGridFunctionSpaceTraits
    {
      enum{
        //! \brief True if this grid function space is composed of others.
        isComposite = 1,
        //! \brief True if it has a dynamic number of childen
        isDynamic = 1
      };

      using EntitySet = G;

      using GridView = typename EntitySet::GridView;

      //! \brief the grid view where grid function is defined upon
      using GridViewType = GridView;

      //! \brief vector backend
      typedef B BackendType;

      typedef B Backend;

      //! \brief mapper
      typedef O MapperType;

      typedef O OrderingTag;

      //! \brief short cut for size type exported by Backend
      typedef typename B::size_type SizeType;
    };

    //! Mixin class providing common functionality of PowerGridFunctionSpace and CompositeGridFunctionSpace
    template<typename GridFunctionSpace, typename GV, typename B, typename O>
    class DynamicTreeGridFunctionSpaceBase
      : public GridFunctionSpaceBase<
                 GridFunctionSpace,
                 DynamicTreeGridFunctionSpaceTraits<GV,B,O>
                 >
    {

#ifndef DOXYGEN

      const GridFunctionSpace& gfs() const
      {
        return static_cast<const GridFunctionSpace&>(*this);
      }

      GridFunctionSpace& gfs()
      {
        return static_cast<GridFunctionSpace&>(*this);
      }

#endif // DOXYGEN

    public:

      //! export traits class
      typedef DynamicTreeGridFunctionSpaceTraits<GV,B,O> Traits;

    private:

      typedef GridFunctionSpaceBase<GridFunctionSpace,Traits> BaseT;

    public:

      typedef O OrderingTag;

      //! extract type for storing constraints
      template<typename E>
      struct ConstraintsContainer
      {
        typedef typename std::conditional<
          std::is_same<
            typename GridFunctionSpace::ChildType::template ConstraintsContainer<E>::Type,
            EmptyTransformation
            >::value,
          EmptyTransformation,
          ConstraintsTransformation<
            typename GridFunctionSpace::Ordering::Traits::DOFIndex,
            typename GridFunctionSpace::Ordering::Traits::ContainerIndex,
            E
            >
          >::type Type;
      };

      //! get grid view
      const typename Traits::GridView& gridView () const
      {
        return gfs().child(0).gridView();
      }

      //! get grid view partition
      const typename Traits::EntitySet& entitySet () const
      {
        return gfs().child(0).entitySet();
      }

      DynamicTreeGridFunctionSpaceBase(const B& backend, const OrderingTag& ordering_tag)
        : BaseT(backend,ordering_tag)
      {}

    };

    //=======================================
    // dynamic power grid function space
    //=======================================

    /** \brief base class for tuples of grid function spaces
        product of identical grid function spaces
        base class that holds implementation of the methods

        PGFS(T,k) = {T}^k

        \tparam T the underlying are all grid function spaces
    */
    template<typename T,
             typename Backend,
             typename OrderingTag = LexicographicOrderingTag>
    class DynamicPowerGridFunctionSpace
      : public TypeTree::DynamicPowerNode<T>
      , public DynamicTreeGridFunctionSpaceBase<
          DynamicPowerGridFunctionSpace<T, Backend, OrderingTag>,
          typename T::Traits::EntitySet,
          Backend,
          OrderingTag>
      , public DataHandleProvider<DynamicPowerGridFunctionSpace<T,Backend,OrderingTag> >
    {

    public:

      typedef DynamicPowerGridFunctionSpaceTag ImplementationTag;

      typedef TypeTree::DynamicPowerNode<T> BaseT;

    private:

      typedef DynamicTreeGridFunctionSpaceBase<
        DynamicPowerGridFunctionSpace,
        typename T::Traits::EntitySet,
        Backend,
        OrderingTag
        > ImplementationBase;

      friend class DynamicTreeGridFunctionSpaceBase<
        DynamicPowerGridFunctionSpace,
        typename T::Traits::EntitySet,
        Backend,
        OrderingTag>;

      template<typename,typename>
      friend class GridFunctionSpaceBase;

      typedef TypeTree::TransformTree<DynamicPowerGridFunctionSpace,
                                      gfs_to_ordering<DynamicPowerGridFunctionSpace>
                                      > ordering_transformation;

    public:

      typedef typename ordering_transformation::Type Ordering;

      //! export traits class
      typedef typename ImplementationBase::Traits Traits;

      /**
       * @brief Construct a new Power Grid Function Space object
       *
       * @param container     array with pointers to child spaces
       * @param backend       backend object
       * @param ordering_tag  ordering tag object
       */
      DynamicPowerGridFunctionSpace(const std::vector<shared_ptr<T>>& container, const Backend& backend = Backend(), const OrderingTag ordering_tag = OrderingTag())
        : BaseT(container)
        , ImplementationBase(backend,ordering_tag)
      {}

      DynamicPowerGridFunctionSpace(T& c, const Backend& backend = Backend(), const OrderingTag ordering_tag = OrderingTag())
        : BaseT(c)
        , ImplementationBase(backend,ordering_tag)
      {}

      DynamicPowerGridFunctionSpace (T& c0,
                              T& c1,
                              const Backend& backend = Backend(),
                              const OrderingTag ordering_tag = OrderingTag())
        : BaseT(c0,c1)
        , ImplementationBase(backend,ordering_tag)
      {}

      DynamicPowerGridFunctionSpace (T& c0,
                              T& c1,
                              T& c2,
                              const Backend& backend = Backend(),
                              const OrderingTag ordering_tag = OrderingTag())
        : BaseT(c0,c1,c2)
        , ImplementationBase(backend,ordering_tag)
      {}

      DynamicPowerGridFunctionSpace (T& c0,
                              T& c1,
                              T& c2,
                              T& c3,
                              const Backend& backend = Backend(),
                              const OrderingTag ordering_tag = OrderingTag())
        : BaseT(c0,c1,c2,c3)
        , ImplementationBase(backend,ordering_tag)
      {}

      DynamicPowerGridFunctionSpace (T& c0,
                              T& c1,
                              T& c2,
                              T& c3,
                              T& c4,
                              const Backend& backend = Backend(),
                              const OrderingTag ordering_tag = OrderingTag())
        : BaseT(c0,c1,c2,c3,c4)
        , ImplementationBase(backend,ordering_tag)
      {}

      DynamicPowerGridFunctionSpace (T& c0,
                              T& c1,
                              T& c2,
                              T& c3,
                              T& c4,
                              T& c5,
                              const Backend& backend = Backend(),
                              const OrderingTag ordering_tag = OrderingTag())
        : BaseT(c0,c1,c2,c3,c4,c5)
        , ImplementationBase(backend,ordering_tag)
      {}

      DynamicPowerGridFunctionSpace (T& c0,
                              T& c1,
                              T& c2,
                              T& c3,
                              T& c4,
                              T& c5,
                              T& c6,
                              const Backend& backend = Backend(),
                              const OrderingTag ordering_tag = OrderingTag())
        : BaseT(c0,c1,c2,c3,c4,c5,c6)
        , ImplementationBase(backend,ordering_tag)
      {}

      DynamicPowerGridFunctionSpace (T& c0,
                              T& c1,
                              T& c2,
                              T& c3,
                              T& c4,
                              T& c5,
                              T& c6,
                              T& c7,
                              const Backend& backend = Backend(),
                              const OrderingTag ordering_tag = OrderingTag())
        : BaseT(c0,c1,c2,c3,c4,c5,c6,c7)
        , ImplementationBase(backend,ordering_tag)
      {}

      DynamicPowerGridFunctionSpace (T& c0,
                              T& c1,
                              T& c2,
                              T& c3,
                              T& c4,
                              T& c5,
                              T& c6,
                              T& c7,
                              T& c8,
                              const Backend& backend = Backend(),
                              const OrderingTag ordering_tag = OrderingTag())
        : BaseT(c0,c1,c2,c3,c4,c5,c6,c7,c8)
        , ImplementationBase(backend,ordering_tag)
      {}

      DynamicPowerGridFunctionSpace (T& c0,
                              T& c1,
                              T& c2,
                              T& c3,
                              T& c4,
                              T& c5,
                              T& c6,
                              T& c7,
                              T& c8,
                              T& c9,
                              const Backend& backend = Backend(),
                              const OrderingTag ordering_tag = OrderingTag())
        : BaseT(c0,c1,c2,c3,c4,c5,c6,c7,c8,c9)
        , ImplementationBase(backend,ordering_tag)
      {}

      template<typename Child0, typename... Children>
      DynamicPowerGridFunctionSpace(std::shared_ptr<Child0> child0, std::shared_ptr<Children>... children)
        : BaseT(child0, children...)
        , ImplementationBase(Backend(),OrderingTag())
      {}

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

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_DYNAMICPOWERGRIDFUNCTIONSPACE_HH
