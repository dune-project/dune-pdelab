// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_VECTORGRIDFUNCTIONSPACE_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_VECTORGRIDFUNCTIONSPACE_HH

#include <algorithm>
#include <cstddef>

#include <dune/common/shared_ptr.hh>

#include <dune/typetree/powernode.hh>

#include <dune/pdelab/gridfunctionspace/powercompositegridfunctionspacebase.hh>
#include <dune/pdelab/gridfunctionspace/tags.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/powergridfunctionspace.hh>

namespace Dune {
  namespace PDELab {

    //=======================================
    // power grid function space
    //=======================================

    /** \brief base class for tuples of grid function spaces
        product of identical grid function spaces
        base class that holds implementation of the methods

        PGFS(T,k) = {T}^k

        \tparam T the underlying are all grid function spaces
        \tparam k power factor
        \tparam Mapper is the ordering parameter. Use e.g.
        \link GridFunctionSpaceLexicographicMapper GridFunctionSpaceLexicographicMapper \endlink
        or \link  GridFunctionSpaceComponentBlockwiseMapper  GridFunctionSpaceComponentBlockwiseMapper \endlink
        or \link  GridFunctionSpaceBlockwiseMapper  GridFunctionSpaceBlockwiseMapper \endlink
        or \link  GridFunctionSpaceDynamicBlockwiseMapper  GridFunctionSpaceDynamicBlockwiseMapper \endlink
    */
    template<typename GV,
             typename FEM,
             std::size_t k,
             typename Backend,
             typename LeafBackend,
             typename Constraints = NoConstraints,
             typename OrderingTag = LexicographicOrderingTag,
             typename LeafOrderingTag = DefaultLeafOrderingTag>
    class VectorGridFunctionSpace
      : public TypeTree::PowerNode<GridFunctionSpace<
                                     impl::EntitySet<GV>,
                                     FEM,
                                     Constraints,
                                     LeafBackend,
                                     LeafOrderingTag
                                     >,
                                   k>
      , public PowerCompositeGridFunctionSpaceBase<VectorGridFunctionSpace<
                                                     GV,
                                                     FEM,
                                                     k,
                                                     Backend,
                                                     LeafBackend,
                                                     Constraints,
                                                     OrderingTag,
                                                     LeafOrderingTag
                                                     >,
                                                   impl::EntitySet<GV>,
                                                   Backend,
                                                   OrderingTag,
                                                   k>

      , public DataHandleProvider<VectorGridFunctionSpace<
                                    GV,
                                    FEM,
                                    k,
                                    Backend,
                                    LeafBackend,
                                    Constraints,
                                    OrderingTag,
                                    LeafOrderingTag
                                    > >

      , public GridFunctionOutputParameters
    {

      typedef GridFunctionSpace<
        impl::EntitySet<GV>,
        FEM,
        Constraints,
        LeafBackend,
        LeafOrderingTag
        > LeafGFS;

      template<typename,typename>
      friend class GridFunctionSpaceBase;

    public:

      typedef VectorGridFunctionSpaceTag ImplementationTag;

      typedef TypeTree::PowerNode<LeafGFS,k> BaseT;

      typedef PowerCompositeGridFunctionSpaceBase<
        VectorGridFunctionSpace,
        impl::EntitySet<GV>,
        Backend,
        OrderingTag,
        k> ImplementationBase;

      friend class PowerCompositeGridFunctionSpaceBase<
        VectorGridFunctionSpace,
        impl::EntitySet<GV>,
        Backend,
        OrderingTag,
        k>;

      typedef TypeTree::TransformTree<VectorGridFunctionSpace,
                                      gfs_to_ordering<VectorGridFunctionSpace>
                                      > ordering_transformation;

    public:

      typedef typename ordering_transformation::Type Ordering;

      //! export traits class
      typedef typename ImplementationBase::Traits Traits;

    private:

      // Preconstruct children - it is important that the children are set before entering the constructor
      // of ImplementationBase!
      static typename BaseT::NodeStorage create_components(const typename Traits::EntitySet& es,
                                                           std::shared_ptr<const FEM> fem_ptr,
                                                           const LeafBackend& leaf_backend,
                                                           const LeafOrderingTag& leaf_ordering_tag)
      {
        typename BaseT::NodeStorage r;
        for (std::size_t i = 0; i < k; ++i)
          r[i] = std::make_shared<LeafGFS>(es,fem_ptr,leaf_backend,leaf_ordering_tag);
        return r;
      }

    public:

      VectorGridFunctionSpace(const typename Traits::GridView& gv, const FEM& fem,
                              const Backend& backend = Backend(), const LeafBackend& leaf_backend = LeafBackend(),
                              const OrderingTag& ordering_tag = OrderingTag(), const LeafOrderingTag& leaf_ordering_tag = LeafOrderingTag())
        : BaseT(create_components(typename Traits::EntitySet(gv),stackobject_to_shared_ptr(fem),leaf_backend,leaf_ordering_tag))
        , ImplementationBase(backend,ordering_tag)
      {}

      VectorGridFunctionSpace(const typename Traits::EntitySet& es, const FEM& fem,
                              const Backend& backend = Backend(), const LeafBackend& leaf_backend = LeafBackend(),
                              const OrderingTag& ordering_tag = OrderingTag(), const LeafOrderingTag& leaf_ordering_tag = LeafOrderingTag())
        : BaseT(create_components(es,stackobject_to_shared_ptr(fem),leaf_backend,leaf_ordering_tag))
        , ImplementationBase(backend,ordering_tag)
      {}

      std::string name() const
      {
        return ImplementationBase::name();
      }

      void name(std::string name)
      {
        ImplementationBase::name(name);
        for (std::size_t i = 0; i < k; ++i)
          {
            std::stringstream ns;
            ns << name << "_" << i;
            this->child(i).name(ns.str());
          }
      }

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
            _ordering->update();
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
            _ordering->update();
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

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_VECTORGRIDFUNCTIONSPACE_HH
