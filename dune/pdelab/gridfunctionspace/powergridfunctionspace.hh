// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_POWERGRIDFUNCTIONSPACE_HH
#define DUNE_PDELAB_POWERGRIDFUNCTIONSPACE_HH

#include <cstddef>

#include <dune/common/shared_ptr.hh>

#include <dune/pdelab/common/typetree/powernode.hh>
#include <dune/pdelab/gridfunctionspace/lexicographicordering.hh>
#include <dune/pdelab/gridfunctionspace/orderingbase.hh>
#include <dune/pdelab/gridfunctionspace/powercompositegridfunctionspacebase.hh>
#include <dune/pdelab/gridfunctionspace/tags.hh>

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
    template<typename T, std::size_t k,
             typename OrderingTag = LexicographicOrderingTag>
    class PowerGridFunctionSpace :
      public TypeTree::PowerNode<T,k>,
      public PowerCompositeGridFunctionSpaceBase<
        PowerGridFunctionSpace<T, k, OrderingTag>,
        typename T::Traits::GridViewType,
        typename T::Traits::BackendType,
        OrderingTag,
        k
      >
    {
      typedef TypeTree::PowerNode<T,k> BaseT;

      typedef PowerCompositeGridFunctionSpaceBase<
        PowerGridFunctionSpace,
        typename T::Traits::GridViewType,
        typename T::Traits::BackendType,
        OrderingTag,
        k
        > ImplementationBase;

      friend class PowerCompositeGridFunctionSpaceBase<
        PowerGridFunctionSpace,
        typename T::Traits::GridViewType,
        typename T::Traits::BackendType,
        OrderingTag,
        k>;

    public:
      typedef PowerGridFunctionSpaceTag ImplementationTag;

      typedef typename TransformPowerGFSToOrdering<OrderingTag>::
        template result<
          typename ImplementationBase::Traits,
          const typename T::Ordering,
          k
        >::type Ordering;

      //! export traits class
      typedef typename ImplementationBase::Traits Traits;

      PowerGridFunctionSpace(T& c)
        : BaseT(c)
      {
        initOrdering();
      }

      PowerGridFunctionSpace (T& c0,
                              T& c1)
        : BaseT(c0,c1)
      {
        initOrdering();
      }

      PowerGridFunctionSpace (T& c0,
                              T& c1,
                              T& c2)
        : BaseT(c0,c1,c2)
      {
        initOrdering();
      }

      PowerGridFunctionSpace (T& c0,
                              T& c1,
                              T& c2,
                              T& c3)
        : BaseT(c0,c1,c2,c3)
      {
        initOrdering();
      }

      PowerGridFunctionSpace (T& c0,
                              T& c1,
                              T& c2,
                              T& c3,
                              T& c4)
        : BaseT(c0,c1,c2,c3,c4)
      {
        initOrdering();
      }

      PowerGridFunctionSpace (T& c0,
                              T& c1,
                              T& c2,
                              T& c3,
                              T& c4,
                              T& c5)
        : BaseT(c0,c1,c2,c3,c4,c5)
      {
        initOrdering();
      }

      PowerGridFunctionSpace (T& c0,
                              T& c1,
                              T& c2,
                              T& c3,
                              T& c4,
                              T& c5,
                              T& c6)
        : BaseT(c0,c1,c2,c3,c4,c5,c6)
      {
        initOrdering();
      }

      PowerGridFunctionSpace (T& c0,
                              T& c1,
                              T& c2,
                              T& c3,
                              T& c4,
                              T& c5,
                              T& c6,
                              T& c7)
        : BaseT(c0,c1,c2,c3,c4,c5,c6,c7)
      {
        initOrdering();
      }

      PowerGridFunctionSpace (T& c0,
                              T& c1,
                              T& c2,
                              T& c3,
                              T& c4,
                              T& c5,
                              T& c6,
                              T& c7,
                              T& c8)
        : BaseT(c0,c1,c2,c3,c4,c5,c6,c7,c8)
      {
        initOrdering();
      }

      PowerGridFunctionSpace (T& c0,
                              T& c1,
                              T& c2,
                              T& c3,
                              T& c4,
                              T& c5,
                              T& c6,
                              T& c7,
                              T& c8,
                              T& c9)
        : BaseT(c0,c1,c2,c3,c4,c5,c6,c7,c8,c9)
      {
        initOrdering();
      }

      //! Direct access to the DOF ordering.
      const Ordering &ordering() const { return *orderingp; }

      //! Direct access to the storage of the DOF ordering.
      shared_ptr<const Ordering> orderingPtr() const { return orderingp; }

    private:
      void initOrdering() {
        typename Ordering::NodeStorage transformedChildren;
        for(std::size_t childIndex = 0; childIndex < BaseT::CHILDREN;
            ++childIndex)
          transformedChildren[childIndex] =
            this->child(childIndex).orderingPtr();
        orderingp = make_shared<Ordering>(*this, transformedChildren);
      }

      shared_ptr<Ordering> orderingp;
    };

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_POWERGRIDFUNCTIONSPACE_HH
