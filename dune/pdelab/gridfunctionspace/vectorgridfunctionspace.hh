// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_VECTORGRIDFUNCTIONSPACE_HH
#define DUNE_PDELAB_VECTORGRIDFUNCTIONSPACE_HH

#include <cstddef>

#include <dune/common/shared_ptr.hh>

#include <dune/pdelab/common/typetree/powernode.hh>
#include <dune/pdelab/gridfunctionspace/lexicographicordering.hh>
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
                                     GV,
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
                                                   GV,
                                                   Backend,
                                                   OrderingTag,
                                                   k>
      , public GridFunctionOutputParameters
    {

      typedef GridFunctionSpace<
        GV,
        FEM,
        Constraints,
        LeafBackend,
        LeafOrderingTag
        > LeafGFS;

    public:

      typedef VectorGridFunctionSpaceTag ImplementationTag;

      typedef TypeTree::PowerNode<LeafGFS,k> BaseT;

      typedef PowerCompositeGridFunctionSpaceBase<
        VectorGridFunctionSpace,
        GV,
        Backend,
        OrderingTag,
        k> ImplementationBase;

      friend class PowerCompositeGridFunctionSpaceBase<
        VectorGridFunctionSpace,
        GV,
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

      VectorGridFunctionSpace(const GV& gv, const FEM& fem, const Backend& backend = Backend(), const LeafBackend& leafBackend = LeafBackend())
        : ImplementationBase(backend)
      {
        shared_ptr<const FEM> fem_ptr = stackobject_to_shared_ptr(fem);
        for (std::size_t i = 0; i < k; ++i)
          this->setChild(i,make_shared<LeafGFS>(gv,fem_ptr,leafBackend));
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

      std::string name() const
      {
        return ImplementationBase::name();
      }

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

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_VECTORGRIDFUNCTIONSPACE_HH
