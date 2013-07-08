// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_ORDERING_ENTITYBLOCKEDLOCALORDERING_HH
#define DUNE_PDELAB_ORDERING_ENTITYBLOCKEDLOCALORDERING_HH

#include <cstddef>
#include <ostream>
#include <string>

#include <dune/common/classname.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/stdstreams.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/deprecated.hh>

#include <dune/pdelab/common/typetree/compositenodemacros.hh>
#include <dune/pdelab/common/typetree/powernode.hh>
#include <dune/pdelab/common/typetree/traversal.hh>
#include <dune/pdelab/common/typetree/visitor.hh>

#include <dune/pdelab/ordering/gridviewordering.hh>

namespace Dune {
  namespace PDELab {

    //! Interface for merging index spaces
    template<typename ChildOrdering, std::size_t k>
    class PowerEntityBlockedLocalOrdering
      : public TypeTree::PowerNode<ChildOrdering,k>
      , public LocalOrderingBase<typename ChildOrdering::Traits::GridView,
                                 typename ChildOrdering::Traits::DOFIndex,
                                 typename ChildOrdering::Traits::ContainerIndex>
    {

      typedef TypeTree::PowerNode<ChildOrdering,k> NodeT;
      typedef LocalOrderingBase<typename ChildOrdering::Traits::GridView,
                                typename ChildOrdering::Traits::DOFIndex,
                                typename ChildOrdering::Traits::ContainerIndex> BaseT;

    public:

      static const bool consume_tree_index = true;

      typedef typename BaseT::Traits Traits;

      PowerEntityBlockedLocalOrdering(const typename NodeT::NodeStorage& child_storage, bool container_blocked)
        : NodeT(child_storage)
        , BaseT(*this,container_blocked,nullptr)
      {}

      const typename Traits::GridView& gridView() const
      {
        return this->child(0).gridView();
      }

    };

    template<typename GFS, typename Transformation>
    struct power_gfs_to_local_ordering_descriptor<GFS,Transformation,EntityBlockedOrderingTag>
    {

      static const bool recursive = true;

      template<typename TC>
      struct result
      {
        typedef PowerEntityBlockedLocalOrdering<TC,GFS::CHILDREN> type;
        typedef shared_ptr<type> storage_type;
      };

      template<typename TC>
      static typename result<TC>::type transform(const GFS& gfs, const Transformation& t, const array<shared_ptr<TC>,GFS::CHILDREN>& children)
      {
        return typename result<TC>::type(children,gfs.backend().blocked());
      }

      template<typename TC>
      static typename result<TC>::storage_type transform_storage(shared_ptr<const GFS> gfs, const Transformation& t, const array<shared_ptr<TC>,GFS::CHILDREN>& children)
      {
        return make_shared<typename result<TC>::type>(children,gfs->backend().blocked());
      }

    };



    template<typename GFS, typename Transformation>
    struct power_gfs_to_ordering_descriptor<GFS,Transformation,EntityBlockedOrderingTag>
    {

      static const bool recursive = false;

      typedef TypeTree::TransformTree<GFS,gfs_to_local_ordering<Transformation> > LocalOrderingTransformation;
      typedef typename LocalOrderingTransformation::Type LocalOrdering;

      typedef GridViewOrdering<LocalOrdering> transformed_type;

      typedef shared_ptr<transformed_type> transformed_storage_type;

      static transformed_type transform(const GFS& gfs, const Transformation& t)
      {
        transformed_type r(make_tuple(make_shared<LocalOrdering>(LocalOrderingTransformation::transform(gfs,gfs_to_local_ordering<Transformation>()))),gfs.backend().blocked(),const_cast<GFS*>(&gfs));
        return std::move(r);
      }

      static transformed_storage_type transform_storage(shared_ptr<const GFS> gfs, const Transformation& t)
      {
        transformed_storage_type r(make_shared<transformed_type>(make_tuple(LocalOrderingTransformation::transform_storage(gfs,gfs_to_local_ordering<Transformation>())),gfs->backend().blocked(),const_cast<GFS*>(gfs.get())));
        return std::move(r);
      }

    };

    // the generic registration for PowerGridFunctionSpace happens in transformations.hh


    //! Interface for merging index spaces
    template<DUNE_TYPETREE_COMPOSITENODE_TEMPLATE_CHILDREN>
    class CompositeEntityBlockedLocalOrdering
      : public DUNE_TYPETREE_COMPOSITENODE_BASETYPE
      , public LocalOrderingBase<typename DUNE_TYPETREE_COMPOSITENODE_FIRST_CHILD::Traits::GridView,
                                 typename DUNE_TYPETREE_COMPOSITENODE_FIRST_CHILD::Traits::DOFIndex,
                                 typename DUNE_TYPETREE_COMPOSITENODE_FIRST_CHILD::Traits::ContainerIndex>
    {

      typedef DUNE_TYPETREE_COMPOSITENODE_BASETYPE Node;
      typedef LocalOrderingBase<typename DUNE_TYPETREE_COMPOSITENODE_FIRST_CHILD::Traits::GridView,
                                typename DUNE_TYPETREE_COMPOSITENODE_FIRST_CHILD::Traits::DOFIndex,
                                typename DUNE_TYPETREE_COMPOSITENODE_FIRST_CHILD::Traits::ContainerIndex> Base;

    public:

      typedef typename Base::Traits Traits;

      static const bool consume_tree_index = true;

      CompositeEntityBlockedLocalOrdering(bool container_blocked, DUNE_TYPETREE_COMPOSITENODE_STORAGE_CONSTRUCTOR_SIGNATURE)
        : Node(DUNE_TYPETREE_COMPOSITENODE_CHILDVARIABLES)
        , Base(*this,container_blocked,nullptr)
      {}

      const typename Traits::GridView& gridView() const
      {
        return this->template child<0>().gridView();
      }

    };


    template<typename GFS, typename Transformation>
    struct composite_gfs_to_local_ordering_descriptor<GFS,Transformation,EntityBlockedOrderingTag>
    {

      static const bool recursive = true;

      template<typename... TC>
      struct result
      {
        typedef CompositeEntityBlockedLocalOrdering<TC...> type;
        typedef shared_ptr<type> storage_type;
      };

      template<typename... TC>
      static typename result<TC...>::type transform(const GFS& gfs, const Transformation& t, shared_ptr<TC>... children)
      {
        return typename result<TC...>::type(gfs.backend().blocked(),children...);
      }

      template<typename... TC>
      static typename result<TC...>::storage_type transform_storage(shared_ptr<const GFS> gfs, const Transformation& t, shared_ptr<TC>... children)
      {
        return make_shared<typename result<TC...>::type>(gfs.backend().blocked(),children...);
      }

    };

    template<typename GFS, typename Transformation>
    struct composite_gfs_to_ordering_descriptor<GFS,Transformation,EntityBlockedOrderingTag>
    {

      static const bool recursive = false;

      typedef TypeTree::TransformTree<GFS,gfs_to_local_ordering<Transformation> > LocalOrderingTransformation;
      typedef typename LocalOrderingTransformation::Type LocalOrdering;

      typedef GridViewOrdering<LocalOrdering> transformed_type;

      typedef shared_ptr<transformed_type> transformed_storage_type;

      static transformed_type transform(const GFS& gfs, const Transformation& t)
      {
        transformed_type r(make_tuple(make_shared<LocalOrdering>(LocalOrderingTransformation::transform(gfs,gfs_to_local_ordering<Transformation>()))),gfs.backend().blocked(),const_cast<GFS*>(&gfs));
        return std::move(r);
      }

      static transformed_storage_type transform_storage(shared_ptr<const GFS> gfs, const Transformation& t)
      {
        transformed_storage_type r(make_shared<transformed_type>(make_tuple(LocalOrderingTransformation::transform_storage(gfs,gfs_to_local_ordering<Transformation>())),gfs->backend().blocked()),const_cast<GFS*>(gfs.get()));
        return std::move(r);
      }

    };

   // the generic registration for CompositeGridFunctionSpace happens in transformations.hh


   //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ORDERING_ENTITYBLOCKEDORDERING_HH
