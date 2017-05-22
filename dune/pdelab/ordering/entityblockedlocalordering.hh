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

#include <dune/typetree/compositenode.hh>
#include <dune/typetree/powernode.hh>
#include <dune/typetree/traversal.hh>
#include <dune/typetree/visitor.hh>
#include <dune/typetree/typetraits.hh>

#include <dune/pdelab/ordering/gridviewordering.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup Ordering
    //! \{

    template<typename ChildOrdering, std::size_t k>
    class PowerEntityBlockedLocalOrdering
      : public TypeTree::PowerNode<ChildOrdering,k>
      , public LocalOrderingBase<typename ChildOrdering::Traits::EntitySet,
                                 typename ChildOrdering::Traits::DOFIndex,
                                 typename ChildOrdering::Traits::ContainerIndex>
    {

      typedef TypeTree::PowerNode<ChildOrdering,k> NodeT;
      typedef LocalOrderingBase<typename ChildOrdering::Traits::EntitySet,
                                typename ChildOrdering::Traits::DOFIndex,
                                typename ChildOrdering::Traits::ContainerIndex> BaseT;

    public:

      static const bool consume_tree_index = true;

      typedef typename BaseT::Traits Traits;

      PowerEntityBlockedLocalOrdering(const typename NodeT::NodeStorage& child_storage, bool container_blocked)
        : NodeT(child_storage)
        , BaseT(*this,container_blocked,nullptr)
      {}

      const typename Traits::EntitySet& entitySet() const
      {
        return this->child(0).entitySet();
      }

    };

    template<typename GFS, typename Transformation>
    struct power_gfs_to_local_ordering_descriptor<GFS,Transformation,EntityBlockedOrderingTag>
    {

      static const bool recursive = true;

      template<typename TC>
      struct result
      {
        typedef PowerEntityBlockedLocalOrdering<TC,TypeTree::StaticDegree<GFS>::value> type;
        typedef std::shared_ptr<type> storage_type;
      };

      template<typename TC>
      static typename result<TC>::type transform(const GFS& gfs, const Transformation& t, const std::array<std::shared_ptr<TC>,TypeTree::StaticDegree<GFS>::value>& children)
      {
        return typename result<TC>::type(children,gfs.backend().blocked(gfs));
      }

      template<typename TC>
      static typename result<TC>::storage_type transform_storage(std::shared_ptr<const GFS> gfs, const Transformation& t, const std::array<std::shared_ptr<TC>,TypeTree::StaticDegree<GFS>::value>& children)
      {
        return std::make_shared<typename result<TC>::type>(children,gfs->backend().blocked(*gfs));
      }

    };



    template<typename GFS, typename Transformation>
    struct power_gfs_to_entityblocked_ordering_descriptor
    {

      static const bool recursive = false;

      typedef TypeTree::TransformTree<GFS,gfs_to_local_ordering<Transformation> > LocalOrderingTransformation;
      typedef typename LocalOrderingTransformation::Type LocalOrdering;

      typedef GridViewOrdering<LocalOrdering> transformed_type;

      typedef std::shared_ptr<transformed_type> transformed_storage_type;

      static transformed_type transform(const GFS& gfs, const Transformation& t)
      {
        transformed_type r(make_tuple(std::make_shared<LocalOrdering>(LocalOrderingTransformation::transform(gfs,gfs_to_local_ordering<Transformation>()))),gfs.backend().blocked(gfs),const_cast<GFS*>(&gfs));
        return std::move(r);
      }

      static transformed_storage_type transform_storage(std::shared_ptr<const GFS> gfs, const Transformation& t)
      {
        transformed_storage_type r(std::make_shared<transformed_type>(make_tuple(LocalOrderingTransformation::transform_storage(gfs,gfs_to_local_ordering<Transformation>())),gfs->backend().blocked(*gfs),const_cast<GFS*>(gfs.get())));
        return std::move(r);
      }

    };

    template<typename GFS, typename Transformation>
    power_gfs_to_entityblocked_ordering_descriptor<GFS,Transformation>
    register_power_gfs_to_ordering_descriptor(GFS*,Transformation*,EntityBlockedOrderingTag*);



    template<typename... Children>
    class CompositeEntityBlockedLocalOrdering
      : public TypeTree::CompositeNode<Children...>
      , public LocalOrderingBase<typename first_type<Children...>::type::Traits::EntitySet,
                                 typename first_type<Children...>::type::Traits::DOFIndex,
                                 typename first_type<Children...>::type::Traits::ContainerIndex>
    {

      typedef TypeTree::CompositeNode<Children...> Node;
      typedef LocalOrderingBase<typename first_type<Children...>::type::Traits::EntitySet,
                                typename first_type<Children...>::type::Traits::DOFIndex,
                                typename first_type<Children...>::type::Traits::ContainerIndex> Base;

    public:

      typedef typename Base::Traits Traits;

      static const bool consume_tree_index = true;

      CompositeEntityBlockedLocalOrdering(bool container_blocked, std::shared_ptr<Children>... children)
        : Node(children...)
        , Base(*this,container_blocked,nullptr)
      {}

      const typename Traits::EntitySet& entitySet() const
      {
        return this->template child<0>().entitySet();
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
        typedef std::shared_ptr<type> storage_type;
      };

      template<typename... TC>
      static typename result<TC...>::type transform(const GFS& gfs, const Transformation& t, std::shared_ptr<TC>... children)
      {
        return typename result<TC...>::type(gfs.backend().blocked(gfs),children...);
      }

      template<typename... TC>
      static typename result<TC...>::storage_type transform_storage(std::shared_ptr<const GFS> gfs, const Transformation& t, std::shared_ptr<TC>... children)
      {
        return std::make_shared<typename result<TC...>::type>(gfs.backend().blocked(*gfs),children...);
      }

    };

    template<typename GFS, typename Transformation>
    struct composite_gfs_to_entityblocked_ordering_descriptor
    {

      static const bool recursive = false;

      typedef TypeTree::TransformTree<GFS,gfs_to_local_ordering<Transformation> > LocalOrderingTransformation;
      typedef typename LocalOrderingTransformation::Type LocalOrdering;

      typedef GridViewOrdering<LocalOrdering> transformed_type;

      typedef std::shared_ptr<transformed_type> transformed_storage_type;

      static transformed_type transform(const GFS& gfs, const Transformation& t)
      {
        transformed_type r(make_tuple(std::make_shared<LocalOrdering>(LocalOrderingTransformation::transform(gfs,gfs_to_local_ordering<Transformation>()))),gfs.backend().blocked(gfs),const_cast<GFS*>(&gfs));
        return std::move(r);
      }

      static transformed_storage_type transform_storage(std::shared_ptr<const GFS> gfs, const Transformation& t)
      {
        transformed_storage_type r(std::make_shared<transformed_type>(make_tuple(LocalOrderingTransformation::transform_storage(gfs,gfs_to_local_ordering<Transformation>())),gfs->backend().blocked(*gfs)),const_cast<GFS*>(gfs.get()));
        return std::move(r);
      }

    };

    template<typename GFS, typename Transformation>
    composite_gfs_to_entityblocked_ordering_descriptor<GFS,Transformation>
    register_composite_gfs_to_ordering_descriptor(GFS*,Transformation*,EntityBlockedOrderingTag*);


   //! \} group Ordering
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ORDERING_ENTITYBLOCKEDLOCALORDERING_HH
