// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_ORDERING_PERMUTEDORDERING_HH
#define DUNE_PDELAB_ORDERING_PERMUTEDORDERING_HH

#include <dune/typetree/typetree.hh>

#include <dune/pdelab/ordering/utility.hh>
#include <dune/pdelab/ordering/orderingbase.hh>
#include <dune/pdelab/ordering/decorator.hh>

namespace Dune {
  namespace PDELab {

    namespace ordering {

#ifndef DOXYGEN // implementation internals

      namespace permuted {

        struct tag_base
        {

          std::vector<std::size_t>& permutation()
          {
            return _permutation;
          }

          const std::vector<std::size_t>& permutation() const
          {
            return _permutation;
          }

        private:

          std::vector<std::size_t> _permutation;

        };

        template<std::size_t i>
        struct base_holder
          : public tag_base
        {};

      } // namespace permuted

#endif // DOXYGEN

      //! Permute the ordering created from the passed-in tag based on a simple lookup table.
      /**
       * This tag modifies the Ordering designed by the passed-in OrderingTag to perform an additional
       * permutation step on the ContainerIndex entries touched by that Ordering. In layman's terms:
       * You can reorder the list of ContainerIndices you would get if the GridFunctionSpace associated
       * with this OrderingTag were the topmost space, but you cannot perform any reordering within those
       * individual blocks.
       *
       *
       *
       * \tparam OrderingTag  The tag describing the Ordering that will be permuted.
       */
      template<typename OrderingTag>
      struct Permuted
        : public permuted::base_holder<decorated_ordering_tag<Permuted<OrderingTag>,OrderingTag>::level>
        , public decorated_ordering_tag<Permuted<OrderingTag>,OrderingTag>
      {

        Permuted()
        {}

        Permuted(const OrderingTag& tag)
          : decorated_ordering_tag<Permuted<OrderingTag>,OrderingTag>(tag)
        {}

        Permuted(OrderingTag&& tag)
          : decorated_ordering_tag<Permuted<OrderingTag>,OrderingTag>(std::move(tag))
        {}

        template<std::size_t i>
        const permuted::base_holder<i>& permuted() const
        {
          return *this;
        }

        template<std::size_t i>
        permuted::base_holder<i>& permuted()
        {
          return *this;
        }

      };

    } // namespace ordering

    //! \addtogroup Ordering
    //! \{

    //! Ordering that permutes top-level ContainerIndex entries.
    template<typename Ordering>
    class PermutedOrdering
      : public TypeTree::CompositeNode<Ordering>
      , public VirtualOrderingBase<typename Ordering::Traits::DOFIndex,
                                   typename Ordering::Traits::ContainerIndex>
      , public OrderingBase<typename Ordering::Traits::DOFIndex,
                            typename Ordering::Traits::ContainerIndex>
    {
    public:
      typedef typename Ordering::Traits Traits;

      static const bool has_dynamic_ordering_children = true;

      static const bool consume_tree_index = false;

    private:

      typedef TypeTree::CompositeNode<Ordering> NodeT;

      typedef OrderingBase<typename Ordering::Traits::DOFIndex,
                           typename Ordering::Traits::ContainerIndex> BaseT;

    public:

      Ordering& ordering()
      {
        return this->template child<0>();
      }

      const Ordering& ordering() const
      {
        return this->template child<0>();
      }


      PermutedOrdering(const typename NodeT::NodeStorage& ordering, const ordering::permuted::tag_base& tag)
        : NodeT(ordering)
        , BaseT(*this,false,nullptr,this)
        , _tag(tag)
      {}

      PermutedOrdering(const PermutedOrdering& r)
        : NodeT(r.nodeStorage())
        , BaseT(r)
        ,_tag(r._tag)
      {
        this->setDelegate(this);
      }

      PermutedOrdering(PermutedOrdering&& r)
        : NodeT(r.nodeStorage())
        , BaseT(std::move(r))
        ,_tag(r._tag)
      {
        this->setDelegate(this);
      }

      virtual void map_index_dynamic(typename Traits::DOFIndexView di, typename Traits::ContainerIndex& ci) const
      {
        ordering().mapIndex(di,ci);
        ci.back() = _tag.permutation()[ci.back()];
      }

      template<typename ItIn, typename ItOut>
      void map_lfs_indices(ItIn in, const ItIn end, ItOut out) const
      {
        for (; in != end; ++in, ++out)
          {
            out->back() = _tag.permutation()[out->back()];
          }
      }

      template<typename CIOutIterator>
      typename Traits::SizeType
      extract_entity_indices(const typename Traits::DOFIndex::EntityIndex& ei,
                             typename Traits::SizeType child_index,
                             CIOutIterator ci_out, const CIOutIterator ci_end) const
      {
        for (; ci_out != ci_end; ++ci_out)
          {
            ci_out->back() = _tag.permutation()[ci_out->back()];
          }
        return 0;
      }

      void update()
      {
        ordering().update();
        BaseT::update();
        if (!_tag.permutation().empty() && _tag.permutation().size() != this->blockCount())
          DUNE_THROW(PermutedOrderingSizeError,
                     "Size of permutation array does not match block count of ordering: "
                     << _tag.permutation().size()
                     << " != "
                     << this->blockCount()
                     );
        else
          {
            auto& mutable_tag = const_cast<ordering::permuted::tag_base&>(_tag);
            mutable_tag.permutation().resize(this->blockCount());
            std::iota(
              mutable_tag.permutation().begin(),
              mutable_tag.permutation().end(),
              0
              );
          }
      }

    private:

      const ordering::permuted::tag_base& _tag;

    };

    namespace ordering {

      namespace permuted {

        template<typename GFS, typename Transformation, typename Undecorated, typename Tag>
        struct gfs_to_permuted
        {

          typedef PermutedOrdering<Undecorated> transformed_type;
          typedef std::shared_ptr<transformed_type> transformed_storage_type;

          static transformed_type transform(const GFS& gfs, const Transformation& t, std::shared_ptr<Undecorated> undecorated)
          {
            return transformed_type(make_tuple(undecorated),gfs.orderingTag().template permuted<Tag::level>());
          }

          static transformed_storage_type transform_storage(std::shared_ptr<const GFS> gfs_pointer, const Transformation& t, std::shared_ptr<Undecorated> undecorated)
          {
            return std::make_shared<transformed_type>(make_tuple(undecorated),gfs_pointer->orderingTag().template permuted<Tag::level>());
          }

        };

        template<typename GFS, typename Transformation, typename Undecorated, typename GlueTag, typename UndecoratedTag>
        gfs_to_permuted<GFS,Transformation,Undecorated,GlueTag>
        register_gfs_to_decorator_descriptor(GFS*,Transformation*,Undecorated*,GlueTag*,Permuted<UndecoratedTag>*);

      } // namespace permuted
    } // namespace ordering


    template<typename GFS, typename Transformation, typename U>
    struct power_gfs_to_local_ordering_descriptor<GFS,Transformation,ordering::Permuted<U> >
      : public power_gfs_to_local_ordering_descriptor<GFS,Transformation,U>
    {};


    template<typename GFS, typename Transformation, typename U>
    struct composite_gfs_to_local_ordering_descriptor<GFS,Transformation,ordering::Permuted<U> >
      : public composite_gfs_to_local_ordering_descriptor<GFS,Transformation,U>
    {};

    //! \} group Ordering
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ORDERING_PERMUTEDORDERING_HH
