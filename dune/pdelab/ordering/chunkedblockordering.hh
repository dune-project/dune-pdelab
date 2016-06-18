// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_ORDERING_CHUNKEDBLOCKORDERING_HH
#define DUNE_PDELAB_ORDERING_CHUNKEDBLOCKORDERING_HH

#include <dune/typetree/typetree.hh>

#include <dune/pdelab/ordering/utility.hh>
#include <dune/pdelab/ordering/orderingbase.hh>
#include <dune/pdelab/ordering/decorator.hh>

namespace Dune {
  namespace PDELab {

    namespace ordering {

#ifndef DOXYGEN // implementation internals

      namespace chunked {

        struct tag_base
        {

          tag_base(std::size_t block_size)
            : _block_size(block_size)
          {}

          std::size_t blockSize() const
          {
            return _block_size;
          }

        private:

          const std::size_t _block_size;

        };

        template<std::size_t i>
        struct base_holder
          : public tag_base
        {

          base_holder(std::size_t block_size)
            : tag_base(block_size)
          {}

        };

      } // namespace chunked

#endif // DOXYGEN

      //! Block the ordering created from the passed-in tag by simply carving it up into blocks of the chunk size
      //! passed to this ordering tag.
      /**
       * This tag modifies the Ordering designed by the passed-in OrderingTag to perform an additional
       * blocking step on the ContainerIndex entries touched by that Ordering. In layman's terms:
       * You can block the list of ContainerIndices by a fixed chunk block size, as long as the block count of the
       * underlying node is the multiple of the chunk size. Otherwise, an exception will be raised.
       *
       * \tparam OrderingTag  The tag describing the Ordering that will be blocked.
       */
      template<typename OrderingTag>
      struct Chunked
        : public chunked::base_holder<decorated_ordering_tag<Chunked<OrderingTag>,OrderingTag>::level>
        , public decorated_ordering_tag<Chunked<OrderingTag>,OrderingTag>
      {

        Chunked(std::size_t block_size)
          : chunked::base_holder<decorated_ordering_tag<Chunked<OrderingTag>,OrderingTag>::level>(block_size)
        {}

        Chunked(std::size_t block_size, const OrderingTag& tag)
          : chunked::base_holder<decorated_ordering_tag<Chunked<OrderingTag>,OrderingTag>::level>(block_size)
          , decorated_ordering_tag<Chunked<OrderingTag>,OrderingTag>(tag)
        {}

        Chunked(std::size_t block_size, OrderingTag&& tag)
          : chunked::base_holder<decorated_ordering_tag<Chunked<OrderingTag>,OrderingTag>::level>(block_size)
          , decorated_ordering_tag<Chunked<OrderingTag>,OrderingTag>(std::move(tag))
        {}

        template<std::size_t i>
        const chunked::base_holder<i>& chunked() const
        {
          return *this;
        }

        template<std::size_t i>
        chunked::base_holder<i>& chunked()
        {
          return *this;
        }

      };

      template<typename Decorated>
      constexpr bool deactivate_standard_blocking_for_ordering(const Chunked<Decorated>&)
      {
        return true;
      }

    } // namespace ordering

    //! \addtogroup Ordering
    //! \{

    //! Ordering that permutes top-level ContainerIndex entries.
    template<typename Ordering>
    class ChunkedBlockOrdering
      : public TypeTree::CompositeNode<Ordering>
      , public VirtualOrderingBase<typename Ordering::Traits::DOFIndex,
                                   typename Ordering::Traits::ContainerIndex>
      , public OrderingBase<typename Ordering::Traits::DOFIndex,
                            typename Ordering::Traits::ContainerIndex>
    {
    public:
      typedef typename Ordering::Traits Traits;

      static const bool consume_tree_index = false;

    private:

      typedef TypeTree::CompositeNode<Ordering> NodeT;

      typedef OrderingBase<
        typename Ordering::Traits::DOFIndex,
        typename Ordering::Traits::ContainerIndex
        > BaseT;

    public:

      Ordering& ordering()
      {
        return this->template child<0>();
      }

      const Ordering& ordering() const
      {
        return this->template child<0>();
      }


      ChunkedBlockOrdering(const typename NodeT::NodeStorage& ordering, const ordering::chunked::tag_base& tag)
        : NodeT(ordering)
        , BaseT(*this,true,nullptr,this)
        , _tag(tag)
      {}

      ChunkedBlockOrdering(const ChunkedBlockOrdering& r)
        : NodeT(r.nodeStorage())
        , BaseT(r)
        ,_tag(r._tag)
      {
        this->setDelegate(this);
      }

      ChunkedBlockOrdering(ChunkedBlockOrdering&& r)
        : NodeT(r.nodeStorage())
        , BaseT(std::move(r))
        ,_tag(r._tag)
      {
        this->setDelegate(this);
      }

      virtual void map_index_dynamic(typename Traits::DOFIndexView di, typename Traits::ContainerIndex& ci) const
      {
        ordering().mapIndex(di,ci);
        std::size_t flat_index = ci.back();
        std::size_t block_index = flat_index / _tag.blockSize();
        std::size_t inner_index = flat_index % _tag.blockSize();
        ci.back() = inner_index;
        ci.push_back(block_index);
      }

      template<typename ItIn, typename ItOut>
      void map_lfs_indices(ItIn in, const ItIn end, ItOut out) const
      {
        for (; in != end; ++in, ++out)
          {
            std::size_t flat_index = out->back();
            std::size_t block_index = flat_index / _tag.blockSize();
            std::size_t inner_index = flat_index % _tag.blockSize();
            out->back() = inner_index;
            out->push_back(block_index);
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
            std::size_t flat_index = ci_out->back();
            std::size_t block_index = flat_index / _tag.blockSize();
            std::size_t inner_index = flat_index % _tag.blockSize();
            ci_out->back() = inner_index;
            ci_out->push_back(block_index);
          }
        return 0;
      }

      void update()
      {
        ordering().update();
        BaseT::update();
        if (ordering().blockCount() % _tag.blockSize() != 0)
          DUNE_THROW(ChunkedBlockOrderingSizeError,
                     "Block size of chunked block ordering does not divide the block count "
                     "of the underlying ordering: "
                     << ordering().blockCount()
                     << " % "
                     << _tag.blockSize()
                     << " != 0"
                     );
        this->_block_count = ordering().blockCount() / _tag.blockSize();
      }

    private:

      const ordering::chunked::tag_base& _tag;

    };

    namespace ordering {

      namespace chunked {

        template<typename GFS, typename Transformation, typename Undecorated, typename Tag>
        struct gfs_to_chunked
        {

          typedef ChunkedBlockOrdering<Undecorated> transformed_type;
          typedef std::shared_ptr<transformed_type> transformed_storage_type;

          static transformed_type transform(const GFS& gfs, const Transformation& t, std::shared_ptr<Undecorated> undecorated)
          {
            return transformed_type(make_tuple(undecorated),gfs.orderingTag().template chunked<Tag::level>());
          }

          static transformed_storage_type transform_storage(std::shared_ptr<const GFS> gfs_pointer, const Transformation& t, std::shared_ptr<Undecorated> undecorated)
          {
            return std::make_shared<transformed_type>(make_tuple(undecorated),gfs_pointer->orderingTag().template chunked<Tag::level>());
          }

        };

        template<typename GFS, typename Transformation, typename Undecorated, typename GlueTag, typename UndecoratedTag>
        gfs_to_chunked<GFS,Transformation,Undecorated,GlueTag>
        register_gfs_to_decorator_descriptor(GFS*,Transformation*,Undecorated*,GlueTag*,Chunked<UndecoratedTag>*);

      } // namespace chunked
    } // namespace ordering


    template<typename GFS, typename Transformation, typename U>
    struct power_gfs_to_local_ordering_descriptor<GFS,Transformation,ordering::Chunked<U> >
      : public power_gfs_to_local_ordering_descriptor<GFS,Transformation,U>
    {};


    template<typename GFS, typename Transformation, typename U>
    struct composite_gfs_to_local_ordering_descriptor<GFS,Transformation,ordering::Chunked<U> >
      : public composite_gfs_to_local_ordering_descriptor<GFS,Transformation,U>
    {};

    //! \} group Ordering
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ORDERING_CHUNKEDBLOCKORDERING_HH
