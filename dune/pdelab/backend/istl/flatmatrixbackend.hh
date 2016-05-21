// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_BACKEND_ISTL_FLATMATRIXBACKEND_HH
#define DUNE_PDELAB_BACKEND_ISTL_FLATMATRIXBACKEND_HH

#include <dune/common/iteratoradapters.hh>

#include <dune/istl/ellmatrix/host.hh>

#include <dune/pdelab/backend/common/tags.hh>
#include <dune/pdelab/backend/common/uncachedmatrixview.hh>
#include <dune/pdelab/backend/istl/matrixhelpers.hh>
#include <dune/pdelab/backend/istl/descriptors.hh>
#include <dune/pdelab/backend/istl/flat/pattern.hh>
#include <dune/pdelab/backend/interface.hh>

namespace Dune {
  namespace PDELab {
    namespace istl {


    template<typename GFSV, typename GFSU, typename C>
    class FlatELLMatrixContainer
      : public Backend::impl::Wrapper<C>
    {

      friend Backend::impl::Wrapper<C>;

    public:

      typedef typename C::value_type ElementType;
      typedef ElementType E;
      typedef C Container;
      typedef C BaseT;
      typedef typename C::value_type field_type;
      typedef typename C::size_type size_type;

      typedef GFSU TrialGridFunctionSpace;
      typedef GFSV TestGridFunctionSpace;

      typedef typename GFSV::Ordering::Traits::ContainerIndex RowIndex;
      typedef typename GFSU::Ordering::Traits::ContainerIndex ColIndex;

      typedef flat::Pattern<typename GFSV::Ordering,typename GFSU::Ordering> Pattern;

    private:

      static const size_type kernel_block_size = Container::kernel_block_size;
      static const size_type block_shift = Container::block_shift;
      static const size_type block_mask = Container::block_mask;

    public:

      template<typename RowCache, typename ColCache>
      using LocalView = UncachedMatrixView<FlatELLMatrixContainer,RowCache,ColCache>;

      template<typename RowCache, typename ColCache>
      using ConstLocalView = ConstUncachedMatrixView<const FlatELLMatrixContainer,RowCache,ColCache>;

    private:

      static typename Container::Layout buildLayout(const Pattern& pattern) {
        const size_type rows = pattern.rowOrdering().blockCount();
        const size_type cols = pattern.colOrdering().blockCount();
        typename Container::LayoutBuilder layout;
        layout.setSize(rows,cols);
        layout.allocateRows();
        pattern.sizes(layout.rowLength());
        layout.allocateCols();
        auto index_streamer = pattern.indexStreamer();
        typedef StridedIterator<size_type*,kernel_block_size> It;
        size_type* col_index = layout.colIndex();
        for (size_type b = 0, blocks = layout.blocks(); b < blocks; ++b)
          {
            size_type block_length = layout.blockLength(b);
            // use signed indices here to avoid underflows in nonpadded_block_size!
            typedef typename Container::LayoutBuilder::Allocator::difference_type diff_t;
            diff_t nonpadded_block_size = std::max(std::min(diff_t(kernel_block_size),diff_t(rows) - diff_t(b * kernel_block_size)),diff_t(0));
            for (diff_t i = 0; i < nonpadded_block_size; ++i)
              {
                It col_begin(col_index+i);
                It col_end(col_begin + block_length);
                It row_end = index_streamer.streamRow(col_begin);
                std::fill(row_end,col_end,(row_end != col_begin ? row_end[-1] : size_type(0)));
              }
            for (diff_t i = nonpadded_block_size; i < diff_t(kernel_block_size); ++i)
              {
                It col_begin(col_index+i);
                It col_end(col_begin + block_length);
                std::fill(col_begin,col_end,size_type(0));
              }
            col_index += kernel_block_size * block_length;
          }
        return layout.layout();
      }

    public:

      template<typename GO>
      explicit FlatELLMatrixContainer (const GO& go)
        : _container(std::make_shared<Container>())
      {
        Pattern pattern(
          go.testGridFunctionSpace().ordering(),
          go.trialGridFunctionSpace().ordering(),
          go.matrixBackend().entriesPerRow()
          );
        go.fill_pattern(pattern);
        _container->setLayout(buildLayout(pattern));
      }

      template<typename GO>
      FlatELLMatrixContainer (const GO& go, const E& e)
        : _container(std::make_shared<Container>())
      {
        Pattern pattern(
          go.testGridFunctionSpace().ordering(),
          go.trialGridFunctionSpace().ordering(),
          go.matrixBackend().entriesPerRow()
          );
        go.fill_pattern(pattern);
        _container->setLayout(buildLayout(pattern));
        (*_container) = e;
      }


      //! Creates an FlatELLMatrixContainer without allocating an underlying ISTL matrix.
      explicit FlatELLMatrixContainer (Dune::PDELab::Backend::unattached_container = Dune::PDELab::Backend::unattached_container())
      {}

      //! Creates an FlatELLMatrixContainer with an empty underlying ISTL matrix.
      explicit FlatELLMatrixContainer (Dune::PDELab::Backend::attached_container)
        : _container(std::make_shared<Container>())
      {}

      FlatELLMatrixContainer(const FlatELLMatrixContainer& rhs)
        : _container(std::make_shared<Container>(*(rhs._container)))
      {}

      FlatELLMatrixContainer& operator=(const FlatELLMatrixContainer& rhs)
      {
        if (this == &rhs)
          return *this;
        if (attached())
          {
            (*_container) = (*(rhs._container));
          }
        else
          {
            _container = std::make_shared<Container>(*(rhs._container));
          }
        return *this;
      }

      void detach()
      {
        _container.reset();
      }

      void attach(std::shared_ptr<Container> container)
      {
        _container = container;
      }

      bool attached() const
      {
        return bool(_container);
      }

      const std::shared_ptr<Container>& storage() const
      {
        return _container;
      }

      size_type N() const
      {
        return _container->N();
      }

      size_type M() const
      {
        return _container->M();
      }

      FlatELLMatrixContainer& operator= (const E& e)
      {
        (*_container) = e;
        return *this;
      }

      FlatELLMatrixContainer& operator*= (const E& e)
      {
        (*_container) *= e;
        return *this;
      }

      E& operator()(const RowIndex& ri, const ColIndex& ci)
      {
        assert(ri.size() == 1);
        assert(ci.size() == 1);
        return (*_container)(ri[0],ci[0]);
      }

      const E& operator()(const RowIndex& ri, const ColIndex& ci) const
      {
        assert(ri.size() == 1);
        assert(ci.size() == 1);
        return (*_container)(ri[0],ci[0]);
      }

      void flush()
      {}

      void finalize()
      {}

      void clear_row(const RowIndex& ri, const E& diagonal_entry)
      {
        assert(ri.size() == 1);
        _container->clearRow(ri[0]);
        (*this)(ri,ri) = diagonal_entry;
      }

    private:

      Container& native()
      {
        return *_container;
      }

      const Container& native() const
      {
        return *_container;
      }

      std::shared_ptr<Container> _container;

    };

    } // namespace istl
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_BACKEND_ISTL_FLATMATRIXBACKEND_HH
