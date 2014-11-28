// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_BACKEND_ISTL_BELLMATRIXBACKEND_HH
#define DUNE_PDELAB_BACKEND_ISTL_BELLMATRIXBACKEND_HH

#include <dune/common/iteratoradapters.hh>

#include <dune/istl/bellmatrix/host.hh>

#include <dune/pdelab/backend/tags.hh>
#include <dune/pdelab/backend/common/uncachedmatrixview.hh>
#include <dune/pdelab/backend/istl/matrixhelpers.hh>
#include <dune/pdelab/backend/istl/descriptors.hh>
#include <dune/pdelab/backend/istl/flat/pattern.hh>

namespace Dune {
  namespace PDELab {
    namespace istl {


    template<typename GFSV, typename GFSU, typename C>
    class BELLMatrixContainer
    {

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
      static const size_type kernel_block_shift = Container::kernel_block_shift;
      static const size_type kernel_block_mask = Container::kernel_block_mask;

    public:

#if HAVE_TEMPLATE_ALIASES

      template<typename RowCache, typename ColCache>
      using LocalView = UncachedMatrixView<BELLMatrixContainer,RowCache,ColCache>;

      template<typename RowCache, typename ColCache>
      using ConstLocalView = ConstUncachedMatrixView<const BELLMatrixContainer,RowCache,ColCache>;

#else

      template<typename RowCache, typename ColCache>
      struct LocalView
        : public UncachedMatrixView<FlatELLMatrixContainer,RowCache,ColCache>
      {

        LocalView()
        {}

        LocalView(FlatELLMatrixContainer& mc)
          : UncachedMatrixView<FlatELLMatrixContainer,RowCache,ColCache>(mc)
        {}

      };

      template<typename RowCache, typename ColCache>
      struct ConstLocalView
        : public ConstUncachedMatrixView<const FlatELLMatrixContainer,RowCache,ColCache>
      {

        ConstLocalView()
        {}

        ConstLocalView(const FlatELLMatrixContainer& mc)
          : ConstUncachedMatrixView<const FlatELLMatrixContainer,RowCache,ColCache>(mc)
        {}

      };

#endif // HAVE_TEMPLATE_ALIASES

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
            // use signed indices here to avoid underflows in nonpadded_kernel_block_size!
            typedef typename Container::LayoutBuilder::Allocator::difference_type diff_t;
            diff_t nonpadded_kernel_block_size = std::max(std::min(diff_t(kernel_block_size),diff_t(rows) - diff_t(b * kernel_block_size)),diff_t(0));
            for (diff_t i = 0; i < nonpadded_kernel_block_size; ++i)
              {
                It col_begin(col_index+i);
                It col_end(col_begin + block_length);
                It row_end = index_streamer.streamRow(col_begin);
                std::fill(row_end,col_end,(row_end != col_begin ? row_end[-1] : size_type(0)));
              }
            // add a single column entry on the diagonal to make diagonal extraction work
            for (diff_t i = nonpadded_kernel_block_size; i < kernel_block_size; ++i)
              {
                layout.rowLength()[b * kernel_block_size + i] = 1;
                It col_begin(col_index+i);
                It col_end(col_begin + block_length);
                std::fill(col_begin,col_end,size_type(b * kernel_block_size + i));
              }
            col_index += kernel_block_size * block_length;
          }
        return layout.layout();
      }

    public:

      template<typename GO, typename Parameters>
      explicit BELLMatrixContainer (const GO& go, Parameters parameters)
        : _container(make_shared<Container>(go.testGridFunctionSpace().backend().blockSize(),go.trialGridFunctionSpace().backend().blockSize()))
      {
        Pattern pattern(
          go.testGridFunctionSpace().ordering(),
          go.trialGridFunctionSpace().ordering(),
          parameters.entriesPerRow()
          );
        go.fill_pattern(pattern);
        _container->setLayout(buildLayout(pattern));
      }

      template<typename GO, typename Parameters>
      BELLMatrixContainer (const GO& go, Parameters parameters, const E& e)
        : _container(make_shared<Container>(go.testGridFunctionSpace().backend().blockSize(),go.trialGridFunctionSpace().backend().blockSize()))
      {
        Pattern pattern(
          go.testGridFunctionSpace().ordering(),
          go.trialGridFunctionSpace().ordering(),
          parameters.entriesPerRow()
          );
        go.fill_pattern(pattern);
        _container->setLayout(buildLayout(pattern));
        (*_container) = e;
      }


      //! Creates an FlatELLMatrixContainer without allocating an underlying ISTL matrix.
      explicit BELLMatrixContainer (Dune::PDELab::tags::unattached_container = Dune::PDELab::tags::unattached_container())
      {}

      //! Creates an FlatELLMatrixContainer with an empty underlying ISTL matrix.
      explicit BELLMatrixContainer (Dune::PDELab::tags::attached_container)
        : _container(make_shared<Container>())
      {}

      BELLMatrixContainer(const BELLMatrixContainer& rhs)
        : _container(make_shared<Container>(*(rhs._container)))
      {}

      BELLMatrixContainer& operator=(const BELLMatrixContainer& rhs)
      {
        if (this == &rhs)
          return *this;
        if (attached())
          {
            (*_container) = (*(rhs._container));
          }
        else
          {
            _container = make_shared<Container>(*(rhs._container));
          }
        return *this;
      }

      void detach()
      {
        _container.reset();
      }

      void attach(shared_ptr<Container> container)
      {
        _container = container;
      }

      bool attached() const
      {
        return bool(_container);
      }

      const shared_ptr<Container>& storage() const
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

      BELLMatrixContainer& operator= (const E& e)
      {
        (*_container) = e;
        return *this;
      }

      BELLMatrixContainer& operator*= (const E& e)
      {
        (*_container) *= e;
        return *this;
      }

      E& operator()(const RowIndex& ri, const ColIndex& ci)
      {
        assert(ri.size() == 2);
        assert(ci.size() == 2);
        return (*_container)(ri[1],ci[1],ri[0],ci[0]);
      }

      const E& operator()(const RowIndex& ri, const ColIndex& ci) const
      {
        assert(ri.size() == 2);
        assert(ci.size() == 2);
        return (*_container)(ri[1],ci[1],ri[0],ci[0]);
      }

      const Container& base() const
      {
        return *_container;
      }

      Container& base()
      {
        return *_container;
      }

      void flush()
      {}

      void finalize()
      {}

      void clear_row(const RowIndex& ri, const E& diagonal_entry)
      {
        assert(ri.size() == 2);
        _container->clearRow(ri[1],ri[0]);
        (*this)(ri,ri) = diagonal_entry;
      }

    private:

      shared_ptr<Container> _container;

    };

    } // namespace istl
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_BACKEND_ISTL_BELLMATRIXBACKEND_HH
