// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_BACKEND_ISTLMATRIXBACKEND_HH
#define DUNE_PDELAB_BACKEND_ISTLMATRIXBACKEND_HH

#include <dune/pdelab/backend/tags.hh>
#include <dune/pdelab/backend/common/uncachedmatrixview.hh>
#include <dune/pdelab/backend/istl/matrixhelpers.hh>
#include <dune/pdelab/backend/istl/descriptors.hh>

namespace Dune {
  namespace PDELab {

    template<typename GFSV, typename GFSU, typename C>
    class ISTLMatrixContainer
    {

    public:

      typedef typename C::field_type ElementType;
      typedef ElementType E;
      typedef C Container;
      typedef C BaseT;
      typedef typename C::field_type field_type;
      typedef typename C::block_type block_type;
      typedef typename C::size_type size_type;

      typedef GFSU TrialGridFunctionSpace;
      typedef GFSV TestGridFunctionSpace;

      typedef typename GFSV::Ordering::Traits::ContainerIndex RowIndex;
      typedef typename GFSU::Ordering::Traits::ContainerIndex ColIndex;

      typedef typename istl::build_pattern_type<C,GFSV,GFSU,typename GFSV::Ordering::ContainerAllocationTag>::type Pattern;

#if HAVE_TEMPLATE_ALIASES

      template<typename RowCache, typename ColCache>
      using LocalView = UncachedMatrixView<ISTLMatrixContainer,RowCache,ColCache>;

      template<typename RowCache, typename ColCache>
      using ConstLocalView = ConstUncachedMatrixView<const ISTLMatrixContainer,RowCache,ColCache>;

#else

      template<typename RowCache, typename ColCache>
      struct LocalView
        : public UncachedMatrixView<ISTLMatrixContainer,RowCache,ColCache>
      {

        LocalView()
        {}

        LocalView(ISTLMatrixContainer& mc)
          : UncachedMatrixView<ISTLMatrixContainer,RowCache,ColCache>(mc)
        {}

      };

      template<typename RowCache, typename ColCache>
      struct ConstLocalView
        : public ConstUncachedMatrixView<const ISTLMatrixContainer,RowCache,ColCache>
      {

        ConstLocalView()
        {}

        ConstLocalView(const ISTLMatrixContainer& mc)
          : ConstUncachedMatrixView<const ISTLMatrixContainer,RowCache,ColCache>(mc)
        {}

      };

#endif // HAVE_TEMPLATE_ALIASES


      template<typename GO>
      ISTLMatrixContainer (const GO& go)
        : _container(make_shared<Container>())
      {
        Pattern pattern(go.testGridFunctionSpace().ordering(),go.trialGridFunctionSpace().ordering());
        go.fill_pattern(pattern);
        allocate_matrix(go.testGridFunctionSpace().ordering(),
                        go.trialGridFunctionSpace().ordering(),
                        pattern,
                        *_container);
      }

      template<typename GO>
      ISTLMatrixContainer (const GO& go, const E& e)
        : _container(make_shared<Container>())
      {
        Pattern pattern(go.testGridFunctionSpace().ordering(),go.trialGridFunctionSpace().ordering());
        go.fill_pattern(pattern);
        allocate_matrix(go.testGridFunctionSpace().ordering(),
                        go.trialGridFunctionSpace().ordering(),
                        pattern,
                        *_container);
        _container = e;
      }


      //! Creates an ISTLMatrixContainer without allocating an underlying ISTL matrix.
      explicit ISTLMatrixContainer (tags::unattached_container = tags::unattached_container())
      {}

      //! Creates an ISTLMatrixContainer with an empty underlying ISTL matrix.
      explicit ISTLMatrixContainer (tags::attached_container)
        : _container(make_shared<Container>())
      {}

      ISTLMatrixContainer(const ISTLMatrixContainer& rhs)
        : _container(make_shared<Container>(*(rhs._container)))
      {}

      ISTLMatrixContainer& operator=(const ISTLMatrixContainer& rhs)
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

      ISTLMatrixContainer& operator= (const E& e)
      {
        (*_container) = e;
        return *this;
      }

      ISTLMatrixContainer& operator*= (const E& e)
      {
        (*_container) *= e;
        return *this;
      }

      E& operator()(const RowIndex& ri, const ColIndex& ci)
      {
        return istl::access_matrix_element(istl::container_tag(*_container),*_container,ri,ci,ri.size()-1,ci.size()-1);
      }

      const E& operator()(const RowIndex& ri, const ColIndex& ci) const
      {
        return istl::access_matrix_element(istl::container_tag(*_container),*_container,ri,ci,ri.size()-1,ci.size()-1);
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
        istl::clear_matrix_row(istl::container_tag(*_container),*_container,ri,ri.size()-1);
        (*this)(ri,ri) = diagonal_entry;
      }

    private:

      shared_ptr<Container> _container;

    };

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_BACKEND_ISTLMATRIXBACKEND_HH
