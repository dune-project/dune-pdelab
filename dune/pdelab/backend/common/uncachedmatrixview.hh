// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_BACKEND_COMMON_UNCACHEDMATRIXVIEW_HH
#define DUNE_PDELAB_BACKEND_COMMON_UNCACHEDMATRIXVIEW_HH

#include <type_traits>

namespace Dune {
  namespace PDELab {


    template<typename M_, typename RowCache, typename ColCache>
    class ConstUncachedMatrixView
    {

    public:

      typedef typename std::remove_const<M_>::type Container;

      static_assert(
        (std::is_same<
           typename RowCache::LocalFunctionSpace::Traits::GridFunctionSpace,
           typename Container::TestGridFunctionSpace
           >::value),
        "The RowCache passed to LocalView must belong to the underlying GFSV"
        );

      static_assert(
        (std::is_same<
           typename ColCache::LocalFunctionSpace::Traits::GridFunctionSpace,
           typename Container::TrialGridFunctionSpace
           >::value),
        "The ColCache passed to LocalView must belong to the underlying GFSU"
        );

    public:

      typedef typename Container::field_type E;
      typedef typename Container::size_type size_type;

      typedef E ElementType;

      typedef RowCache RowIndexCache;
      typedef ColCache ColIndexCache;

      typedef typename RowCache::LocalFunctionSpace LFSV;
      typedef typename ColCache::LocalFunctionSpace LFSU;

      typedef typename LFSV::Traits::DOFIndex RowDOFIndex;
      typedef typename LFSV::Traits::GridFunctionSpace::Ordering::Traits::ContainerIndex RowContainerIndex;

      typedef typename LFSU::Traits::DOFIndex ColDOFIndex;
      typedef typename LFSU::Traits::GridFunctionSpace::Ordering::Traits::ContainerIndex ColContainerIndex;

      ConstUncachedMatrixView()
        : _container(nullptr)
        , _row_cache(nullptr)
        , _col_cache(nullptr)
      {}

      ConstUncachedMatrixView(M_& container)
        : _container(&container)
        , _row_cache(nullptr)
        , _col_cache(nullptr)
      {}

      const RowIndexCache& rowIndexCache() const
      {
        assert(_row_cache);
        return *_row_cache;
      }

      const ColIndexCache& colIndexCache() const
      {
        assert(_col_cache);
        return *_col_cache;
      }

      void attach(M_& container)
      {
        _container = &container;
      }

      void detach()
      {
        _container = nullptr;
      }

      void bind(const RowCache& row_cache, const ColCache& col_cache)
      {
        _row_cache = &row_cache;
        _col_cache = &col_cache;
      }

      void unbind()
      {}

      size_type N() const
      {
        return rowIndexCache().size();
      }

      size_type M() const
      {
        return colIndexCache().size();
      }

      template<typename LC>
      void read(LC& local_container) const
      {
        for (size_type i = 0; i < N(); ++i)
          for (size_type j = 0; j < M(); ++j)
            local_container.getEntry(i,j) = container()(rowIndexCache().containerIndex(i),colIndexCache().containerIndex(j));
      }



      const ElementType& operator()(size_type i, size_type j) const
      {
        return container()(rowIndexCache().containerIndex(i),colIndexCache().containerIndex(j));
      }

      // disable this function if DOFIndex and ContainerIndex have the same type - required for interoperability
      // with function spaces based on dune-functions bases
      template<typename RDI, typename CDI>
      std::enable_if_t<
        (std::is_same<RDI,RowDOFIndex>{} and std::is_same<CDI,ColDOFIndex>{} and not
         (std::is_same<RDI,RowContainerIndex>{} and std::is_same<CDI,ColContainerIndex>{})),
        const ElementType&
        >
      operator()(const RDI& i, const CDI& j) const
      {
        return container()(rowIndexCache().containerIndex(i),colIndexCache().containerIndex(j));
      }

      const ElementType& operator()(const RowContainerIndex& i, const ColContainerIndex& j) const
      {
        return container()(i,j);
      }

      const ElementType& operator()(const RowContainerIndex& i, size_type j) const
      {
        return container()(i,colIndexCache().containerIndex(j));
      }

      const ElementType& operator()(size_type i, const ColContainerIndex& j) const
      {
        return container()(rowIndexCache().containerIndex(i),j);
      }

      const Container& container() const
      {
        return *_container;
      }

    protected:

      M_* _container;
      const RowCache* _row_cache;
      const ColCache* _col_cache;

    };


    template<typename M_, typename RowCache, typename ColCache>
    class UncachedMatrixView
      : public ConstUncachedMatrixView<M_,RowCache,ColCache>
    {

      typedef ConstUncachedMatrixView<M_,RowCache,ColCache> BaseT;

    public:

      typedef M_ Container;
      typedef typename Container::ElementType ElementType;
      typedef typename Container::size_type size_type;

      typedef RowCache RowIndexCache;
      typedef ColCache ColIndexCache;

      typedef typename RowCache::LocalFunctionSpace LFSV;
      typedef typename ColCache::LocalFunctionSpace LFSU;

      typedef typename LFSV::Traits::DOFIndex RowDOFIndex;
      typedef typename LFSV::Traits::GridFunctionSpace::Ordering::Traits::ContainerIndex RowContainerIndex;

      typedef typename LFSU::Traits::DOFIndex ColDOFIndex;
      typedef typename LFSU::Traits::GridFunctionSpace::Ordering::Traits::ContainerIndex ColContainerIndex;

      using BaseT::rowIndexCache;
      using BaseT::colIndexCache;
      using BaseT::N;
      using BaseT::M;

      // Explicitly pull in operator() from the base class to work around a problem
      // with clang not finding the const overloads of the operator from the base class.
      using BaseT::operator();

      UncachedMatrixView()
      {}

      UncachedMatrixView(Container& container)
        : BaseT(container)
      {}

      void commit()
      {}

      template<typename LC>
      void write(const LC& local_container)
      {
        for (size_type i = 0; i < N(); ++i)
          for (size_type j = 0; j < M(); ++j)
            container()(rowIndexCache().containerIndex(i),colIndexCache().containerIndex(j)) = local_container.getEntry(i,j);
      }

      template<typename LC>
      void add(const LC& local_container)
      {
        for (size_type i = 0; i < N(); ++i)
          for (size_type j = 0; j < M(); ++j)
            container()(rowIndexCache().containerIndex(i),colIndexCache().containerIndex(j)) += local_container.getEntry(i,j);
      }



      ElementType& operator()(size_type i, size_type j)
      {
        return container()(rowIndexCache().containerIndex(i),colIndexCache().containerIndex(j));
      }

      // disable this function if DOFIndex and ContainerIndex have the same type - required for interoperability
      // with function spaces based on dune-functions bases
      template<typename RDI, typename CDI>
      std::enable_if_t<
        (std::is_same<RDI,RowDOFIndex>{} and std::is_same<CDI,ColDOFIndex>{} and not
         (std::is_same<RDI,RowContainerIndex>{} and std::is_same<CDI,ColContainerIndex>{})),
        const ElementType&
        >
      operator()(const RDI& i, const CDI& j)
      {
        return container()(rowIndexCache().containerIndex(i),colIndexCache().containerIndex(j));
      }

      ElementType& operator()(const RowContainerIndex& i, const ColContainerIndex& j)
      {
        return container()(i,j);
      }

      ElementType& operator()(const RowContainerIndex& i, size_type j)
      {
        return container()(i,colIndexCache().containerIndex(j));
      }

      ElementType& operator()(size_type i, const ColContainerIndex& j)
      {
        return container()(rowIndexCache().containerIndex(i),j);
      }

      void add(size_type i, size_type j, const ElementType& v)
      {
        container()(rowIndexCache().containerIndex(i),colIndexCache().containerIndex(j)) += v;
      }

      // disable this function if DOFIndex and ContainerIndex have the same type - required for interoperability
      // with function spaces based on dune-functions bases
      template<typename RDI, typename CDI>
      std::enable_if_t<
        (std::is_same<RDI,RowDOFIndex>{} and std::is_same<CDI,ColDOFIndex>{} and not
         (std::is_same<RDI,RowContainerIndex>{} and std::is_same<CDI,ColContainerIndex>{}))
        >
      add(const RDI& i, const CDI& j, const ElementType& v)
      {
        container()(rowIndexCache().containerIndex(i),colIndexCache().containerIndex(j)) += v;
      }

      void add(const RowContainerIndex& i, const ColContainerIndex& j, const ElementType& v)
      {
        container()(i,j) += v;
      }

      void add(const RowContainerIndex& i, size_type j, const ElementType& v)
      {
        container()(i,colIndexCache().containerIndex(j)) += v;
      }

      void add(size_type i, const ColContainerIndex& j, const ElementType& v)
      {
        container()(rowIndexCache().containerIndex(i),j) += v;
      }

      Container& container()
      {
        return *(this->_container);
      }

    };


  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_BACKEND_COMMON_UNCACHEDMATRIXVIEW_HH
