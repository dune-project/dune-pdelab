// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_BACKEND_ISTLMATRIXBACKEND_HH
#define DUNE_PDELAB_BACKEND_ISTLMATRIXBACKEND_HH

#include <dune/common/typetraits.hh>
#include <dune/pdelab/backend/tags.hh>
#include <dune/pdelab/backend/common/uncachedmatrixview.hh>
#include <dune/pdelab/backend/istl/matrixhelpers.hh>
#include <dune/pdelab/backend/istl/descriptors.hh>

namespace Dune {
  namespace PDELab {

    template<typename GFSV, typename GFSU, typename C, typename Stats>
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

      typedef Stats PatternStatistics;

#ifndef DOXYGEN

      // some trickery to avoid exposing average users to the fact that there might
      // be multiple statistics objects
      typedef typename conditional<
        (C::blocklevel > 2),
        std::vector<PatternStatistics>,
        PatternStatistics
        >::type StatisticsReturnType;

#endif // DOXYGEN

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
        : _container(std::make_shared<Container>())
      {
        _stats = go.matrixBackend().buildPattern(go,*this);
      }

      /** \brief Construct matrix container using an externally given matrix as storage
       *
       * \tparam GO GridOperator type used to assemble into the matrix
       *
       * \param go GridOperator object used to assemble into the matrix
       * \param container ISTL matrix type that stores the actual data
       *
       * This ISTLMatrixContainer constructor will reassemble the matrix occupation pattern.
       */
      template<typename GO>
      ISTLMatrixContainer (const GO& go, Container& container)
        : _container(Dune::stackobject_to_shared_ptr(container))
      {
        _stats = go.matrixBackend().buildPattern(go,*this);
      }

      template<typename GO>
      ISTLMatrixContainer (const GO& go, const E& e)
        : _container(std::make_shared<Container>())
      {
        _stats = go.matrixBackend().buildPattern(go,*this);
        (*_container) = e;
      }

      //! Creates an ISTLMatrixContainer without allocating an underlying ISTL matrix.
      explicit ISTLMatrixContainer (tags::unattached_container = tags::unattached_container())
      {}

      //! Creates an ISTLMatrixContainer with an empty underlying ISTL matrix.
      explicit ISTLMatrixContainer (tags::attached_container)
        : _container(std::make_shared<Container>())
      {}

      ISTLMatrixContainer(const ISTLMatrixContainer& rhs)
        : _container(std::make_shared<Container>(*(rhs._container)))
      {}

      ISTLMatrixContainer& operator=(const ISTLMatrixContainer& rhs)
      {
        if (this == &rhs)
          return *this;
        _stats.clear();
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

      //! Returns pattern statistics for all contained BCRSMatrix objects.
      const StatisticsReturnType& patternStatistics() const
      {
        return patternStatistics(integral_constant<bool,(C::blocklevel > 2)>());
      }

#ifndef DOXYGEN

    private:

      const PatternStatistics& patternStatistics(false_type multiple) const
      {
        if (_stats.empty())
          DUNE_THROW(InvalidStateException,"no pattern statistics available");
        return _stats[0];
      }

      const std::vector<PatternStatistics>& patternStatistics(true_type multiple) const
      {
        if (_stats.empty())
          DUNE_THROW(InvalidStateException,"no pattern statistics available");
        return _stats;
      }

    public:

#endif

      void detach()
      {
        _container.reset();
        _stats.clear();
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

      std::shared_ptr<Container> _container;
      std::vector<PatternStatistics> _stats;

    };

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_BACKEND_ISTLMATRIXBACKEND_HH
