// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_ASSEMBLER_FINITEELEMENTWRAPPER_HH
#define DUNE_PDELAB_ASSEMBLER_FINITEELEMENTWRAPPER_HH

#include <cassert>
#include <type_traits>
#include <optional>

#include <dune/common/concept.hh>

#include <dune/localfunctions/common/interfaceswitch.hh>

#include <dune/typetree/leafnode.hh>

#include <dune/pdelab/assembler/utility.hh>

#ifndef DUNE_PDELAB_DEBUG_RANGE_PROXY
#ifdef NDEBUG
#define DUNE_PDELAB_DEBUG_RANGE_PROXY 0
#else
#define DUNE_PDELAB_DEBUG_RANGE_PROXY 1
#endif
#endif

namespace Dune {
  namespace PDELab {

    namespace Concept {

      struct LocalBasis
      {
        template<class B>
        auto require(B&& b) -> decltype(
          Dune::Concept::requireTrue<(B::Traits::dimDomain > 0)>()
        );
      };

    }

    struct OnlyMovable
    {
      OnlyMovable() = default;
      OnlyMovable(const OnlyMovable&) = delete;
      OnlyMovable(OnlyMovable&&) = default;
      OnlyMovable& operator=(const OnlyMovable&) = delete;
      OnlyMovable& operator=(OnlyMovable&&) = delete;
    };

    template<typename FE, typename Cell>
    class FiniteElementWrapper;


    template<typename Container>
    class ExclusiveRangeHolder
    {

    public:

      class Proxy {

      public:

        using value_type = const typename Container::value_type;
        using const_reference = typename Container::const_reference;
        using reference = const_reference;
        using const_iterator = typename Container::const_iterator;
        using iterator   = const_iterator;
        using size_type      = typename Container::size_type;

        size_type size() const
        {
          assert(_range_holder);
          return _range_holder->_container.size();
        }

        iterator begin() const
        {
          assert(_range_holder);
          return _range_holder->_container.begin();
        }

        iterator end() const
        {
          assert(_range_holder);
          return _range_holder->_container.end();
        }

        const_reference operator[](size_type i) const
        {
          assert(_range_holder);
          return _range_holder->_container[i];
        }

        Proxy() noexcept
          : _range_holder(nullptr)
        {}

        Proxy(ExclusiveRangeHolder& range_holder)
          : _range_holder(&range_holder)
        {
#if DUNE_PDELAB_DEBUG_RANGE_PROXY
          ++_range_holder->_proxy_count;
#endif
        }

#if DUNE_PDELAB_DEBUG_RANGE_PROXY
        ~Proxy() noexcept
        {
          if (_range_holder)
            --_range_holder->_proxy_count;
        }

        Proxy(const Proxy& proxy) noexcept
          : _range_holder(proxy._range_holder)
        {
          if (_range_holder)
            ++_range_holder->_proxy_count;
        }

        friend void swap(Proxy& a, Proxy& b) noexcept
        {
          using std::swap;
          swap(a._range_holder,b._range_holder);
        }

        Proxy(Proxy&& proxy) noexcept
          : Proxy()
        {
          swap(*this,proxy);
        }

        Proxy& operator=(const Proxy& proxy) noexcept
        {
          swap(*this,Proxy(proxy));
          return *this;
        }

        Proxy& operator=(Proxy&& proxy) noexcept
        {
          swap(*this,proxy);
          return *this;
        }
#endif

      private:

        const ExclusiveRangeHolder* _range_holder;

      };

      Proxy proxy()
      {
        return Proxy{*this};
      }

      Container& container()
      {
#if DUNE_PDELAB_DEBUG_RANGE_PROXY
        if (_proxy_count > 0)
          DUNE_THROW(Dune::Exception,"Attempt to access underlying container while there are still " << _proxy_count << " proxies around.");
#endif
        return _container;
      }

    private:

      Container _container;

#if DUNE_PDELAB_DEBUG_RANGE_PROXY
      mutable std::size_t _proxy_count;
#endif

    };

    template<typename T, bool = models<Concept::LocalBasis,T>()>
    struct ReferenceGradientType
    {
      using type = int;
    };

    template<typename T>
    struct ReferenceGradientType<T,true>
    {
      using type = typename T::Traits::JacobianType;
    };

    template<typename Basis_,typename Context>
    class BasisWrapper
      : public OnlyMovable
    {

      using Switch = BasisInterfaceSwitch<Basis_>;

      template<typename, typename>
      friend class FiniteElementWrapper;

    public:

      using Basis = Basis_;
      using Native = Basis;
      using DomainField = typename Switch::DomainField;
      static const int dimDomainLocal = Switch::dimDomainLocal;
      using DomainLocal = typename Switch::DomainLocal;
      using RangeField = typename Switch::RangeField;
      static const int dimRange = Switch::dimRange;
      using Range = typename Switch::Range;
      using size_type = std::size_t;
      using Jacobian = typename Native::Traits::JacobianType;
      using ReferenceGradient = typename ReferenceGradientType<Basis_>::type;
      using Gradient = FieldMatrix<RangeField,dimRange,Context::Geometry::coorddimension>;

      const Basis& native() const
      {
        return *_basis;
      }

      size_type size() const
      {
        return _basis->size();
      }

      size_type order() const
      {
        return _basis->order();
      }

      template<typename In, typename Out>
      void evaluateFunction(const In& in, Out& out) const
      {
        _basis->evaluateFunction(in,out);
      }

      template<typename QP>
      typename ExclusiveRangeHolder<std::vector<Range>>::Proxy operator()(const QP& qp)
      {
        _values.container().resize(size());
        _basis->evaluateFunction(Context::Flavor::quadratureCoordinate(qp),_values.container());
        return _values.proxy();
      }

      template<typename QP>
      std::enable_if_t<
        models<Concept::LocalBasis,Native>() and Std::to_true_type<QP>(),
        typename ExclusiveRangeHolder<std::vector<ReferenceGradient>>::Proxy
        >
      referenceGradients(const QP& qp)
      {
        _reference_gradients.container().resize(size());
        _basis->evaluateJacobian(Context::Flavor::quadratureCoordinate(qp),_reference_gradients.container());
        return _reference_gradients.proxy();
      }

      template<typename QP>
      std::enable_if_t<
        models<Concept::LocalBasis,Native>() and dimRange == 1 and Std::to_true_type<QP>(),
        typename ExclusiveRangeHolder<std::vector<Gradient>>::Proxy
        >
      gradients(const QP& qp)
      {
        size_type size = this->size();
        auto& gradients = _gradients.container();
        gradients.resize(size);

        auto reference_gradients = referenceGradients(qp);

        auto jac = _ctx->geometry().jacobianInverseTransposed(Context::Flavor::quadratureCoordinate(qp));

        for(std::size_t i = 0 ; i < size ; ++i)
          jac.mv(reference_gradients[i][0],gradients[i][0]);

        return _gradients.proxy();
      }

      template<typename QP>
      std::enable_if_t<
        not models<Concept::LocalBasis,Native>() and Std::to_true_type<QP>(),
        typename ExclusiveRangeHolder<std::vector<Gradient>>::Proxy
        >
      gradients(const QP& qp)
      {
        size_type size = this->size();
        auto& gradients = _gradients.container();
        gradients.resize(size);

        _basis->evaluateJacobian(Context::Flavor::quadratureCoordinate(qp),gradients);

        return _gradients.proxy();
      }

      BasisWrapper() noexcept
        : _basis(nullptr)
        , _ctx(nullptr)
      {}

    private:

      void setContext(Context& ctx)
      {
        _ctx = &ctx;
      }

      void setBasis(const Native& basis)
      {
        _basis = &basis;
      }

      const Basis* _basis;
      Context* _ctx;
      ExclusiveRangeHolder<std::vector<Range>> _values;
      ExclusiveRangeHolder<std::vector<Gradient>> _gradients;
      ExclusiveRangeHolder<std::vector<ReferenceGradient>> _reference_gradients;

    };


    struct set_finite_elements
      : public TypeTree::TreePairVisitor
      , public TypeTree::DynamicTraversal
    {

      template<typename LFS, typename FE, typename TreePath>
      void leaf(const LFS& lfs, FE& fe, TreePath treePath) const
      {
        fe.setFiniteElement(lfs.finiteElement());
      }

    };

    template<typename Context>
    struct set_context
      : public TypeTree::TreeVisitor
      , public TypeTree::DynamicTraversal
    {

      template<typename FE, typename TreePath>
      void leaf(FE& fe, TreePath treePath) const
      {
        fe.setContext(_ctx);
      }

      set_context(Context& ctx)
        : _ctx(ctx)
      {}

      Context& _ctx;

    };

    template<typename FE, typename Context>
    class FiniteElementWrapper
      : public TypeTree::LeafNode
      , public OnlyMovable
    {

      using Switch = FiniteElementInterfaceSwitch<FE>;
      friend struct set_finite_elements;

      friend struct set_context<Context>;

    public:

      using Basis         = BasisWrapper<typename Switch::Basis,Context>;
      using Interpolation = typename Switch::Interpolation;
      using Coefficients  = typename Switch::Coefficients;

      using FiniteElement = FE;
      using Native = FiniteElement;

      const FiniteElement& native() const
      {
        return *_fe;
      }

      Basis& basis()
      {
        return _basis;
      }

      const Basis& basis() const
      {
        return _basis;
      }

      const Interpolation& interpolation() const
      {
        return Switch::interpolation(*_fe);
      }

      const Coefficients& coefficients() const
      {
        return Switch::coefficients(*_fe);
      }

    private:

      void setContext(Context& ctx)
      {
        _basis.setContext(ctx);
      }

      void setFiniteElement(const Native& fe)
      {
        _fe = &fe;
        _basis.setBasis(Switch::basis(fe));
      }

      const FiniteElement* _fe;
      Basis _basis;

    };

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ASSEMBLER_FINITEELEMENTWRAPPER_HH
