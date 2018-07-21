// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_ASSEMBLER_FINITEELEMENTWRAPPER_HH
#define DUNE_PDELAB_ASSEMBLER_FINITEELEMENTWRAPPER_HH

#include <cassert>
#include <type_traits>
#include <optional>

#include <dune/common/concept.hh>
#include <dune/common/hash.hh>

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


  template<typename F, int size>
  inline std::size_t hash_value(const FieldVector<F,size>& v)
  {
    std::size_t seed = 0;
    for (auto x : v)
      hash_combine(seed,x);
    return seed;
  }

}

DUNE_DEFINE_HASH(DUNE_HASH_TEMPLATE_ARGS(typename F, int size),DUNE_HASH_TYPE(Dune::FieldVector<F,size>))

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


    template<typename Implementation, typename Container>
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
          return _range_holder->accessContainer().size();
        }

        iterator begin() const
        {
          assert(_range_holder);
          return _range_holder->accessContainer().begin();
        }

        iterator end() const
        {
          assert(_range_holder);
          return _range_holder->accessContainer().end();
        }

        const_reference operator[](size_type i) const
        {
          assert(_range_holder);
          return _range_holder->accessContainer()[i];
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
          Proxy tmp(proxy);
          swap(*this,tmp);
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
        return accessContainer();
      }

      void assertWriteAccess() const
      {
#if DUNE_PDELAB_DEBUG_RANGE_PROXY
        if (_proxy_count > 0)
          DUNE_THROW(Dune::Exception,"Attempt to access underlying container while there are still " << _proxy_count << " proxies around.");
#endif
      }

    private:

      Container& accessContainer()
      {
        return static_cast<Implementation*>(this)->accessContainer();
      }

      const Container& accessContainer() const
      {
        return static_cast<const Implementation*>(this)->accessContainer();
      }

#if DUNE_PDELAB_DEBUG_RANGE_PROXY
      mutable std::size_t _proxy_count;
#endif

    };

    struct EvaluationProvider
    {

      template<typename BasisWrapper>
      void bind(BasisWrapper& basis_wrapper)
      {}

      template<typename QR>
      void beginQuadrature(QR& qr)
      {}

      template<typename QR>
      void endQuadrature(QR& qr)
      {}

    };

    template<typename Context, typename Evaluation>
    struct UncachedEvaluationStore
      : public ExclusiveRangeHolder<UncachedEvaluationStore<Context,Evaluation>,std::vector<typename Evaluation::Range>>
      , public EvaluationProvider
    {

      using Container = std::vector<typename Evaluation::Range>;
      using Base = ExclusiveRangeHolder<UncachedEvaluationStore,Container>;
      using Base::container;

      template<typename Point>
      void update(const Point& x)
      {
        container().resize(_evaluate.size());
        _evaluate(x,container());
      }

      UncachedEvaluationStore() = default;

      UncachedEvaluationStore(const Evaluation& evaluate)
        : _evaluate(evaluate)
      {}

      template<typename BasisWrapper>
      void setBasisWrapper(BasisWrapper& basis_wrapper)
      {
        _evaluate.setBasisWrapper(basis_wrapper);
      }

      Container& accessContainer()
      {
        return _container;
      }

      const Container& accessContainer() const
      {
        return _container;
      }

      Evaluation _evaluate;
      Container _container;

    };


    template<typename Context, typename Evaluation>
    struct CoordinateCachedEvaluationStore
      : public ExclusiveRangeHolder<CoordinateCachedEvaluationStore<Context,Evaluation>,std::vector<typename Evaluation::Range>>
      , public EvaluationProvider
    {

      using Container = std::vector<typename Evaluation::Range>;
      using Base = ExclusiveRangeHolder<CoordinateCachedEvaluationStore,Container>;
      using Base::container;
      using Base::assertWriteAccess;

      using RuleCache = std::unordered_map<typename Evaluation::Domain,Container>;
      using Key = std::size_t;
      using Cache = std::unordered_map<Key,RuleCache>;

      Container& accessContainer()
      {
        assert(_current);
        return *_current;
      }

      const Container& accessContainer() const
      {
        assert(_current);
        return *_current;
      }

      RuleCache& ruleCache()
      {
        assert(_rule_cache);
        return *_rule_cache;
      }

      template<typename Point>
      void update(const Point& x)
      {
        assertWriteAccess();
        if (auto it = ruleCache().find(Context::Flavor::quadratureCoordinate(x)) ; it != ruleCache().end())
          _current = &it->second;
        else
          {
            _current = &ruleCache()[Context::Flavor::quadratureCoordinate(x)];
            container().resize(_evaluate.size());
            _evaluate(x,container());
          }
      }

      CoordinateCachedEvaluationStore()
        : _current(nullptr)
        , _rule_cache(nullptr)
      {}

      CoordinateCachedEvaluationStore(const Evaluation& evaluate)
        : _evaluate(evaluate)
        , _current(nullptr)
        , _rule_cache(nullptr)
      {}

      template<typename BasisWrapper>
      void setBasisWrapper(BasisWrapper& basis_wrapper)
      {
        _evaluate.setBasisWrapper(basis_wrapper);
      }

      template<typename BasisWrapper>
      void bind(BasisWrapper& basis_wrapper)
      {
        _rule_cache = &_cache[Dune::GlobalGeometryTypeIndex::index(basis_wrapper.context().embedding().global().type())];
        _current = nullptr;
      }

      Evaluation _evaluate;

      Container* _current;
      RuleCache* _rule_cache;
      Cache _cache;

    };


    template<typename Context, typename Evaluation>
    struct IndexCachedEvaluationStore
      : public ExclusiveRangeHolder<IndexCachedEvaluationStore<Context,Evaluation>,std::vector<typename Evaluation::Range>>
      , public EvaluationProvider
    {

      using Container = std::vector<typename Evaluation::Range>;
      using Base = ExclusiveRangeHolder<IndexCachedEvaluationStore,Container>;
      using Base::container;
      using Base::assertWriteAccess;

      using RuleCache = std::vector<Container>;

      struct Key
      {
        std::size_t type;
        std::size_t index;
        std::size_t order;

        inline friend std::size_t hash_value(const Key& key)
        {
          std::size_t seed = 0;
          hash_combine(seed,key.type);
          hash_combine(seed,key.index);
          hash_combine(seed,key.order);
          return seed;
        }

        inline friend bool operator==(const Key& a, const Key& b)
        {
          return a.type == b.type and a.index == b.index and a.order == b.order;
        }

      };

      struct Hasher
      {
        std::size_t operator()(const Key& key) const
        {
          return hash_value(key);
        }
      };

      using Cache = std::unordered_map<Key,RuleCache,Hasher>;

      Container& accessContainer()
      {
        assert(_current);
        return *_current;
      }

      const Container& accessContainer() const
      {
        assert(_current);
        return *_current;
      }

      RuleCache& ruleCache()
      {
        assert(_rule_cache);
        return *_rule_cache;
      }

      template<typename Point>
      void update(const Point& x)
      {
        assertWriteAccess();
        _current = &ruleCache()[x.index()];
      }

      IndexCachedEvaluationStore()
        : _current(nullptr)
        , _rule_cache(nullptr)
      {}

      IndexCachedEvaluationStore(const Evaluation& evaluate)
        : _evaluate(evaluate)
        , _current(nullptr)
        , _rule_cache(nullptr)
      {}

      template<typename BasisWrapper>
      void setBasisWrapper(BasisWrapper& basis_wrapper)
      {
        _evaluate.setBasisWrapper(basis_wrapper);
      }

      template<typename QR>
      void beginQuadrature(QR& qr)
      {
        auto key = Key{
          Dune::GlobalGeometryTypeIndex::index(qr.type()),
          Context::Flavor::embeddingDescriptor(qr.embedding()),
          std::size_t(qr.order())
        };

        if (auto it = _cache.find(key) ; it != _cache.end())
          {
            _rule_cache = &it->second;
          }
        else
          {
            _rule_cache = &_cache[key];
            ruleCache().resize(qr.size());
            auto size = _evaluate.size();
            std::transform(qr.util_begin(),qr.end(),ruleCache().begin(),[&](const auto& x)
              {
                Container temp(size);
                _evaluate(x,temp);
                return std::move(temp);
              });
          }
        _current = nullptr;
      }

      template<typename QR>
      void endQuadrature(QR& qr)
      {
        assertWriteAccess();
        _current = nullptr;
        _rule_cache = nullptr;
      }

      Evaluation _evaluate;

      Container* _current;
      RuleCache* _rule_cache;
      Cache _cache;

    };


    template<typename BasisWrapper, typename Basis_>
    struct ValueEvaluator
    {

      using Basis  = Basis_;
      using Switch = BasisInterfaceSwitch<Basis>;

      using Domain = typename Switch::DomainLocal;
      using Range  = typename Switch::Range;

      template<typename X, typename Container>
      void operator()(const X& x, Container& y) const
      {
        _basis_wrapper->native().evaluateFunction(BasisWrapper::Context::Flavor::quadratureCoordinate(x),y);
      }

      std::size_t size() const
      {
        return _basis_wrapper->size();
      }

      ValueEvaluator() noexcept
        : _basis_wrapper(nullptr)
      {}

      void setBasisWrapper(BasisWrapper& basis_wrapper)
      {
        _basis_wrapper = &basis_wrapper;
      }

    private:

      BasisWrapper* _basis_wrapper;

    };


    template<typename BasisWrapper, typename Basis_>
    struct ReferenceGradientEvaluator
    {

      using Basis  = Basis_;
      using Switch = BasisInterfaceSwitch<Basis>;

      using Domain = typename Switch::DomainLocal;
      using Range  = typename ReferenceGradientType<Basis>::type;

      template<typename X, typename Container>
      void operator()(const X& x, Container& y) const
      {
        static_assert(Std::to_true_v<X> and not std::is_same_v<Range,int>,"basis does not support reference gradients");
        _basis_wrapper->native().evaluateJacobian(BasisWrapper::Context::Flavor::quadratureCoordinate(x),y);
      }

      std::size_t size() const
      {
        return _basis_wrapper->size();
      }

      ReferenceGradientEvaluator() noexcept
        : _basis_wrapper(nullptr)
      {}

      void setBasisWrapper(BasisWrapper& basis_wrapper)
      {
        _basis_wrapper = &basis_wrapper;
      }

    private:

      BasisWrapper* _basis_wrapper;

    };


    template<typename BasisWrapper, typename Basis_, bool haveReferenceGradient = models<Concept::LocalBasis,Basis_>()>
    struct GradientEvaluator
    {

      using Basis  = Basis_;
      using Switch = BasisInterfaceSwitch<Basis>;

      using Domain = typename Switch::DomainLocal;
      using Range  = std::conditional_t<
        Switch::dimRange == 1,
        FieldVector<typename Switch::RangeField,BasisWrapper::Context::Geometry::coorddimension>,
        FieldMatrix<typename Switch::RangeField,Switch::dimRange,BasisWrapper::Context::Geometry::coorddimension>
        >;

      template<typename X, typename Container>
      void operator()(const X& x, Container& y) const
      {
        std::size_t size = _basis_wrapper->size();
        auto reference_gradients = _basis_wrapper->referenceGradients(x);

        auto jac = BasisWrapper::Context::Flavor::cellJacobianInverseTransposed(x);

        for(std::size_t i = 0 ; i < size ; ++i)
          jac.mv(reference_gradients[i][0],y[i]);
      }

      std::size_t size() const
      {
        return _basis_wrapper->size();
      }

      GradientEvaluator() noexcept
        : _basis_wrapper(nullptr)
      {}

      void setBasisWrapper(BasisWrapper& basis_wrapper)
      {
        _basis_wrapper = &basis_wrapper;
      }

    private:

      BasisWrapper* _basis_wrapper;

    };



    template<typename BasisWrapper, typename Basis_>
    struct GradientEvaluator<BasisWrapper,Basis_,false>
    {

      using Basis  = Basis_;
      using Switch = BasisInterfaceSwitch<Basis>;

      using Domain = typename Switch::DomainLocal;
      using Range  = std::conditional_t<
        Switch::dimRange == 1,
        FieldVector<typename Switch::RangeField,BasisWrapper::Context::Geometry::coorddimension>,
        FieldMatrix<typename Switch::RangeField,Switch::dimRange,BasisWrapper::Context::Geometry::coorddimension>
        >;

      template<typename X, typename Container>
      void operator()(const X& x, Container& y) const
      {
        _basis_wrapper->native().evaluateJacobian(BasisWrapper::Context::Flavor::quadratureCoordinate(x),y);
      }

      std::size_t size() const
      {
        return _basis_wrapper->size();
      }

      GradientEvaluator() noexcept
        : _basis_wrapper(nullptr)
      {}

      void setBasisWrapper(BasisWrapper& basis_wrapper)
      {
        _basis_wrapper = &basis_wrapper;
      }

      void setBasis(const Basis& basis)
      {}

    private:

      BasisWrapper* _basis_wrapper;

    };


    template<typename Basis_,typename Context_>
    class BasisWrapper
      : public OnlyMovable
    {

      using Switch = BasisInterfaceSwitch<Basis_>;

      template<typename, typename>
      friend class FiniteElementWrapper;

      using ValueEvaluator             = Dune::PDELab::ValueEvaluator<BasisWrapper,Basis_>;
      using ReferenceGradientEvaluator = Dune::PDELab::ReferenceGradientEvaluator<BasisWrapper,Basis_>;
      using GradientEvaluator          = Dune::PDELab::GradientEvaluator<BasisWrapper,Basis_>;

      //using ValueProvider             = UncachedEvaluationStore<Context_,ValueEvaluator>;
      //using ReferenceGradientProvider = UncachedEvaluationStore<Context_,ReferenceGradientEvaluator>;
      //using ValueProvider             = CoordinateCachedEvaluationStore<Context_,ValueEvaluator>;
      //using ReferenceGradientProvider = CoordinateCachedEvaluationStore<Context_,ReferenceGradientEvaluator>;
      using ValueProvider             = IndexCachedEvaluationStore<Context_,ValueEvaluator>;
      using ReferenceGradientProvider = IndexCachedEvaluationStore<Context_,ReferenceGradientEvaluator>;
      using GradientProvider          = UncachedEvaluationStore<Context_,GradientEvaluator>;

    public:

      using            Basis               = Basis_;
      using            Native              = Basis;
      using            Context             = Context_;
      using            DomainField         = typename Switch::DomainField;
      static const int dimDomainLocal      = Switch::dimDomainLocal;
      using            DomainLocal         = typename Switch::DomainLocal;
      using            RangeField          = typename Switch::RangeField;
      static const int dimRange            = Switch::dimRange;
      using            Range               = typename Switch::Range;
      using            size_type           = std::size_t;
      using            Jacobian            = typename Native::Traits::JacobianType;
      using            RefeferenceGradient = typename ReferenceGradientEvaluator::Range;
      using            Gradient            = typename GradientEvaluator::Range;
      using            Values              = typename ValueProvider::Proxy;
      using            ReferenceGradients  = typename ReferenceGradientProvider::Proxy;
      using            Gradients           = typename GradientProvider::Proxy;

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
        _values.evaluator().evaluateFunction(in,out);
      }


      template<typename QP>
      Values operator()(const QP& qp)
      {
        return values(qp);
      }

      template<typename QP>
      Values values(const QP& qp)
      {
        if (qp.index() != _values_qp_index)
          {
            _values.update(qp);
            _values_qp_index = qp.index();
          }
        return _values.proxy();
      }

      template<typename QP>
      std::enable_if_t<
        models<Concept::LocalBasis,Native>() and Std::to_true_type<QP>(),
        ReferenceGradients
        >
      referenceGradients(const QP& qp)
      {
        if (qp.index() != _reference_gradients_qp_index)
          {
            _reference_gradients.update(qp);
            _reference_gradients_qp_index = qp.index();
          }
        return _reference_gradients.proxy();
      }

      template<typename QP>
      Gradients gradients(const QP& qp)
      {
        if (qp.index() != _gradients_qp_index)
          {
            _gradients.update(qp);
            _gradients_qp_index = qp.index();
          }
        return _gradients.proxy();
      }

      BasisWrapper() noexcept
        : _basis(nullptr)
        , _ctx(nullptr)
        , _values_qp_index(invalid_index)
        , _gradients_qp_index(invalid_index)
        , _reference_gradients_qp_index(invalid_index)
      {}

      Context& context()
      {
        assert(_ctx);
        return *_ctx;
      }

    private:

      void setContext(Context& ctx)
      {
        _ctx = &ctx;
        _values.setBasisWrapper(*this);
        _reference_gradients.setBasisWrapper(*this);
        _gradients.setBasisWrapper(*this);
      }

      void bind(const Native& basis)
      {
        _basis = &basis;
        _values_qp_index = invalid_index;
        _gradients_qp_index = invalid_index;
        _reference_gradients_qp_index = invalid_index;
        _values.bind(*this);
        _reference_gradients.bind(*this);
        _gradients.bind(*this);
      }

      template<typename QR>
      void beginQuadrature(QR& qr)
      {
        _values.beginQuadrature(qr);
        _reference_gradients.beginQuadrature(qr);
        _gradients.beginQuadrature(qr);
      }

      template<typename QR>
      void endQuadrature(QR& qr)
      {
        _values.endQuadrature(qr);
        _reference_gradients.endQuadrature(qr);
        _gradients.endQuadrature(qr);
      }

      static constexpr size_type invalid_index = ~size_type(0);

      const Basis* _basis;
      Context* _ctx;
      ValueProvider _values;
      ReferenceGradientProvider _reference_gradients;
      GradientProvider _gradients;
      size_type _values_qp_index;
      size_type _gradients_qp_index;
      size_type _reference_gradients_qp_index;

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


    template<typename QR>
    struct begin_quadrature
      : public TypeTree::TreeVisitor
      , public TypeTree::DynamicTraversal
    {

      template<typename FE, typename TreePath>
      void leaf(FE& fe, TreePath treePath) const
      {
        fe.beginQuadrature(_qr);
      }

      begin_quadrature(QR& qr)
        : _qr(qr)
      {}

      QR& _qr;

    };

    template<typename QR>
    struct end_quadrature
      : public TypeTree::TreeVisitor
      , public TypeTree::DynamicTraversal
    {

      template<typename FE, typename TreePath>
      void leaf(FE& fe, TreePath treePath) const
      {
        fe.endQuadrature(_qr);
      }

      end_quadrature(QR& qr)
        : _qr(qr)
      {}

      QR& _qr;

    };


    template<typename FE, typename Context>
    class FiniteElementWrapper
      : public TypeTree::LeafNode
      , public OnlyMovable
    {

      using Switch = FiniteElementInterfaceSwitch<FE>;
      friend struct set_finite_elements;

      friend struct set_context<Context>;

      template<typename>
      friend struct begin_quadrature;

      template<typename>
      friend struct end_quadrature;

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
        _basis.bind(Switch::basis(fe));
      }

      template<typename QR>
      void beginQuadrature(QR& qr)
      {
        _basis.beginQuadrature(qr);
      }

      template<typename QR>
      void endQuadrature(QR& qr)
      {
        _basis.endQuadrature(qr);
      }

      const FiniteElement* _fe;
      Basis _basis;

    };

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ASSEMBLER_FINITEELEMENTWRAPPER_HH
