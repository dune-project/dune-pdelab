#ifndef DUNE_PDELAB_BACKEND_ISTL_PRECONDITIONERS_HH
#define DUNE_PDELAB_BACKEND_ISTL_PRECONDITIONERS_HH

#include <memory>
#include <functional>
#include <unordered_map>

#include <dune/common/parametertree.hh>

#include <dune/istl/preconditioners.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/istl/paamg/pinfo.hh>
#include <dune/pdelab/backend/interface.hh>
#include <dune/pdelab/backend/istl/sequentialcomponents.hh>

namespace Dune::Std {

  template<typename T>
  struct type_identity
  {
    using type = T;
  };

  template<typename T>
  using type_identity_t = typename type_identity<T>::type;

}

namespace Dune::PDELab::ISTL::Experimental {

  constexpr inline auto seqILU =
    [](auto provider, const ParameterTree& params, auto old) {
      using Provider = typename decltype(provider)::element_type;
      using Dune::PDELab::Backend::Native;
      using Dune::PDELab::Backend::native;

      using Real = typename Provider::Real;

      auto relaxation_factor = params.get<Real>("relaxation",1.0);
      bool resort = params.get("ilu.resort",false);
      int level = params.get("ilu.level",0);

      return std::make_shared<
        Dune::SeqILU<
          Native<typename Provider::Matrix>,
          Native<typename Provider::Domain>,
          Native<typename Provider::Range>
          >
        >(native(provider->matrix()),level,relaxation_factor,resort);
    };

  constexpr inline auto seqSSOR =
    [](auto provider, const ParameterTree& params, auto old) {
      using Provider = typename decltype(provider)::element_type;
      using Dune::PDELab::Backend::Native;
      using Dune::PDELab::Backend::native;

      using Real = typename Provider::Real;

      auto relaxation_factor = params.get<Real>("relaxation",1.0);
      int iterations = params.get("iterations",3);

      return std::make_shared<
        Dune::SeqSSOR<
          Native<typename Provider::Matrix>,
          Native<typename Provider::Domain>,
          Native<typename Provider::Range>
          >
        >(native(provider->matrix()),iterations,relaxation_factor);
    };

  constexpr inline auto seqSOR =
    [](auto provider, const ParameterTree& params, auto old) {
      using Provider = typename decltype(provider)::element_type;
      using Dune::PDELab::Backend::Native;
      using Dune::PDELab::Backend::native;

      using Real = typename Provider::Real;

      auto relaxation_factor = params.get<Real>("relaxation",1.0);
      int iterations = params.get("iterations",3);

      return std::make_shared<
        Dune::SeqSOR<
          Native<typename Provider::Matrix>,
          Native<typename Provider::Domain>,
          Native<typename Provider::Range>
          >
        >(native(provider->matrix()),iterations,relaxation_factor);
    };

  constexpr inline auto seqGaussSeidel =
    [](auto provider, const ParameterTree& params, auto old) {
      using Provider = typename decltype(provider)::element_type;
      using Dune::PDELab::Backend::Native;
      using Dune::PDELab::Backend::native;

      using Real = typename Provider::Real;

      auto relaxation_factor = params.get<Real>("relaxation",1.0);
      int iterations = params.get("iterations",3);

      return std::make_shared<
        Dune::SeqGS<
          Native<typename Provider::Matrix>,
          Native<typename Provider::Domain>,
          Native<typename Provider::Range>
          >
        >(native(provider->matrix()),iterations,relaxation_factor);
    };

  constexpr inline auto seqJacobi =
    [](auto provider, const ParameterTree& params, auto old) {
      using Provider = typename decltype(provider)::element_type;
      using Dune::PDELab::Backend::Native;
      using Dune::PDELab::Backend::native;

      using Real = typename Provider::Real;

      auto relaxation_factor = params.get<Real>("relaxation",1.0);
      int iterations = params.get("iterations",3);

      return std::make_shared<
        Dune::SeqJac<
          Native<typename Provider::Matrix>,
          Native<typename Provider::Domain>,
          Native<typename Provider::Range>
          >
        >(native(provider->matrix()),iterations,relaxation_factor);
    };

  constexpr inline auto seqRichardson =
    [](auto provider, const ParameterTree& params, auto old) {
      using Provider = typename decltype(provider)::element_type;
      using Dune::PDELab::Backend::Native;
      using Dune::PDELab::Backend::native;

      using Real = typename Provider::Real;

      auto relaxation_factor = params.get<Real>("relaxation",1.0);

      return std::make_shared<
        Dune::Richardson<
          Native<typename Provider::Domain>,
          Native<typename Provider::Range>
          >
        >(relaxation_factor);
    };

  namespace detail {

    template<typename AMG_>
    class AMGWithOperator
      : public Dune::Preconditioner<typename AMG_::Domain,typename AMG_::Domain>
    {

    public:

      using AMG      = AMG_;
      using Matrix   = typename AMG::Operator::matrix_type;
      using Vector   = typename AMG::Domain;
      using Operator = Dune::MatrixAdapter<Matrix,Vector,Vector>;

      void pre(Vector& domain, Vector& range) override
      {
        _amg.pre(domain,range);
      }

      void apply(Vector& domain, const Vector& range) override
      {
        _amg.apply(domain,range);
      }

      void post (Vector& domain) override
      {
        _amg.post(domain);
      }

      SolverCategory::Category category() const override
      {
        return _amg.category();
      }

      template<typename... T>
      explicit AMGWithOperator(const Matrix& matrix, T&&... args)
        : _operator(matrix)
        , _amg(_operator,std::forward<T>(args)...)
      {}

    private:

      Operator _operator;
      AMG _amg;

    };


    constexpr inline auto buildAMGForSmoother =
      [](auto criterion_type, auto smoother_type, auto provider, auto parameters, const ParameterTree& config, auto old) {
        using Provider     = typename decltype(provider)::element_type;
        using Criterion    = typename decltype(criterion_type)::type;
        using Smoother     = typename decltype(smoother_type)::type;

        using Matrix       = Backend::Native<typename Provider::Matrix>;
        using Vector       = Backend::Native<typename Provider::Domain>;
        using Operator     = Dune::MatrixAdapter<Matrix,Vector,Vector>;
        using SmootherArgs = typename Dune::Amg::SmootherTraits<Smoother>::Arguments;
        using AMG          = Dune::Amg::AMG<Operator,Vector,Smoother>;

        SmootherArgs smoother_args;
        smoother_args.iterations = config.get("smoother.iterations",1);
        smoother_args.relaxationFactor = config.get("smoother.relaxation",1.0);

        Criterion criterion(parameters);

        return std::make_shared<AMGWithOperator<AMG>>(Backend::native(provider->matrix()),criterion,smoother_args);
      };

    constexpr inline auto buildAMGForCriterion =
      [](auto criterion_type, auto provider, auto parameters, const ParameterTree& config, auto old)
      -> std::shared_ptr<
        Dune::Preconditioner<
          Backend::Native<typename decltype(provider)::element_type::Domain>,
          Backend::Native<typename decltype(provider)::element_type::Domain>
          >
        >{
        using Criterion = Dune::Amg::CoarsenCriterion<typename decltype(criterion_type)::type>;
        using Provider = typename decltype(provider)::element_type;

        using Matrix = Backend::Native<typename Provider::Matrix>;
        using Vector = Backend::Native<typename Provider::Domain>;

        auto smoother = config.get("smoother.type","ssor");

        if (smoother == "ssor")
          return buildAMGForSmoother(
            Std::type_identity<Criterion>{},
            Std::type_identity<Dune::SeqSSOR<Matrix,Vector,Vector,1>>{},
            provider,
            parameters,
            config,
            old
            );
        else if (smoother == "sor")
          return buildAMGForSmoother(
            Std::type_identity<Criterion>{},
            Std::type_identity<Dune::SeqSOR<Matrix,Vector,Vector,1>>{},
            provider,
            parameters,
            config,
            old
            );
        else
          DUNE_THROW(Exception,"Unknown AMG smoother type: " << smoother);
      };

  }

  constexpr inline auto seqAMG =
    [](auto provider, const ParameterTree& top_config, auto old) {
      using Provider = typename decltype(provider)::element_type;
      using Dune::PDELab::Backend::Native;
      using Dune::PDELab::Backend::native;

      using Matrix       = Native<typename Provider::Matrix>;
      using Vector       = Native<typename Provider::Domain>;
      //using Operator     = Dune::MatrixAdapter<Matrix,Vector,Vector>;
      //using Smoother     = Dune::SeqSSOR<Matrix,Vector,Vector,1>;
      //using SmootherArgs = typename Dune::Amg::SmootherTraits<Smoother>::Arguments;
      //using AMG          = Dune::Amg::AMG<Operator,Vector,Smoother>;
      using Parameters   = Dune::Amg::Parameters;

      static_assert(std::is_same_v<Vector,Native<typename Provider::Range>>,"range and domain must be the same");

      constexpr int dim = Provider::Domain::GridFunctionSpace::Traits::EntitySet::dimension;

      auto& config = top_config.sub("amg");

      Parameters parameters;

      auto diameter = config.get("amg.diameter",2);

      if (config.hasKey("preset"))
      {
        auto preset = config["preset"];
        if (preset == "isotropic")
          parameters.setDefaultValuesIsotropic(dim,diameter);
        else if (preset == "anisotropic")
          parameters.setDefaultValuesAnisotropic(dim,diameter);
        else
          DUNE_THROW(Exception,"Unknown AMG preset: " << preset);
      }

      if (config.hasKey("max-distance"))
        parameters.setMaxDistance(config.get<std::size_t>("max-distance"));

      if (config.hasKey("skip-isolated"))
        parameters.setMaxDistance(config.get<bool>("skip-isolated"));

      if (config.hasKey("min-aggregate-size"))
        parameters.setMinAggregateSize(config.get<std::size_t>("min-aggregate-size"));

      if (config.hasKey("max-aggregate-size"))
        parameters.setMinAggregateSize(config.get<std::size_t>("max-aggregate-size"));

      if (config.hasKey("max-connectivity"))
        parameters.setMinAggregateSize(config.get<std::size_t>("max-connectivity"));

      if (config.hasKey("alpha"))
        parameters.setAlpha(config.get<double>("alpha"));

      if (config.hasKey("beta"))
        parameters.setBeta(config.get<double>("beta"));

      if (config.hasKey("max-level"))
        parameters.setMaxLevel(config.get<int>("max-level"));

      if (config.hasKey("coarsen-target"))
        parameters.setCoarsenTarget(config.get<int>("coarsen-target"));

      if (config.hasKey("min-coarsen-rate"))
        parameters.setMinCoarsenRate(config.get<double>("min-coarsen-rate"));

      // TODO: accumulation mode!

     if (config.hasKey("prolongation-damping-factor"))
        parameters.setProlongationDampingFactor(config.get<double>("prolongation-damping-factor"));

     if (config.hasKey("debug-level"))
        parameters.setDebugLevel(config.get<int>("debug-level"));

     if (config.hasKey("pre-smooth-steps"))
       parameters.setNoPreSmoothSteps(config.get<std::size_t>("pre-smooth-steps"));

     if (config.hasKey("post-smooth-steps"))
       parameters.setNoPostSmoothSteps(config.get<std::size_t>("post-smooth-steps"));

     if (config.hasKey("gamma"))
       parameters.setGamma(config.get<std::size_t>("gamma"));

     if (config.hasKey("additive"))
       parameters.setAdditive(config.get<bool>("additive"));

     auto criterion = config.get("criterion","symmetric");

     if (criterion == "symmetric")
       return detail::buildAMGForCriterion(
         Std::type_identity<Dune::Amg::SymmetricCriterion<Matrix,Dune::Amg::FirstDiagonal>>{},
         provider,
         parameters,
         config,
         old
         );
     else if (criterion == "unsymmetric")
       return detail::buildAMGForCriterion(
         Std::type_identity<Dune::Amg::UnSymmetricCriterion<Matrix,Dune::Amg::FirstDiagonal>>{},
         provider,
         parameters,
         config,
         old
         );
     else
       DUNE_THROW(Exception,"Unknown criterion: " << criterion);
    };


  template<typename MatrixProvider_>
  class PreconditionerRepository
  {

  public:

    using MatrixProvider  = MatrixProvider_;
    using Matrix          = typename MatrixProvider::Matrix;
    using Domain          = typename MatrixProvider::Domain;
    using Range           = typename MatrixProvider::Range;
    using Preconditioner  = Dune::Preconditioner<Backend::Native<Domain>,Backend::Native<Range>>;

    using Factory         = std::function<
      std::shared_ptr<Preconditioner>(
        std::shared_ptr<MatrixProvider>,
        const ParameterTree&,
        std::shared_ptr<Preconditioner>
        )>;

    using FactoryWrapper  = std::function<
      std::shared_ptr<Preconditioner>(
        std::shared_ptr<MatrixProvider>,
        std::shared_ptr<Preconditioner>
        )>;

    using Repository = std::unordered_map<std::string,Factory>;

    static Repository& repository()
    {
      // create repository with default methods
      static Repository repo {
        { "ilu" , seqILU },
        { "ssor" , seqSSOR },
        { "sor" , seqSOR },
        { "gauss-seidel" , seqGaussSeidel },
        { "jacobi" , seqJacobi },
        { "richardson" , seqRichardson },
        { "amg" , seqAMG }
      };
      return repo;
    }

  public:

    void add(const std::string& name, Factory factory) const
    {
      repository().emplace(name,std::move(factory));
    }

    Factory get(const std::string& name) const
    {
      try {
        return repository().at(name);
      } catch (std::out_of_range&) {
        DUNE_THROW(Exception, "Could not find preconditioner: " << name);
      }
    }

    FactoryWrapper get(const ParameterTree& config) const
    {
      return [&config,factory=get(config["type"])](std::shared_ptr<MatrixProvider> provider, std::shared_ptr<Preconditioner> old)
             {
               return factory(std::move(provider),config.sub(config["type"]),old);
             };
    }

  };


} // namespace Dune::PDELab::ISTL::Experimental

#endif // DUNE_PDELAB_BACKEND_ISTL_PRECONDITIONERS_HH
