#ifndef DUNE_PDELAB_BACKEND_ISTL_SOLVERS_HH
#define DUNE_PDELAB_BACKEND_ISTL_SOLVERS_HH

#include <memory>
#include <unordered_map>

#include <dune/common/parametertree.hh>

#include <dune/logging.hh>

#include <dune/istl/solvers.hh>
#include <dune/pdelab/backend/interface.hh>
#include <dune/pdelab/backend/istl/seqistlsolverbackend.hh>
#include <dune/pdelab/backend/istl/preconditioners.hh>
#include <dune/pdelab/backend/istl/utility.hh>

namespace Dune::PDELab::ISTL::Experimental {

  constexpr inline auto cg =
    [](auto lin_op, auto prec, const ParameterTree& config, Logging::Logger log) {
      using LinearOperator = typename decltype(lin_op)::element_type;

      double reduction = config.get("reduction",1e-10);
      int max_iterations = config.get("max-iterations",1000);
      int verbose = logLevelToVerbosity(log.level());

      log.info("Creating CG solver, reduction={}, max_iterations={}, verbose={}"_fmt,
               reduction,max_iterations,verbose);

      return std::make_shared<
        Dune::CGSolver<
          typename LinearOperator::Domain
          >
        >(*lin_op,*prec,reduction,max_iterations,verbose);
    };

  constexpr inline auto biCGSTAB =
    [](auto lin_op, auto prec, const ParameterTree& config, Logging::Logger log) {
      using LinearOperator = typename decltype(lin_op)::element_type;

      double reduction = config.get("reduction",1e-10);
      int max_iterations = config.get("max-iterations",1000);
      int verbose = logLevelToVerbosity(log.level());

      log.info("Creating BiCGStab solver, reduction={}, max_iterations={}, verbose={}"_fmt,
               reduction,max_iterations,verbose);

      return std::make_shared<
        Dune::BiCGSTABSolver<
          typename LinearOperator::Domain
          >
        >(*lin_op,*prec,reduction,max_iterations,verbose);
    };

  constexpr inline auto loop =
    [](auto lin_op, auto prec, const ParameterTree& config, Logging::Logger log) {
      using LinearOperator = typename decltype(lin_op)::element_type;

      double reduction = config.get("reduction",1e-10);
      int max_iterations = config.get("max-iterations",1000);
      int verbose = logLevelToVerbosity(log.level());

      log.info("Creating loop solver, reduction={}, max_iterations={}, verbose={}"_fmt,
               reduction,max_iterations,verbose);


      return std::make_shared<
        Dune::LoopSolver<
          typename LinearOperator::Domain
          >
        >(*lin_op,*prec,reduction,max_iterations,verbose);
    };

  constexpr inline auto minRes =
    [](auto lin_op, auto prec, const ParameterTree& config, Logging::Logger log) {
      using LinearOperator = typename decltype(lin_op)::element_type;

      double reduction = config.get("reduction",1e-10);
      int max_iterations = config.get("max-iterations",1000);
      int verbose = logLevelToVerbosity(log.level());

      log.info("Creating MinRes solver, reduction={}, max_iterations={}, verbose={}"_fmt,
               reduction,max_iterations,verbose);

      return std::make_shared<
        Dune::MINRESSolver<
          typename LinearOperator::Domain
          >
        >(*lin_op,*prec,reduction,max_iterations,verbose);
    };

  constexpr inline auto restartedGMRes =
    [](auto lin_op, auto prec, const ParameterTree& config, Logging::Logger log) {
      using LinearOperator = typename decltype(lin_op)::element_type;

      double reduction = config.get("reduction",1e-10);
      int max_iterations = config.get("max-iterations",1000);
      int verbose = logLevelToVerbosity(log.level());
      int restart = config.get("gmres.restart",30);

      log.info("Creating GMRes solver, reduction={}, max_iterations={}, restart={}. verbose={}"_fmt,
               reduction,max_iterations,restart,verbose);

      return std::make_shared<
        Dune::RestartedGMResSolver<
          typename LinearOperator::Domain
          >
        >(*lin_op,*prec,reduction,restart,max_iterations,verbose);
    };


  template<typename Domain_, typename Range_ = Domain_>
  class LinearSolverRepository
  {

  public:

    using Domain          = Domain_;
    using Range           = Range_;
    using Solver          = Dune::IterativeSolver<Domain,Range>;
    using Preconditioner  = Dune::PDELab::ISTL::Preconditioner<Domain,Range>;
    using LinearOperator  = Dune::PDELab::ISTL::LinearOperator<Domain,Range>;
    using Factory         = std::function<
      std::shared_ptr<Solver>(
        std::shared_ptr<LinearOperator>,
        std::shared_ptr<Preconditioner>,
        const ParameterTree&,
        Logging::Logger
        )>;

  private:

    using Repository = std::unordered_map<std::string,Factory>;

    static Repository& repository()
    {
      // create repository with default methods
      static Repository repo {
        { "cg" , cg },
        { "bicgstab" , biCGSTAB },
        { "loop" , loop },
        { "minres"     , minRes },
        { "gmres"           , restartedGMRes }
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
        DUNE_THROW(Exception, "Could not find one step method \"" << name <<  "\"");
      }
    }

  };


  template<typename LO>
  auto makeIterativeLinearSolver(
    std::shared_ptr<LO> lin_op,
    const ParameterTree& config,
    Logging::Logger log = Logging::Logger(),
    std::shared_ptr<ISTL::Preconditioner<typename LO::Domain,typename LO::Range>> prec = nullptr
    )
    -> std::shared_ptr<Backend::LinearSolver<typename LO::Domain,typename LO::Range>>
  {
    if (not log.attached())
      log = Logging::componentLogger(config,"linalg");

    if (not prec)
    {
      if (not config.hasSub("preconditioner"))
      {
        ParameterTree empty_config;
        prec = makeMatrixFreeSequentialPreconditioner(lin_op,seqRichardson,empty_config);
      }
      else
      {
        using MatrixProvider = Backend::MatrixBasedLinearSolverComponent<typename LO::Matrix>;
        auto provider = std::dynamic_pointer_cast<MatrixProvider>(lin_op);
        if (not provider)
          DUNE_THROW(Exception,"Cannot construct preconditioner for matrix-free operator");

        prec = makeMatrixBasedSequentialPreconditioner(
          provider,
          PreconditionerRepository<MatrixProvider>().get(config.sub("preconditioner")),
          config.sub("preconditioner")
          );
      }
    }

    using Repository = LinearSolverRepository<typename LO::Domain,typename LO::Range>;
    using Factory = typename Repository::Factory;

    return std::make_shared<
      ISTL::IterativeLinearSolver<
        typename LO::Domain,
        typename LO::Range,
        Factory
        >
      >(lin_op,prec,Repository().get(config["type"]),config,log);

  }

} // namespace Dune::PDELab::ISTL::Experimental

#endif // DUNE_PDELAB_BACKEND_ISTL_SOLVERS_HH
