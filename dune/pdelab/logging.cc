#include "config.h"

#include <algorithm>
#include <cassert>
#include <functional>
#include <iterator>
#include <memory>
#include <optional>
#include <regex>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <dune/pdelab/common/utility.hh>
#include <dune/pdelab/logging.hh>
#include <dune/pdelab/logging/loggerbackend.hh>
#include <dune/pdelab/logging/filesinks.hh>

namespace Dune::PDELab {

  ////////////////////////////////////////////////////////////////////////////////
  // sink factories
  ////////////////////////////////////////////////////////////////////////////////

  // Factory function for file sinks
  // Handles splitting into multiple files on parallel runs
  static auto fileSinkFactory(
    std::string_view name,
    LogLevel level,
    int widest_logger,
    const ParameterTree& params
    )
  {
    if (not params.hasKey("file"))
      DUNE_THROW(LoggingError,"You must specify an output file name for file sink: " << name);
    auto file_name = params["file"];
    if (file_name.empty())
      DUNE_THROW(LoggingError,"You must specify an output file name for file sink: " << name);

    auto mode = std::ios::trunc;

    if (params.hasKey("mode"))
    {
      auto mode_name = params["mode"];
      if (mode_name == "truncate")
        mode = std::ios::trunc;
      else if (mode_name == "append")
        mode = std::ios::ate;
      else
        DUNE_THROW(LoggingError,"Unknown file open mode " << mode_name << ": " << name);
    }

    auto comm = Logging::comm();
    auto size = comm.size();

    auto pattern_in_name = file_name.find("{}") != std::string::npos;

    if (comm.size() > 1 or pattern_in_name)
    {
      auto size_digits = 1;
      while (size / 10 > 0)
        ++size_digits;

      if (file_name.find("{}") != std::string::npos)
        file_name = fmt::format(file_name,fmt::format("{0:{1}}",comm.rank(),size_digits));
      else
        file_name = fmt::format("{1:0>{2}}-{0}"_fmt,file_name,comm.rank(),size_digits);
    }

    auto sink = std::make_shared<FileSink>(name,level,widest_logger,file_name,mode);

    PatternFormatSink::setParameters(*sink,params);
    return sink;
  }

  ////////////////////////////////////////////////////////////////////////////////
  // internal structs
  ////////////////////////////////////////////////////////////////////////////////

  struct Logging::SinkFactoryRepository
    : public std::unordered_map<std::string,Logging::SinkFactory>
  {
    using Base = std::unordered_map<std::string,Logging::SinkFactory>;
    using Base::Base;
  };

  struct Logging::State
  {

    using BackendRegistry = std::unordered_map<std::string_view,std::unique_ptr<LoggerBackend>>;
    using SinkRegistry = std::unordered_map<std::string_view,std::shared_ptr<Sink>>;

    using NameStorage = std::unordered_set<std::string>;

    // In order to use string_view in the map keys, we need to make sure the underlying string
    // doesn't go away, so we store all of them in a set.
    std::string_view internalize(std::string_view name)
    {
      auto [internalized,inserted] = name_storage.emplace(name);
      return *internalized;
    }

    State(const State&) = delete;
    State& operator=(const State&) = delete;

    // The state constructor makes sure that the initialization happens from a call
    // to state() that got passed a pointer to a comm object, which only happens in
    // init(). This allows us to catch uninitialized usage of the system.
    State(const CollectiveCommunication* comm_)
    {
      if (not comm_)
        DUNE_THROW(LoggingError,"You must call Dune::PDELab::Logging::init() before using the logging system");
      comm.emplace(*comm_);
    }

    NameStorage name_storage;
    BackendRegistry backends;
    SinkRegistry sinks;
    LoggerBackend* default_backend = nullptr;
    std::size_t widest_logger = 0;
    std::shared_ptr<ConsoleSink> cout;
    std::shared_ptr<ConsoleSink> cerr;
    bool muted = false;
    LogLevel unmuted_cout = LogLevel::all;
    LogLevel unmuted_cerr = LogLevel::all;
    std::optional<CollectiveCommunication> comm;
    LogMessage::Time startup_time = LogMessage::Clock::now();

  };

  ////////////////////////////////////////////////////////////////////////////////
  // internal helper functionality
  ////////////////////////////////////////////////////////////////////////////////

  // helper function that reaises an exception if passed an invalid name
  static void validateName(std::string_view name)
  {
    if (not Logging::isValidName(name))
      DUNE_THROW(LoggingError,"Invalid name for logging component: " << name);
  }

  // This function does the one-time setup for the sink factory, mostly
  // registering default factories.
  Logging::SinkFactoryRepository Logging::makeSinkFactoryRepository()
  {
    Logging::SinkFactoryRepository repo;
    repo["file"] = fileSinkFactory;
    return repo;
  }

  Logging::SinkFactoryRepository& Logging::sinkFactoryRepository()
  {
    static SinkFactoryRepository repository = makeSinkFactoryRepository();
    return repository;
  }

  Logging::SinkFactory& Logging::sinkFactory(const std::string& name)
  {
    return sinkFactoryRepository().at(name);
  }

  Logging::State& Logging::state(const CollectiveCommunication* comm)
  {
    static State state(comm);
    return state;
  }

  LoggerBackend& Logging::backend(std::string_view name)
  {
    using namespace std::literals;

    // The empty name gets forwarded to "default"
    if (name.empty())
      name = "default"sv;
    try {
      return *state().backends.at(name);
    } catch (std::out_of_range&) {
      DUNE_THROW(LoggingError,"Could not find backend in registry: " << name);
    }
  }

  ////////////////////////////////////////////////////////////////////////////////
  // initialization
  ////////////////////////////////////////////////////////////////////////////////

  void Logging::init(const CollectiveCommunication& comm, const ParameterTree& params)
  {
    // initialize singleton
    auto& s = state(&comm);

    // create default sinks for stdout and stderr
    {
      std::string default_pattern = "{0}";
      auto        default_level   = LogLevel::all;

      auto level = default_level;
      if (params.hasKey("sink.stdout.level"))
        level = parseLogLevel(params["sink.stdout.level"]);

      auto pattern = params.get("sink.stdout.pattern",default_pattern);

      s.cout = std::make_shared<ConsoleSink>("stdout",stdout,pattern,level,0);
      registerSink(s.cout);

      level = default_level;
      if (params.hasKey("sink.stderr.level"))
        level = parseLogLevel(params["sink.stderr.level"]);

      pattern = params.get("sink.stderr.pattern",default_pattern);

      s.cerr = std::make_shared<ConsoleSink>("stderr",stderr,pattern,level,0);
      registerSink(s.cerr);
    }

    // create custom sinks
    if (params.hasSub("sink"))
    {
      auto& sinks = params.sub("sink");
      for (auto& name : sinks.getSubKeys())
      {
        if (name == "stdout" or name == "stderr")
          continue;
        validateName(name);
        auto& config = sinks.sub(name);
        makeSink(name,config);
      }
    }

    // create default logger backend
    {
      auto level = LogLevel::notice;
      if (params.hasKey("default.level"))
        level = parseLogLevel(params["default.level"]);

      auto enabled = params.get("default.enabled",true);
      auto indent  = params.get("default.indent",0);

      auto internalized_name = s.internalize("default");
      s.default_backend = s.backends.emplace(
        internalized_name,
        new LoggerBackend(
          "default",
          s.startup_time,
          enabled,
          level,
          indent
          )
        ).first->second.get();

      if (params.hasKey("default.sinks"))
      {
        auto sinks = parseConfigList(params["default.sinks"]);
        for (auto name : sinks)
          s.default_backend->_sinks.push_back(sink({begin(name),end(name)}));
      }
      else
      {
        // log to stdout by default
        s.default_backend->_sinks.push_back(cout());
      }

      if (params.hasKey("default.extra_sinks")){
        auto sinks = parseConfigList(params["default.extra_sinks"]);
        for (auto name : sinks)
          s.default_backend->_sinks.push_back(sink({begin(name),end(name)}));
      }
    }

    // create logger backends
    if (params.hasSub("backend"))
    {
      auto& backends = params.sub("backend");

      // create backends with default configuration (inherited from default logger)
      for (auto& name : backends.getValueKeys())
      {

        validateName(name);

        if (backends[name] != "default")
          DUNE_THROW(LoggingError,"When declaring a logger with default configuration, its assigned "
                                  "value must be \"default\": " << name);

        if (name == "default")
          DUNE_THROW(LoggingError,"You cannot declare the default logger as \"default\"");

        registerBackend(name, s.default_backend->_default_level, true);
      }

      // create backends with custom configuration
      for (auto& name : backends.getSubKeys())
      {
        validateName(name);

        if (name == "default")
          DUNE_THROW(
            LoggingError,
            "You cannot create a custom logger backend with the reserved name \"default\". Change "
            "the default backend configuration using the key group \"default\" without prepending "
            "\"backend.\"");

        auto& config = backends.sub(name);

        LogLevel level = s.default_backend->_default_level;
        if (config.hasKey("level"))
          level = parseLogLevel(config["level"]);

        bool enabled = config.get("enabled",s.default_backend->_enabled);
        int indent   = config.get("indent",s.default_backend->_default_indent);

        auto& backend = *registerBackend(name,level,not config.hasKey("sinks"))._backend;
        backend._default_indent = indent;
        backend._enabled = enabled;

        if (config.hasKey("sinks"))
        {
          auto sinks = parseConfigList(config["sinks"]);
          for (auto name : sinks)
            backend._sinks.push_back(sink({begin(name),end(name)}));
        }

        if (config.hasKey("extra_sinks")){
          auto sinks = parseConfigList(params["extra_sinks"]);
          for (auto name : sinks)
            s.default_backend->_sinks.push_back(sink({begin(name),end(name)}));
        }
      }
    }

    if (params.get("muted",s.comm->rank() > 0))
    {
      mute();
      if (s.comm->size() > 0)
        logger().info("Muted console log sinks on MPI ranks > 0"_fmt);
    }

    logger().notice("Logging system initialized"_fmt);
  }

  // Make not to use state() in here, this function must work before init().
  void Logging::registerSinkFactory(const std::string& name, SinkFactory sink_factory)
  {
    validateName(name);
    auto factories = sinkFactoryRepository();
    if (factories.count(name) > 0)
      DUNE_THROW(LoggingError, "Cannot register sink factory, name already used: " << name);
    factories[name] = sink_factory;
  }

  ////////////////////////////////////////////////////////////////////////////////
  // Sinks
  ////////////////////////////////////////////////////////////////////////////////

  std::shared_ptr<Sink> Logging::makeSink(const std::string& name, const ParameterTree& params)
  {
    validateName(name);
    auto& sinks = state().sinks;
    if (sinks.count(name) > 0)
      DUNE_THROW(LoggingError, "Cannot register sink, name already used: " << name);

    auto level = LogLevel::all;
    if (params.hasKey("level"))
      level = parseLogLevel(params["level"]);

    auto [sink,_] = sinks.emplace(name,sinkFactory(params["type"])(name,level,state().widest_logger,params));
    return sink->second;;
  }

  void Logging::registerSink(std::shared_ptr<Sink> sink)
  {
    validateName(sink->name());
    auto &sinks = state().sinks;
    if (sinks.count(sink->name()) > 0)
      DUNE_THROW(LoggingError,
                 "Cannot register sink, name already used: " << sink->name());
    sinks[sink->name()] = sink;
    sink->setWidestLogger(state().widest_logger);
  }

  std::shared_ptr<Sink> Logging::sink(const std::string& name)
  {
    try {
      return state().sinks.at(name);
    } catch (std::out_of_range&) {
      DUNE_THROW(LoggingError,"Could not find sink in registry: " << name);
    }
  }

  bool Logging::retireSink(std::string_view sink)
  {
    return state().sinks.erase(sink) > 0;
  }

  std::shared_ptr<ConsoleSink> Logging::cout()
  {
    return state().cout;
  }

  std::shared_ptr<ConsoleSink> Logging::cerr()
  {
    return state().cerr;
  }

  ////////////////////////////////////////////////////////////////////////////////
  // Loggers
  ////////////////////////////////////////////////////////////////////////////////

  Logger Logging::logger()
  {
    auto& backend = *state().default_backend;
    return {backend,backend._default_level,backend._default_indent};
  }

  Logger Logging::logger(std::string_view name)
  {
    using namespace std::literals;
    if (name.empty())
      name = "default"sv;
    auto& backend = *state().backends.at(name);
    return {backend,backend._default_level,backend._default_indent};
  }

  Logger Logging::logger(const ParameterTree& params)
  {
    auto log = logger(params.get<std::string>("log.backend",""));
    if (params.hasKey("log.level"))
      log.setLevel(parseLogLevel(params["log.level"]));
    if (params.hasKey("log.indent"))
      log.setIndent(params.get<int>("log.indent"));
    return log;
  }

  ////////////////////////////////////////////////////////////////////////////////
  // Backends
  ////////////////////////////////////////////////////////////////////////////////

  Logger Logging::registerBackend(
    std::string_view name,
    LogLevel level,
    bool attach_default_sinks
    )
  {
    auto& s = state();

    validateName(name);

    if (name == "default")
      DUNE_THROW(LoggingError, "You cannot create a custom logger backend with the reserved name \"default\".");

    auto internalized_name = s.internalize(name);
    auto& backend = *s.backends.emplace(
      internalized_name,
      new LoggerBackend(
        internalized_name,
        s.startup_time,
        s.default_backend->_enabled,
        level,
        s.default_backend->_default_indent
        )
      ).first->second;

    if (attach_default_sinks)
    {
      // inherit the sink configuration from the default backend
      auto& default_sinks = s.default_backend->_sinks;
      std::copy(begin(default_sinks),end(default_sinks),std::back_inserter(backend._sinks));
    }

    if (name.length() > s.widest_logger)
    {
      s.widest_logger = name.length();
      // inform all sinks about the new longest logger name
      for (auto& [_,sink] : s.sinks)
        sink->setWidestLogger(s.widest_logger);
    }

    return {backend,backend._default_level,backend._default_indent};
  }

  bool Logging::attachSink(std::string_view backend_name, std::string_view sink_name)
  {
    auto& b = backend(backend_name);
    auto  s = sink({begin(sink_name),end(sink_name)});

    if (std::any_of(begin(b._sinks),end(b._sinks),[=](auto& sink) { return sink->name() == sink_name; }))
      return false;

    b._sinks.push_back(s);

    return true;
  }

  bool Logging::detachSink(std::string_view backend_name, std::string_view sink_name)
  {
    auto& b = backend(backend_name);
    auto  s = sink({begin(sink_name),end(sink_name)});

    if (auto it = std::remove_if(
          begin(b._sinks),
          end(b._sinks),
          [=](auto& sink) { return sink->name() == sink_name; }
          ) ;
        it != end(b._sinks)
      )
    {
      b._sinks.erase(it,end(b._sinks));
      return true;
    }

    return false;
  }

  ////////////////////////////////////////////////////////////////////////////////
  // Muting
  ////////////////////////////////////////////////////////////////////////////////

  bool Logging::muted()
  {
    return state().muted;
  }

  void Logging::mute()
  {
    auto& s = state();
    s.unmuted_cout = s.cout->level();
    s.cout->setLevel(LogLevel::off);
    s.unmuted_cerr = s.cerr->level();
    s.cerr->setLevel(LogLevel::off);
    s.muted = true;
  }

  void Logging::unmute()
  {
    auto& s = state();
    s.cout->setLevel(s.unmuted_cout);
    s.cerr->setLevel(s.unmuted_cerr);
    s.muted = false;
  }

  ////////////////////////////////////////////////////////////////////////////////
  // Utilities
  ////////////////////////////////////////////////////////////////////////////////

  bool Logging::isValidName(std::string_view name)
  {
    static const std::regex valid_pattern("[a-z](-?[a-z0-9])*",std::regex::nosubs);
    return std::regex_match(begin(name),end(name),valid_pattern);
  }

  const Logging::CollectiveCommunication& Logging::comm()
  {
    return *state().comm;
  }




} // end namespace Dune::PDELab
