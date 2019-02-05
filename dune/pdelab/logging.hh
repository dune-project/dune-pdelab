#ifndef DUNE_PDELAB_LOGGING_HH
#define DUNE_PDELAB_LOGGING_HH

#include <memory>
#include <string_view>

#include <dune/common/parametertree.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/pdelab/logging/logmessage.hh>
#include <dune/pdelab/logging/sink.hh>
#include <dune/pdelab/logging/consolesink.hh>
#include <dune/pdelab/logging/logger.hh>
#include <dune/pdelab/logging/logger.hh>

namespace Dune::PDELab {

/**
 * \addtogroup logging Logging Infrastructure
 * \brief Logging of status and error messages to different targets with
 * ParameterTree-based configuration.
 *
 * # Overview
 *
 * PDELab has a powerful, flexible logging system for reporting the current
 * status of a running program. The system revolves around a small number of
 * components:
 *
 * - The Logger is the most important class for normal users. Messages are
 * logged with the help of a Logger. Loggers are lightweight objects with value
 * semantics, so you can (and should) copy them around. Loggers can also be
 * default-constructed into an invalid state and then assigned a valid logger at
 * a later time.
 *
 * - The logger backend is a mostly invisible component to which a Logger
 * forwards its log messages. Logger backends are long-lived singletons stored
 * in the logging system that carry information like the list of attached sinks
 * and defaults for new Loggers. Many Loggers can share a single logger backend,
 * and every Logger is attached to exactly one backend.
 *
 * - The Sink is a class that is responsible for actually storing a log message
 * somewhere. It gets passed preprocessed log messages in the form of a
 * LogMessage object and has to write them to the terminal, a file, over the
 * network or whatever you come up with. Like backends, Sinks are managed by the
 * logging system, but you can obtain a shared_ptr to a registered sink and
 * interact with it, e.g. to change its configuration. Unlike backends, Sinks
 * can also be removed from the logging system by calling Logging::retireSink().
 * Even after retiring a sink, backends can keep using it safely until they
 * don't it any longer. There are two default Sinks connected to stdout and
 * stderr, which you can obtain via Logging::cout() and Logging::cerr().
 *
 * - The configuration of the logging system is handled by calling static
 * functions on the central class Logging. This class contains the central
 * registries for sinks and logger backends. It must be initialized by calling
 * Logging::init() before using any part of the logging system. During
 *     initialization, it parses a ParameterTree and creates sinks and backends
 * according to the given configuration. It is also possible to register new
 * types of sinks with this configuration system. Moreover, Logging also lets
 * you centrally disable all console output, which is very convenient for
 *     parallel programs. The system mutes itself by default for all parallel
 * ranks except 0.
 *
 * # Usage
 *
 * - Before you can use the logging system, you first need to initialize it:
 * ~~~
 * auto& parallel_helper = ...;
 * Dune::ParameterTree params = ...;
 * Dune::PDELab::Logging::init(parallel_helper.getCollectiveCommunication(),params.sub("logging"));
 * ~~~
 *
 * - In your normal program (outside of the namespace Dune), you need to import
 * the namespace Dune::Literals to make the user-defined literal `_fmt`
 * available:
 * ~~~
 * using namespace Dune::Literals;
 * ~~~
 *     This is not necessary if your code is in namespace Dune or a nested
 * namespace.
 *
 * - Now you can obtain the default Logger and log some messages:
 * ~~~
 * auto log = Dune::PDELab::logger();
 * log.notice("Hello, {}! The answer is {}."_fmt,"world",42);
 * log(Dune::PDELab::LogLevel::warning,"This seems complicated"_fmt);
 * ~~~
 *     \note The logging system uses the [{fmt}](https://fmtlib.net) library for
 * formatting the log message. This system works like printf(): You first supply
 * a format string, and placeholders in the format string are then replaced by
 * textual representions of the remaining arguments. {fmt} has a powerful syntax
 * for controlling the exact formatting of the arguments, see
 *           http://fmtlib.net/latest/syntax.html for details.
 *     \note The format string must be specially marked with the user-defined
 * literal `_fmt`, otherwise your code will not compile.
 *
 * - If you write a component that can be configured with a ParameterTree, you
 * can let the logging system configure a logger for you according to some
 * configuration keys (see Logging::logger(const Dune::ParameterTree&)):
 * ~~~
 * class MyComponent
 * {
 *   Dune::PDELab::Logger _log;
 *   ...
 * public:
 *   MyComponent(const Dune::ParameterTree& config)
 *     : _log(Dune::PDELab::logger(config))
 *   {}
 * };
 * ~~~
 *
 * # Configuration
 *
 * While it is possible to configure the logging system using function calls,
 * the recommended configuration style is via the ParameterTree, most likely
 * driven by an INI file.
 *
 * Logging::init() expects to be passed a ParameterTree subtree that only
 * contains logging-related configuration. You can configure sinks, backends and
 * global settings through this ParameterTree.
 *
 * * Sinks can be configured by creating a subtree with the prefix "sink." and
 * putting the sink-specific configuration into that subtree, e.g.:
 * ~~~{.ini}
 * [logging.sink.logfile]
 * type = file
 * file = production-{}.log
 * ~~~
 *     This configuration will create a sink named "logfile" that writes to the
 * file "production-X.log", where X will be replaced by the MPI rank of the
 * process, zero-padded on the left as necessary.
 *
 * * Backends can be configured by placing information under the prefix
 *     "backend.":
 *     * If you want to create a backend with the same configuration as the default
 *         backend, just a different name, create a key with value "default". The new
 *         backend will be named after the key:
 * ~~~{.ini}
 * [logging.backend]
 * newton = default
 * timestep = default
 * ~~~
 *     * If you want to apply custom configuration to the new backend, create a subtree
 *         named after the new backend and place the configuration keys for the backend
 *         in the subtree:
 * ~~~{.ini}
 * [logging.backend.timestep]
 * level = debug   # default level of loggers attached to this backend
 * sinks = logfile # only write to the logfile, not to the console
 * ~~~
 *     It is possible to combine these two approaches.
 *
 * * You can also configure the components that are automatically created by the system.
 *     In order to change the configuration of one of the default sinks stdout or stderr,
 *     place the parameters in the corresponding entries unter "sink.":
 * ~~~{.ini}
 * [logging.sink.stdout]
 * level = info
 *
 * [logging.sink.stderr]
 * level = waring
 * pattern = "STDERR [{backend}] {msg}"
 * ~~~
 *     The configuration for the default backend is in the subtree "default" at the root of
 *     the ParameterTree, not in "backend.default":
 * ~~~{.ini}
 * [logging]
 * default.level = info
 * default.sinks = stderr, logfile
 * ~~~
 *
 * ## Sink Configuration
 *
 * All sinks support the following configuration keys:
 *
 * | Key         | Description                                                                      |
 * |-------------|----------------------------------------------------------------------------------|
 * | type        | The type of the sink, selects the factory used to create the sink.               |
 * | level       | The maximum LogLevel processed, messages with a higher level are ignored.        |
 *
 * When configuring a sink, you **must** specify the "type" of the sink. Currently, supported
 * types include
 *
 * | Name          | Description                                                                    |
 * |---------------|--------------------------------------------------------------------------------|
 * | null          | A NullSink that ignores all input.                                             |
 * | file-per-rank | A FileSink that logs to a different file per rank.                             |
 * | rank-0-file   | A FileSink that logs to a file on rank 0 and ignores input everywhere else.    |
 *
 * See the documentation of the created class for supported type-specific configuration keys.
 *
 * Most sinks inherit from PatternFormatSink and support its pattern-based formatting of the log
 * output. For these sinks, you can set the configuration key "pattern" to a {fmt} format string
 * that supports the following named arguments:
 *
 * | Parameter   |Description                                                      |Notes         |
 * |-------------|-----------------------------------------------------------------|--------------|
 * | msg         | The message submitted by the user                               |              |
 * | level       | The log level of the message                                    |              |
 * | paddedlevel | The log level of the message, right-padded to the longest level |              |
 * | reltime     | The relative time since program start                           |              |
 * | reldays     | The number of full days since program start                     |              |
 * | abstime     | The absolute system time                                        | expensive    |
 * | backend     | The name of the backend used to log the message                 | right-padded |
 * | sink        | The name of the sink currently processing the message           |              |
 *
 * All of these parameters have the types returned by the corresponding member functions of
 * LogMessage, exept for the log levels, which are converted to a string representation, and the
 * logger names, which are right-padded to the width of the longest logger name.
 *
 * If your format string has to contain literal "{" or "}", escape them by doubling to "{{" or
 * "}}".
 *
 * You can employ additional formatting with the standard {fmt} format specification language. This
 * is especially important for "reltime" and "abstime".
 *
 * If the input format string lacks a trailing newline character, the sink will append one.
 *
 * For example, the pattern "[{reldays:0>2}-{reltime:12%T}] [{logger}] {msg}" causes messages to be
 * logged like
 * ~~~{.txt}
 * [02-08:32:51.941] [default] This is the actual message
 * ~~~
 *
 * ## Backend Configuration
 *
 * Backends support the following configuration keys;
 *
 * | Key         | Description                                                                      |
 * |-------------|----------------------------------------------------------------------------------|
 * | level       | The default maximum LogLevel of attached Loggers.                                |
 * | indent      | The default indentation of attached Loggers.                                     |
 * | enabled     | Flag for globally enabling or disabling all Loggers attached to this backend.    |
 * | sinks       | Comma-separated list of sinks for this backend. Replaces the default sinks.      |
 * | extra_sinks | Comma-separated list of sinks for this backend in addition to the default sinks. |
 *
 * ## Logger Configuration
 *
 * See Logging::logger(const ParameterTree&) for a list of configuration keys supported by Loggers.
 *
 * ## Global Configuration
 *
 * The logging system itself supports the following configuration keys:
 *
 * | Key         | Description                                                                      |
 * |-------------|----------------------------------------------------------------------------------|
 * | muted       | Overrides the default behavior of muting all MPI ranks with `comm.rank() > 0`    |
 *
 * \{
 */

  //! Central configuration of the logging system.
  class Logging {

#ifndef DOXYGEN

    // Internal class that stores all sink facories
    struct SinkFactoryRepository;

    // Internal data structure that stores the global state of the logging system
    struct State;

#endif // DOXYGEN

  public:

#if HAVE_MPI
    using CollectiveCommunication = Dune::CollectiveCommunication<MPI_Comm>;
#else
    using CollectiveCommunication = Dune::CollectiveCommunication<No_Comm>;
#endif

    //! Type-erased holder for a sink factory;
    using SinkFactory = std::function<std::shared_ptr<Sink>(std::string_view,LogLevel,std::size_t,const ParameterTree&)>;

  private:


    // Create initial sink factory repository and register default sinks
    static SinkFactoryRepository makeSinkFactoryRepository();

    // Access the sink factory repository
    static SinkFactoryRepository& sinkFactoryRepository();

    // Get sink factory with given name
    static SinkFactory& sinkFactory(const std::string& name);

    // Access the internal logging system state
    static State& state(const CollectiveCommunication* = nullptr);

    // Get a reference to the logger backend with given name
    static LoggerBackend& backend(std::string_view name);

  public:

    /**
     * \name Initialization
     *
     * Before the logging system can be used, its central configuration must be initialized.
     *
     * \{
     */


    //! Initializes the logging system.
    /**
     * This function must be called before any interaction with the logging system, except for registering new
     * `SinkFactory`s. As part of the setup, the passed-in ParameterTree will be used to populate the logging system
     * with sinks and backends.
     *
     * The logging system will automatically mute() itself if its rank in the passed in CollectiveCommunication object
     * is not zero.
     *
     * See the general description of the logging system for an explanation of the accepted keys in the ParameterTree.
     *
     * \param comm    The collective communication used by MPI-aware parts of the logging system.
     * \param params  Parameters for externally driven setup of the logging system.
     *
     */
    static void init(const CollectiveCommunication& comm, const ParameterTree& params = {});

    //! Registers a new SinkFactory.
    /**
     * Registering a new SinkFactory makes sinks of the produced type available to the ParameterTree-driven
     * configuration via the "type" parameter of a new sink.
     *
     * \note This function may be called before init() to make additional sink types available in the configuration
     *       parsed by init().
     */
    static void registerSinkFactory(const std::string& name, SinkFactory sink_factory);

    /**
     * \}
     */

    /**
     * \name Sinks
     *
     * Functionality for working with `Sink`s.
     *
     * \{
     */

    //! Makes a new sink with the given name and parameters.
    /**
     * This function will invoke the sink factory given by the key "type" in params and store the result under the given
     * name in the sink registry.
     *
     * If the name is already taken, the function raises an exception.
     *
     * \returns A shared_ptr to the newly created sink.
     */
    static std::shared_ptr<Sink> makeSink(const std::string& name,const ParameterTree& params);

    //! Registers a sink with the registry.
    /**
     * This function will store the given sink under the given name in the global sink registry.
     *
     * If the name is already taken, the function raises an exception.
     *
     */
    static void registerSink(std::shared_ptr<Sink> sink);

    //! Returns the sink with the given name.
    /**
     * As the function returns a shared_ptr, users of this function can safely continue using the sink even after it has
     * been retired from the registry.
     *
     * If there is no sink with the given name, the function raises an exception.
     *
     */
    static std::shared_ptr<Sink> sink(const std::string& name);

    //! Retires the named sink from the sink registry.
    /**
     * A retired sink cannot be used in configuring new backends anymore, but any backend already using the sink
     * can continue to do so. The sink only gets destroyed when its last user goes away.
     *
     * \note The logging system will not update the widestLogger() property of retired sinks, which can lead to
     *       format inconsistencies if backends with longer names are added after retiring the sink.
     */
    static bool retireSink(std::string_view sink);

    //! Returns the ConsoleSink connected to stdout.
    static std::shared_ptr<ConsoleSink> cout();

    //! Returns the ConsoleSink connected to stderr.
    static std::shared_ptr<ConsoleSink> cerr();

    /**
     * \}
     */

    /**
     * \name Loggers
     *
     * Functionality for obtaining `Logger`s from the system.
     *
     * \{
     */

    //! Returns a logger connected to the default backend.
    static Logger logger();

    //! Returns a logger connected to the given backend.
    /**
     * This function throws an exception if the backend could not be found.
     */
    static Logger logger(std::string_view backend);

    //! Returns a logger configured according to the ParameterTree.
    /**
     * This function inspects the given ParameterTree for its configuration. It understands the following keys:
     *
     * | Key         | Description                                                                |
     * |-------------|----------------------------------------------------------------------------|
     * | log.backend | The backend for the new logger. Uses the default backend if not specified. |
     * | log.level   | The maximum LogLevel that this logger will forward to the backend.         |
     * | log.indent  | The default indentation for messages logged with this logger.              |
     *
     * When not given, the "level" and "indent" parameters are set to the default values of the backend. If there is a
     * problem with any of the parameters, this function will throw an exception.
     */
    static Logger logger(const Dune::ParameterTree& params);

    //! Returns a logger configured according to the ParameterTree with a non-default fallback
    /**
     * This function is mostly intended for components that want to default to a named logger
     * without forcing the user to create the corresponding backend. The function first tries to get
     * a logger according to the configuration in the ParameterTree; if that configuration does not
     * name a backend, the system tries to get the backed named in `preferred`, and if that does not
     * exist either, it falls back to the default backend.
     */
    static Logger tryLogger(const Dune::ParameterTree& params, std::string_view preferred);

    /**
     * \}
     */

    /**
     * \name Backends
     *
     * Functionality for creating and configuring backends.
     *
     * \{
     */

    //! Registers a new backend with the logging system.
    /**
     * \param name                  The name of the backend.
     * \param level                 The default level of new loggers created for this backend.
     * \param attach_default_sinks  If this parameter is true, the backend will be connected to the same sinks as the
     *                              default backen, otherwise it will not be connected to any sinks.
     *
     * \returns A logger with default configuration attached to the new backend.
     */
    static Logger registerBackend(
      std::string_view name,
      LogLevel level,
      bool attach_default_sinks = true
      );

    //! Attaches a new sink to a backend.
    /**
     * Both the backend and the sink must be registered in the global logging system registry.
     *
     * \return  true if the sink was added to the backend, false if it was already added and nothing happenend.
     */
    static bool attachSink(std::string_view backend, std::string_view sink);

    //! Removes a sink from a backend.
    /**
     * Both the backend and the sink must be registered in the global logging system registry.
     *
     * \return  true if the sink was removed from the backend, false if it was not connected to the sink.
     */
    static bool detachSink(std::string_view backend, std::string_view sink);

    /**
     * \}
     */

    /**
     * \name Muting and unmuting the logging system
     *
     * The logging system can be centrally muted. In a muted logging system, all output to the standard console sinks is
     * disabled. This is a convenient way of avoiding multiple ranks in a parallel program writing over each other.
     *
     * \{
     */

    //! Returns whether the logging system has been muted.
    static bool muted();

    //! Mutes the logging system.
    static void mute();

    //! Unmutes the logging system.
    static void unmute();

    /**
     * \}
     */

    /**
     * \name Utility Functions
     *
     * Assorted utility functions and helpers.
     *
     * \{
     */

    //! Returns the collective communication object used by the logging system.
    static const CollectiveCommunication& comm();

    //! Returns whether the given name is a valid name for a logging system component.
    /**
     * The logging system enforces the following rules for its components:
     *
     * - All names consist of lowercase letters, digits, and the symbol '-'.
     * - All names must start with a letter.
     * - There cannot be more than one consecutive occurence of '-'.
     */
    static bool isValidName(std::string_view name);

    /**
     * \}
     */


  };


  //! \copydoc Logging::logger()
  inline Logger logger()
  {
    return Logging::logger();
  }

  //! \copydoc Logging::logger(std::string_view)
  inline Logger logger(std::string_view name)
  {
    return Logging::logger(name);
  }

  //! \copydoc Logging::logger(const Dune::ParameterTree&)
  inline Logger logger(const Dune::ParameterTree& params)
  {
    return Logging::logger(params);
  }

  //! Logs a message to the default logger.
  template<typename... Args>
  void log(Args&&... args)
  {
    logger()(std::forward<Args>(args)...);
  }

  /**
   * \}
   */

} // end namespace Dune::PDELab


#endif // DUNE_PDELAB_LOGGING_HH
