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
   * \brief Logging of status and error messages to different targets with ParameterTree-based configuration.
   *
   * Some more weird text.
   * \{
   */

  //! Central configuration of the logging system.
  class Logging
  {

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
