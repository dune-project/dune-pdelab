#ifndef DUNE_PDELAB_LOGGING_LOGMESSAGE_HH
#define DUNE_PDELAB_LOGGING_LOGMESSAGE_HH

#include <chrono>
#include <ctime>
#include <optional>
#include <string_view>

#include <dune/pdelab/logging/fmt.hh>

namespace Dune::PDELab {

  /**
   * \addtogroup logging
   * \{
   */

  //! Severity level for log messages.
  enum class LogLevel
  {

    //! Designates that log messages are disabled.
    /**
     * \note This flag is only intended for setting the log level of a logging component, **never**
     *       use this level when logging a message!
     */
    off = 0,

    //! A critical error during program execution, forcing the program to abort.
    critical = 3,

    //! A recoverable error, from which the program might recover and continue.
    error = 6,

    //! A warning about some problem, but which does not necessarily require an immediate reaction.
    warning = 9,

    //! Information about the high-level program flow.
    notice = 12,

    //! Somewhat finer information about the program flow, still interesting in most cases.
    info = 15,

    //! Detailed information about the program flow, like per-iteration information in solvers.
    detail = 18,

    //! Debug-level information that is normally only interesting when debugging a program.
    debug = 21,

    //! Fine-grained debug information that should not be turned on for the whole program.
    /**
     * Log messages at this level could for example be used to trace the entry into and exit out of
     * a given function.
     */
    trace = 24,

    //! Designates that all log messages are enabled.
    /**
     * \note This flag is only intended for setting the log level of a logging component, **never**
     *       use this level when logging a message!
     */
    all = 30

  };

  //! Returns the name of a given log level.
  /**
   * \note This function will throw an exception of type LoggingError if level is not a valid
   *       LogLevel.
   */
  std::string_view name(LogLevel level);

  //! Parses the given string into the corresponding LogLevel.
  /**
   * \note This function will throw an exception of type LoggingError if name is not the name of a
   *       valid LogLevel.
   */
  LogLevel parseLogLevel(std::string_view name);

  class Logger;

  //! A log message as created by a Logger and passed on to its connected `Sink`s.
  class LogMessage
  {

    friend class LoggerBackend;

  public:

    //! The std::chrono clock used by the logging system.
    using Clock    = std::chrono::system_clock;

    //! The type used to represent the time point at which the message was logged.
    using Time     = Clock::time_point;

    //! The type used to represent the duration since program start.
    using Duration = Clock::duration;

    //! LogMessage cannot be copied.
    LogMessage(const LogMessage&) = delete;
    //! LogMessage cannot be copied.
    LogMessage& operator=(const LogMessage&) = delete;

  private:

    // LogMessage can only be constructed by LoggerBackend
    LogMessage(
      const Logger& logger,
      LogLevel level,
      int indent,
      std::string_view payload,
      Time time,
      Duration relative_time
      )
      : _payload(payload)
      , _level(level)
      , _indent(indent)
      , _logger(logger)
      , _time(time)
      , _relative_time(relative_time)
    {}

  public:

    //! Returns the formatted message logged by the user.
    std::string_view payload() const
    {
      return _payload;
    }

    //! Returns the log level specified by the user.
    LogLevel level() const
    {
      return _level;
    }

    //! Returns the requested indentation of the message.
    int indent() const
    {
      return _indent;
    }

    //! Returns a reference to the Logger that was used to log this message.
    const Logger& logger() const
    {
      return _logger;
    }

    //! Returns the absolute point in time at which this message was logged.
    Time time() const
    {
      return _time;
    }

    /**
     * \brief Returns the duration since program start (more precisely, since setup of the logging
     * system) at which this message was logged.
     */
    Duration relativeTime() const
    {
      return _relative_time;
    }

    //! Returns the sub-second part of time().
    std::chrono::nanoseconds absoluteSecondFraction() const
    {
      if (_absolute_second_fraction == Duration::zero())
      {
        auto seconds = std::chrono::floor<std::chrono::seconds>(_time);
        _absolute_second_fraction = _time - seconds;
      }
      return _absolute_second_fraction;
    }

    //! Returns the sub-second part of relativeTime().
    std::chrono::nanoseconds relativeSecondFraction() const
    {
      if (_relative_second_fraction == Duration::zero())
      {
        auto seconds = std::chrono::floor<std::chrono::seconds>(_relative_time);
        _relative_second_fraction = _relative_time - seconds;
      }
      return _relative_second_fraction;
    }

    //! Returns local time of the point in time  at which this message was logged.
    /**
     * \warning The first call to this method is rather expensive, avoid it if you don't need to
     *          know the absolute point im time at which log messages occured. It might be a better
     *          alternative to record the relative time since program start instead. The return
     *          value is cached internally, so subsequent calls to this method on a given LogMessage
     *          are cheap.
     */
    const std::tm& localTime() const;

  private:

    std::string_view _payload;
    LogLevel _level;
    int _indent;
    const Logger& _logger;

    Time _time;
    Duration _relative_time;

    mutable std::chrono::nanoseconds _absolute_second_fraction = Duration::zero();
    mutable std::chrono::nanoseconds _relative_second_fraction = Duration::zero();

    mutable std::optional<std::tm> _local_time;

  };

  /**
   * \} logging
   */

} // end namespace Dune::PDELab

#endif // DUNE_PDELAB_LOGGING_LOGMESSAGE_HH
