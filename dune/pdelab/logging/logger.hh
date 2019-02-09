#ifndef DUNE_PDELAB_LOGGING_LOGGER_HH
#define DUNE_PDELAB_LOGGING_LOGGER_HH

#include <cstdio>
#include <cstring>
#include <iostream>
#include <memory>
#include <optional>
#include <string_view>
#include <unordered_map>
#include <unordered_set>

#include <dune/pdelab/common/checks.hh>
#include <dune/pdelab/common/typetraits.hh>
#include <dune/pdelab/logging/fmt.hh>
#include <dune/pdelab/logging/sink.hh>

namespace Dune::PDELab {

  /**
   * \addtogroup logging
   * \{
   */

  class LoggerBackend;

  //! A Logger is a lightweight value object used to submit messages to the logging system.
  /**
   * A Logger is the main interaction point between user code and the logging system.
   *
   * * The Logger class has value semantics. You should **never** store a Logger by reference;
   *      copying it is extremely cheap. Due to the value semantics, Logger is
   *      default-constructible, but before you can use a default-constructed Logger, you must
   *      assign to it from a valid Logger.
   *
   * * A Logger instance is attached to a LoggerBackend, which controls its global state, like the
   *     list of attached `Sink`s and its name. A LoggerBackend can be globally disabled for all
   *     attached Loggers, which will completely turn off all attached Loggers. You cannot use a
   *     Logger that is not currently attached.
   *
   * * Each Logger instance carries its own trigger level, which limits which messages will be
   *     logged, and its own default indentation. The default values for these are taken from the
   *     LoggerBackend when creating the Logger, but changes to these values are local to the
   *     current Logger.
   *
   * In order to use the Logger, you call one of its logging methods. The generic logging methods,
   * which are not named after a LogLevel, take the following arguments:
   *
   * * An optional LogLevel.
   * * An optional indentation, which will be added to the Logger's default indentation.
   * * A format string, followed by the arguments referenced in the format string.
   *
   * \note The format string **must** be a compile-time literal (you cannot pass a string that you
   *       have constructed at runtime). Moreover, this fixed string must be followed by the
   *       user-defined literal **_fmt**. If your code does not live in namespace Dune, you must
   *       import this literal by placing
   * ~~~
   * using namespace Dune::Literals;
   * ~~~
   *       in your code.
   *
   * \warning Passing a normal string to the logging functions will result in a compile error!
   *
   * When you enable the macro `DUNE_PDELAB_CHECK_FORMAT_STRINGS`, the {fmt} library will check the
   * syntax of your format string at compile time and will make sure that you have provided all the
   * arguments mentioned in the format string. This is enabled by default in debug builds.
   *
   * Instead of the generic logging methods, you can (and mostly should) use the methods named after
   * the LogLevel that they will log to. Note that while the signature of these functions looks very
   * generic, they actually accept the same arguments as the generic versions described above, minus
   * the optional LogLevel, which is now given through the function name.
   */
  class Logger
  {

    friend class Logging;

  public:

    /**
     * \name Logging
     * Methods for submitting log messages
     *
     * Calling the logging methods is only allowed if `attached() == true`, otherwise you either get
     * an exception or undefined behavior, depending on DUNE_PDELAB_CHECK_LOGGER.
     *
     * \{
     */

    //! Logs the given log message with the default level of the Logger.
    template<typename... Args>
    void operator()(format_string_view format, Args&&... args)
    {
      if (_default_level <= _level) {
        handle(LogLevel::notice,_indent,format,fmt::make_format_args(std::forward<Args>(args)...));
      }
    }

    //! Logs the given log message.
    template<typename... Args>
    void operator()(LogLevel level, format_string_view format, Args&&... args)
    {
      if (level <= _level) {
        handle(level,_indent,format,fmt::make_format_args(std::forward<Args>(args)...));
      }
    }

    //! Logs the given log message and requests additional indentation that will be added to the Logger's default indentation.
    template<typename... Args>
    void operator()(LogLevel level, int indent, format_string_view format, Args&&... args)
    {
      if (level <= _level) {
        handle(level,_indent + indent,format,fmt::make_format_args(std::forward<Args>(args)...));
      }
    }

#ifndef DOXYGEN

    // The following methods are just alternative versions of the ones above, but these get selected
    // when compile-time verification of the format string is enabled.

    template<typename FS, typename... Args>
    std::enable_if_t<is_format_string_v<FS>>
    operator()(FS format, Args&&... args)
    {
      // Instantiate the code that performs the compile-time verification of the format message,
      // but never run it.
      if (false)
        fmt::format(format,std::forward<Args>(args)...);
      if (_default_level <= _level)
      {
        std::string_view raw_format(format);
        handle(_default_level,_indent,raw_format,fmt::make_format_args(std::forward<Args>(args)...));
      }
    }

    template<typename FS, typename... Args>
    std::enable_if_t<is_format_string_v<FS>>
    operator()(LogLevel level, FS format, Args&&... args)
    {
      // Instantiate the code that performs the compile-time verification of the format message,
      // but never run it.
      if (false)
        fmt::format(format,std::forward<Args>(args)...);
      if (level <= _level)
      {
        std::string_view raw_format(format);
        handle(level,_indent,raw_format,fmt::make_format_args(std::forward<Args>(args)...));
      }
    }

    template<typename FS, typename... Args>
    std::enable_if_t<is_format_string_v<FS>>
    operator()(LogLevel level, int indent, FS format, Args&&... args)
    {
      // Instantiate the code that performs the compile-time verification of the format message,
      // but never run it.
      if (false)
        fmt::format(format,std::forward<Args>(args)...);
      if (level <= _level)
      {
        std::string_view raw_format(format);
        handle(level,_indent + indent,raw_format,fmt::make_format_args(std::forward<Args>(args)...));
      }
    }


    // the following methods are fallback versions that get selected when the user forgot to append
    // _fmt to the format string. Instead of failing with an obscure error message, they tell the
    // user how to fix the problem.

    template<typename FS, typename... Args>
    std::enable_if_t<not is_format_string_v<FS>>
    operator()(FS format, Args&&... args)
    {
      static_assert(not Std::to_true_type_v<FS>,"You need to tag your format string with the user defined literal suffix _fmt");
    }

    template<typename FS, typename... Args>
    std::enable_if_t<not is_format_string_v<FS>>
    operator()(LogLevel level, FS format, Args&&... args)
    {
      static_assert(not Std::to_true_type_v<FS>,"You need to tag your format string with the user defined literal suffix _fmt");     }

    template<typename FS, typename... Args>
    std::enable_if_t<not is_format_string_v<FS>>
    operator()(LogLevel level, int indent, FS format, Args&&... args)
    {
      static_assert(not Std::to_true_type_v<FS>,"You need to tag your format string with the user defined literal suffix _fmt");     }

    #endif // DOXYGEN

    //! Logs the given message at LogLevel::critical.
    template<typename... Args>
    void critical(Args&&... args)
    {
      (*this)(LogLevel::critical,std::forward<Args>(args)...);
    }

    //! Logs the given message at LogLevel::error.
    template<typename... Args>
    void error(Args&&... args)
    {
      (*this)(LogLevel::error,std::forward<Args>(args)...);
    }

    //! Logs the given message at LogLevel::warning.
    template<typename... Args>
    void warning(Args&&... args)
    {
      (*this)(LogLevel::warning,std::forward<Args>(args)...);
    }

    //! Logs the given message at LogLevel::notice.
    template<typename... Args>
    void notice(Args&&... args)
    {
      (*this)(LogLevel::notice,std::forward<Args>(args)...);
    }

    //! Logs the given message at LogLevel::info.
    template<typename... Args>
    void info(Args&&... args)
    {
      (*this)(LogLevel::info,std::forward<Args>(args)...);
    }

    //! Logs the given message at LogLevel::detail.
    template<typename... Args>
    void detail(Args&&... args)
    {
      (*this)(LogLevel::detail,std::forward<Args>(args)...);
    }

    //! Logs the given message at LogLevel::debug.
    template<typename... Args>
    void debug(Args&&... args)
    {
      (*this)(LogLevel::debug,std::forward<Args>(args)...);
    }

    //! Logs the given message at LogLevel::trace.
    template<typename... Args>
    void trace(Args&&... args)
    {
      (*this)(LogLevel::trace,std::forward<Args>(args)...);
    }

    /**
     * \}
     */

    /**
     * \name Observers
     * Functions for inspection the Logger's state
     * \{
     */

    //! Returns the maximum log level at which this Logger will actually process log messages.
    LogLevel level() const
    {
      return _level;
    }

    //! Returns the default log level at which this Logger will log messages.
    LogLevel defaultLevel() const
    {
      return _default_level;
    }

    //! Returns the default indentation of messages logged with this Logger.
    int indent() const
    {
      return _indent;
    }

    //! Returns whether this Logger is currently attached to a LoggerBackend.
    bool attached() const
    {
      return _backend;
    }

    //! Returns the name of the LoggerBackend that this Logger is attached to.
    /**
     * \note This method will throw a LoggingError if `attached() == false`.
     */
    std::string_view name() const;

    //! Returns whether the LoggerBackend that this logger is attached to is enabled.
    /**
     * \note This method will throw a LoggingError if `attached() == false`.
     */
    bool backendEnabled() const;

    /**
     * \} observers
     */


    /**
     * \name Modifiers
     * Functions for changing the Logger's state.
     * \{
     */

    //! Increases the default indentation of this logger by the given value.
    void indent(int additional_indent)
    {
      DUNE_PDELAB_CHECK_LOGGER(_indent + additional_indent >= 0);
      _indent += additional_indent;
    }

    //! Creates a new logger with additional indentation
    Logger indented(int additional_indent)
    {
      DUNE_PDELAB_CHECK_LOGGER(_indent + additional_indent >= 0);
      return {*_backend,_level,_indent + additional_indent};
    }

    //! Sets the maximum log level at which this Logger will actually process log messages.
    void setLevel(LogLevel level)
    {
      _level = level;
    }

    //! Sets tthe default log level at which this Logger will log messages.
    void setDefaultLevel(LogLevel default_level)
    {
      _default_level = default_level;
    }

    //! Sets the default indentation of messages logged with this Logger.
    void setIndent(int indent)
    {
      DUNE_PDELAB_CHECK_LOGGER(indent >= 0);
      _indent = indent;
    }

    //! Enables the backend that this logger is attached to.
    /**
     * \note This method will throw a LoggingError if `attached() == false`.
     */
    void enableBackend();

    //! Disables the backend that this logger is attached to.
    /**
     * \note This method will throw a LoggingError if `attached() == false`.
     */
    void disableBackend();

    /**
     * \}
     */


    /**
     * \name Constructors
     * \{
     */

    //! Constructs an empty logger, which cannot be used before initializing it.
    Logger() = default;

    //! Creates a new logger that shares the backend of the given logger.
    /**
     * \param logger  The logger whose backend the new logger will share.
     * \param indent  Value added to the default indentation of the passed-in logger.
     */
    Logger(const Logger& logger, int indent)
      : _level(logger._level)
      , _indent(logger._indent + indent)
      , _backend(logger._backend)
    {}

    /**
     * \}
     */

  private:

    //! Internal constructor from a backend.
    Logger(LoggerBackend& backend, LogLevel level, int indent, LogLevel default_level = LogLevel::notice)
      : _level(level)
      , _indent(indent)
      , _backend(&backend)
      , _default_level(default_level)
    {}

    //! Method for handling the actual logging.
    /**
     * This method should never be inlined, as it might expand to a rather large amount of code that
     * we don't want to duplicate at every log site, as this would unnecessarily bloat the code and
     * make disabled log calls very expensive.
     */
    void handle(LogLevel level, int indent, std::string_view format, fmt::format_args args);

    LogLevel _level = LogLevel::all;
    int _indent = 0;
    LoggerBackend* _backend = nullptr;
    LogLevel _default_level = LogLevel::notice;

  };

  /**
   * \}
   */

} // end namespace Dune::PDELab

#endif // DUNE_PDELAB_LOGGING_LOGGER_HH
