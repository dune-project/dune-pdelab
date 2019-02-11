#ifndef DUNE_PDELAB_LOGGING_PATTERNFORMATSINK_HH
#define DUNE_PDELAB_LOGGING_PATTERNFORMATSINK_HH

#include <cstdio>
#include <functional>
#include <string_view>

#include <dune/common/parametertree.hh>

#include <dune/pdelab/logging/fmt.hh>
#include <dune/pdelab/logging/sink.hh>

namespace Dune::PDELab {

  /**
   * \addtogroup logging
   * \{
   */

  //! Base class for sinks that offer pattern-based message line formatting.
  /**
   * PatternFormatSink makes it easy to customize the logging format of a sink. It stores a {fmt}
   * pattern that will be used to format the log entry. In order to improve readability, this
   * pattern can access a number of named parameters:
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
   * If the input format string lacks a trailing newline character, the setPattern() method will
   * append it.
   *
   * For example, the pattern "[{reldays:0>2}-{reltime:12%T}] [{logger}] {msg}" causes messages to be
   * logged like
   * ~~~{.txt}
   * [02-08:32:51.941] [default] This is the actual message
   * ~~~
   *
   * You can set the pattern via the ParameterTree key "pattern".
   *
   * As a programmer, you use this class by inheriting from it. In your `process()` implementation,
   * create a named instance of `Arguments` and pass it to a version of `fmt::vformat()` or
   * `fmt::vprint()`, together with the pattern obtained by calling `pattern()`.
   */
  class PatternFormatSink
    : public Sink
  {

  protected:

    static constexpr std::size_t ArgumentDataBufferSize = 200;
    static constexpr std::size_t ArgumentArgsBufferSize = 200;

    //! Stores the required arguments for the pattern.
    /**
     * This class provides opaque storage for the format arguments required by
     * the current pattern. An object of this type should be constructed in
     * derived classes before formatting the pattern. This type implicitly casts
     * to `fmt::format_args`, which can be used passed to the {fmt} formatting
     * functions.
     *
     * \note You **must** create a named object of this type; trying to feed a
     *       temporary object to {}fmt} will cause a compiler error.
     */
    class Arguments
    {

      fmt::format_args _args;
      std::function<void()> _cleanup;
      char* _data_buffer[ArgumentDataBufferSize];
      char* _args_buffer[ArgumentArgsBufferSize];

    public:

      //! Create argument store for the given LogMessage and PatternFormatSink.
      Arguments(const LogMessage&, const PatternFormatSink&);

      //! Arguments cannot be copied or moved.
      Arguments(const Arguments&) = delete;

      //! Arguments cannot be copied or moved.
      Arguments& operator=(const Arguments&) = delete;

      //! Free any resources allocated by Arguments.
      ~Arguments()
      {
        if (_cleanup)
          _cleanup();
      }

      //! Named variables of type Arguments cast implicitly into `fmt::format_args`.
      operator fmt::format_args&() &
      {
        return _args;
      }

    };

    //! Constructs a PatternFormatSink with the default pattern "{reltime:9%M:%S} {msg}".
    PatternFormatSink(
      std::string_view name,
      LogLevel level,
      std::size_t widest_logger
      );

    //! Constructs a {fmt} argument list from the LogMessage that can be used to format the stored pattern.
    Arguments arguments(const LogMessage& msg) const;

  public:

    //! Returns the {fmt} format pattern for formatting the log messages for this sink.
    /**
     * This function returns the transformed pattern that can be used to actually format messages.
     * In this pattern, keywords like "reltime" have been replaced by the indices of the corresponding
     * data items.
     */
    const std::string& pattern() const
    {
      return _pattern;
    }

    //! Returns the original format pattern set by the user.
    /**
     * This function returns the original pattern, which is more readable but cannot be used to format
     * messages.
     */
    const std::string& userPattern() const
    {
      return _input_pattern;
    }

    //! Sets a new format pattern for this sink.
    /**
     * This method sets the given pattern as the new format string for this sink. The pattern is
     * transformed before storing it, any numeric keywords as described in the general class
     * documentation is replaced by the corresponding index as returned by itemIndex().
     *
     * This process will detect some, but not all syntax problems in the format string and throw
     * an exception in that case.
     */
    void setPattern(const std::string& pattern);

    //! Returns the index of a named log message item.
    static std::size_t itemIndex(std::string_view item);

    //! Parses the ParameterTree for applicable parameters and applies then to this sink.
    void setPatternFormatParameters(const ParameterTree& params);

  private:

    std::string _input_pattern;
    std::string _pattern;

  };

  /**
   * \}
   */

} // end namespace Dune::PDELab

#endif // DUNE_PDELAB_LOGGING_PATTERNFORMATSINK_HH
