#ifndef DUNE_PDELAB_LOGGING_CONSOLESINK_HH
#define DUNE_PDELAB_LOGGING_CONSOLESINK_HH

#include <cstdio>
#include <string_view>

#include <dune/pdelab/logging/patternformatsink.hh>

namespace Dune::PDELab {

  /**
   * \addtogroup logging
   * \{
   */

  //! A sink for writing to the console, typically stdout or stderr.
  /**
   * The logging system wraps stdout and stderr in two instances of this class, which are available
   * as Logging::cout() and Logging::cerr();
   */
  class ConsoleSink
    : public PatternFormatSink
  {

  public:

    ConsoleSink(
      std::string_view name,
      std::FILE* stream,
      const std::string& pattern,
      LogLevel level,
      std::size_t widest_logger
      )
      : PatternFormatSink(name,level,widest_logger)
      , _stream(stream)
    {
      setPattern(pattern);
    }

    void process(const LogMessage& msg) override
    {
      fmt::vprint(_stream,pattern(),arguments(msg));
    }

  private:

    std::FILE* _stream;

  };

  /**
   * \}
   */

} // end namespace Dune::PDELab

#endif // DUNE_PDELAB_LOGGING_CONSOLESINK_HH
