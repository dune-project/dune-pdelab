#ifndef DUNE_PDELAB_LOGGING_LOGGERBACKEND_HH
#define DUNE_PDELAB_LOGGING_LOGGERBACKEND_HH

#include <memory>
#include <vector>

#include <dune/pdelab/logging/fmt.hh>
#include <dune/pdelab/logging/sink.hh>

namespace Dune::PDELab {

  /**
   * \addtogroup logging
   * \{
   */

  //! Internal backend class for loggers that holds their common state.
  /**
   * \note This class is an implementation detail of the logging system. While it is used by the
   *       system to provide certain properties and behaviors described in the general
   *       documentation, it can change at any time without further notice. For this reason, the
   *       entire class contents is `private`.
   */
  class LoggerBackend
  {

    friend class Logger;
    friend class Logging;

    // This method is defined in logger.cc to allow inlining into Logger::handle()
    void handle(const Logger& logger, LogLevel level, int indent, std::string_view format, fmt::format_args args);

    bool _enabled = true;
    LogMessage::Time _startup_time;
    std::vector<std::shared_ptr<Sink>> _sinks;
    LogLevel _default_level = LogLevel::notice;
    int _default_indent = 0;
    std::string_view _name;
    fmt::basic_memory_buffer<char, 200> _buffer;

    // This method is defined in logger.cc
    LoggerBackend(
      std::string_view name,
      LogMessage::Time startup_time,
      bool enabled,
      LogLevel default_level = LogLevel::notice,
      int default_indent = 0
      );

  };

  /**
   * \}
   */

} // namespace Dune::PDELab

#endif //  DUNE_PDELAB_LOGGING_LOGGERBACKEND_HH
