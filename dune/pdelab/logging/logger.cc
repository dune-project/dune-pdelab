#include "config.h"

#include <cassert>
#include <string>
#include <string_view>
#include <vector>

#include <dune/pdelab/common/checks.hh>
#include <dune/pdelab/logging/loggerbackend.hh>
#include <dune/pdelab/logging/logger.hh>

namespace Dune::PDELab {

  LoggerBackend::LoggerBackend(
    std::string_view name,
    LogMessage::Time startup_time,
    bool enabled,
    LogLevel default_level,
    int default_indent
    )
    : _enabled(enabled)
    , _startup_time(startup_time)
    , _default_level(default_level)
    , _default_indent(default_indent)
    , _name(name)
  {
    if (_default_indent < 0)
      DUNE_THROW(
        LoggingError,
        "Cannot create logger backend with negative indent " << _default_indent << ": " << _name
        );
  }

  void LoggerBackend::handle(const Logger& logger, LogLevel level, int indent, std::string_view format, fmt::format_args args)
  {
    if (not _enabled or _sinks.empty())
      return;

    assert(indent >= 0);
    assert(indent <= 40);

    auto time = LogMessage::Clock::now();

    fmt::basic_memory_buffer<char, 200> buffer;

    using namespace std::literals;
    constexpr auto indent_template = "                                        "sv;
    buffer.append(begin(indent_template),begin(indent_template) + indent);

    fmt::vformat_to(buffer, format, args);

    auto msg = LogMessage(logger,level,indent,{buffer.data(),buffer.size()},time,time - _startup_time);

    for (auto& sink : _sinks)
    {
      if (sink->level() >= level)
        sink->process(msg);
    }
  }

  bool Logger::backendEnabled() const
  {
    DUNE_PDELAB_CHECK_LOGGER(_backend && "This logger has no backend attached");
    return _backend->_enabled;
  }

  void Logger::enableBackend()
  {
    DUNE_PDELAB_CHECK_LOGGER(_backend && "This logger has no backend attached");
    _backend->_enabled = true;
  }

  void Logger::disableBackend()
  {
    DUNE_PDELAB_CHECK_LOGGER(_backend && "This logger has no backend attached");
    _backend->_enabled = false;
  }

  void Logger::handle(LogLevel level, int indent, std::string_view format, fmt::format_args args)
  {
    DUNE_PDELAB_CHECK_LOGGER(level > LogLevel::off);
    DUNE_PDELAB_CHECK_LOGGER(_backend && "This logger has no backend attached");
    _backend->handle(*this,level,indent,format,args);
  }

  std::string_view Logger::name() const
  {
    return _backend->_name;
  }

} // end namespace Dune::PDELab
