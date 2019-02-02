#include "config.h"

#include <cassert>
#include <string>
#include <string_view>
#include <vector>

#include <dune/pdelab/common/checks.hh>
#include <dune/pdelab/logging/loggerbackend.hh>
#include <dune/pdelab/logging/logger.hh>

namespace Dune::PDELab {

  void LoggerBackend::handle(const Logger& logger, LogLevel level, int indent, std::string_view format, fmt::format_args args)
  {
    if (not _enabled or _sinks.empty())
      return;

    auto time = LogMessage::Clock::now();

    _buffer.clear();
    assert(indent >= 0);
    assert(indent <= 40);

    using namespace std::literals;
    constexpr auto indent_template = "                                        "sv;
    _buffer.append(begin(indent_template),begin(indent_template) + indent);

    fmt::vformat_to(_buffer, format, args);

    auto msg = LogMessage(logger,level,indent,{_buffer.data(),_buffer.size()},time,time - _startup_time);

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
