#include "config.h"

#include <cstring>

#include <dune/pdelab/common/exceptions.hh>
#include <dune/pdelab/logging/logmessage.hh>

namespace Dune::PDELab {

  std::string_view name(LogLevel level)
  {
    using namespace std::literals;
    switch (level)
    {
    case LogLevel::off:
      return "off"sv;
    case LogLevel::critical:
      return "critical"sv;
    case LogLevel::error:
      return "error"sv;
    case LogLevel::warning:
      return "warning"sv;
    case LogLevel::notice:
      return "notice"sv;
    case LogLevel::info:
      return "info"sv;
    case LogLevel::detail:
      return "detail"sv;
    case LogLevel::debug:
      return "debug"sv;
    case LogLevel::trace:
      return "trace"sv;
    case LogLevel::all:
      return "all"sv;
    default:
      DUNE_THROW(LoggingError,"Unknown log level: " << static_cast<int>(level));
    }
  }

  std::string_view paddedName(LogLevel level)
  {
    using namespace std::literals;
    switch (level)
    {
    case LogLevel::off:
      return "off     "sv;
    case LogLevel::critical:
      return "critical"sv;
    case LogLevel::error:
      return "error   "sv;
    case LogLevel::warning:
      return "warning "sv;
    case LogLevel::notice:
      return "notice  "sv;
    case LogLevel::info:
      return "info    "sv;
    case LogLevel::detail:
      return "detail  "sv;
    case LogLevel::debug:
      return "debug   "sv;
    case LogLevel::trace:
      return "trace   "sv;
    case LogLevel::all:
      return "all     "sv;
    default:
      DUNE_THROW(LoggingError,"Unknown log level: " << static_cast<int>(level));
    }
  }

  LogLevel parseLogLevel(std::string_view name)
  {
    using namespace std::literals;
    if (name == "off"sv)
      return LogLevel::off;
    if (name == "critical"sv)
      return LogLevel::critical;
    if (name == "error"sv)
      return LogLevel::error;
    if (name == "warning"sv)
      return LogLevel::warning;
    if (name == "notice"sv)
      return LogLevel::notice;
    if (name == "info"sv)
      return LogLevel::info;
    if (name == "detail"sv)
      return LogLevel::detail;
    if (name == "debug"sv)
      return LogLevel::debug;
    if (name == "trace"sv)
      return LogLevel::trace;
    if (name == "all"sv)
      return LogLevel::all;
    DUNE_THROW(LoggingError,"Cannot parse log level name: " << name);
  }

  const std::tm& LogMessage::localTime() const
  {
    if (not _local_time)
    {
      _local_time.emplace();
      std::time_t now = std::chrono::system_clock::to_time_t(time());
      std::tm* tm = std::localtime(&now);
      std::memcpy(&(*_local_time),tm,sizeof(std::tm));
    }
    return *_local_time;
  }

} // end namespace Dune::PDELab
