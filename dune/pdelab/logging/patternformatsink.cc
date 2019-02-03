#include "config.h"

#include <iterator>

#include <dune/pdelab/common/exceptions.hh>
#include <dune/pdelab/logging/patternformatsink.hh>
#include <dune/pdelab/logging/sinkmessageitems.hh>

namespace Dune::PDELab {

  PatternFormatSink::PatternFormatSink(
      std::string_view name,
      LogLevel level,
      std::size_t widest_logger
    )
    : Sink(name,level,widest_logger)
  {
    setPattern("{reltime:9%M:%S} {msg}");
  }


  fmt::format_args PatternFormatSink::arguments(const LogMessage& msg) const
  {
    // just package up the promised arguments. We use lazy evaluation
    // for the logger name and the local time, as those are rather
    // expensive.

    // IMPORTANT: Update itemIndex() when making changes to this list!
    return fmt::make_format_args(
      msg.payload(),
      Dune::PDELab::name(msg.level()),
      paddedName(msg.level()),
      msg.relativeTime(),
      SinkMessageItems::localTime(msg),
      SinkMessageItems::backend(msg,_logger_buffer,widestLogger()),
      name()
      );
  }

  std::size_t PatternFormatSink::itemIndex(std::string_view item)
  {
    // Translates the documented keywords to the correct position in the
    // list of arguments in argumetns()
    if (item == "payload" or item == "msg")
      return 0;
    if (item == "level")
      return 1;
    if (item == "paddedlevel")
      return 2;
    if (item == "reltime")
      return 3;
    if (item == "abstime")
      return 4;
    if (item == "backend")
      return 5;
    if (item == "sink")
      return 6;
    DUNE_THROW(LoggingError,"Unknown log sink pattern item: " << item);
  }


  void PatternFormatSink::setPattern(const std::string& pattern)
  {

    // We need to transform the pattern and replace the named keywords with their numeric
    // equivalents. This parser is a little more elaborate than it should be, but we have to cope
    // with escaped { and }.

    fmt::basic_memory_buffer<char,50> buf;
    fmt::basic_memory_buffer<char,10> arg;

    constexpr char start = '{';
    constexpr char stop  = '}';
    constexpr char none  = '\0';

    enum class Transition { none, start, stop };

    auto transition = Transition::none;

    bool in_argument = false;
    bool in_name = false;
    bool is_number = true;
    bool first = true;

    for(auto c : pattern)
    {

      if (transition == Transition::start)
      {
        transition = Transition::none;
        if (c == start)
        {
          buf.push_back(c);
          continue;
        }
        if (in_argument)
          DUNE_THROW(LoggingError,"Invalid pattern: { inside arguments not supported");
        in_argument = true;
        in_name = true;
        is_number = true;
        first = true;
        arg.clear();
      }

      if (transition == Transition::stop)
      {
        transition = Transition::none;
        if (c == stop)
        {
          buf.push_back(c);
          continue;
        }
        if (in_name)
        {
          if (is_number)
            buf.append(arg.data(),arg.data() + arg.size());
          else
          {
            fmt::format_to(std::back_inserter(buf),"{}",itemIndex({arg.data(),arg.size()}));
          }
        }
        in_argument = false;
        in_name = false;
      }

      if (in_name)
      {
        if (c == ':' or c == stop)
        {
          if (is_number)
            buf.append(arg.data(),arg.data() + arg.size());
          else
          {
            fmt::format_to(std::back_inserter(buf),"{}",itemIndex({arg.data(),arg.size()}));
          }
          in_name = false;
        }
        else
        {
          if (not ('0' <= c and c <= '9'))
            is_number = false;
          arg.push_back(c);
          first = false;
          continue;
        }
      }

      if (c == start)
        transition = Transition::start;

      if (c == stop)
        transition = Transition::stop;

      buf.push_back(c);

    }

    // Append a newline character if the user didn't.
    if (buf[buf.size()-1] != '\n')
      buf.push_back('\n');

    _input_pattern = pattern;
    _pattern = std::string(begin(buf),end(buf));
  }


  void PatternFormatSink::setParameters(PatternFormatSink& sink, const ParameterTree& params)
  {
    if (params.hasKey("pattern"))
      sink.setPattern(params["pattern"]);
  }


} // end namespace Dune::PDELab
