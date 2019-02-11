#include "config.h"

#include <iterator>
#include <tuple>

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

  // helper functions for building a fmt::format_arg_store from a tuple
  template<typename T, std::size_t... i>
  static auto _extract_format_args(T& args, std::index_sequence<i...>)
  {
    return fmt::make_format_args(std::get<i>(args)...);
  }

  // helper functions for building a fmt::format_arg_store from a tuple
  template<typename T>
  static auto extract_format_args(T& args)
  {
    return _extract_format_args(args,std::make_index_sequence<std::tuple_size_v<T>>{});
  }

  PatternFormatSink::Arguments::Arguments(const LogMessage& msg, const PatternFormatSink& sink)
  {
    // As fmt::format_args only stores a pointer to an array that must be held somewhere for the
    // duration of related {fmt} calls, and as the format args stored in that array only hold
    // references to the actual data, we must store both the objects returned by the accessors of
    // the log message etc. and a fmt::format_arg_store for the array of format arguments.


    // Extract and store data items in a tuple within an internal buffer.

    // We use lazy evaluation for the logger name and the local time, as those are rather expensive.
    // IMPORTANT: Update itemIndex() when making changes to this list!

    using Data = std::tuple<
      std::string_view,
      std::string_view,
      std::string_view,
      LogMessage::Duration,
      decltype(SinkMessageItems::relativeDays(msg)),
      decltype(SinkMessageItems::localTime(msg)),
      decltype(SinkMessageItems::backend(msg,sink.widestLogger())),
      std::string_view
      >;

    // Make sure our data actually fits into the buffer
    static_assert(
      sizeof(Data) <= ArgumentDataBufferSize,
      "Pattern format data does not fit into buffer"
      );

    auto data = new(_data_buffer) Data(
      msg.payload(),
      Dune::PDELab::name(msg.level()),
      paddedName(msg.level()),
      msg.relativeTime(),
      SinkMessageItems::relativeDays(msg),
      SinkMessageItems::localTime(msg),
      SinkMessageItems::backend(msg,sink.widestLogger()),
      sink.name()
      );

    // If Data is not trivially destructible, we have to arrange for its destructor to be called
    // in our destructor
    if (not std::is_trivially_destructible_v<Data>)
    {
      _cleanup = [data]()
      {
        data->~Data();
      };
    }

    // Make sure the format_arg_store actually fits into its bufffer
    static_assert(
      sizeof(decltype(extract_format_args(*data))) <= ArgumentArgsBufferSize,
     "Pattern format_arg_store does not fit into buffer"
      );

    // fmt::format_arg_store is always trivially destructible
    auto args = new(_args_buffer) auto(extract_format_args(*data));

    // store type-erased data
    _args = *args;

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
    if (item == "reldays")
      return 4;
    if (item == "abstime")
      return 5;
    if (item == "backend")
      return 6;
    if (item == "sink")
      return 7;
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

    enum class Transition { none, start, stop };

    auto transition = Transition::none;

    bool in_argument = false;
    bool in_name = false;
    bool is_number = true;

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


  void PatternFormatSink::setPatternFormatParameters(const ParameterTree& params)
  {
    if (params.hasKey("pattern"))
      setPattern(params["pattern"]);
  }


} // end namespace Dune::PDELab
