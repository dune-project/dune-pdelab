#ifndef DUNE_PDELAB_LOGGING_CONSOLESINK_HH
#define DUNE_PDELAB_LOGGING_CONSOLESINK_HH

#include <cstdio>
#include <string_view>

#include <dune/pdelab/logging/sink.hh>

namespace Dune::PDELab {


  class ConsoleSink
    : public Sink
  {

  public:

    ConsoleSink(std::string_view name, std::FILE* stream, LogLevel level, std::size_t widest_logger)
      : Sink(name,level,widest_logger)
      , _stream(stream)
    {}

    void process(const LogMessage& msg) override;

  private:

    std::FILE* _stream;
    std::size_t _widest_logger = 0;
    std::string _format = "[{1}] {0}";
    fmt::basic_memory_buffer<char,200> _msg_buffer;
    fmt::basic_memory_buffer<char,20> _logger_buffer;
  };

} // end namespace Dune::PDELab

#endif // DUNE_PDELAB_LOGGING_CONSOLESINK_HH
