#ifndef DUNE_PDELAB_LOGGING_FILESINKS_HH
#define DUNE_PDELAB_LOGGING_FILESINKS_HH

#include <iterator>
#include <fstream>

#include <dune/pdelab/logging/patternformatsink.hh>
#include <dune/pdelab/logging/fmt.hh>

namespace Dune::PDELab {

  class FileSink
    : public PatternFormatSink
  {

  public:

    FileSink(
      std::string_view name,
      LogLevel level,
      std::size_t widest_logger,
      const std::string& file_name,
      std::ios_base::openmode mode
      )
      : PatternFormatSink(name,level,widest_logger)
      , _file_name(file_name)
      , _file(file_name,mode)
    {}

    void process(const LogMessage& msg)
    {
      fmt::vprint(_file, fmt::to_string_view(pattern()), arguments(msg));
    }

  private:

    std::string _file_name;
    std::ofstream _file;

  };


} // namespace Dune::PDELab

#endif // DUNE_PDELAB_LOGGING_FILESINKS_HH
