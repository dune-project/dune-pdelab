#ifndef DUNE_PDELAB_LOGGING_FILESINKS_HH
#define DUNE_PDELAB_LOGGING_FILESINKS_HH

#include <iterator>
#include <fstream>

#include <dune/pdelab/logging/patternformatsink.hh>
#include <dune/pdelab/logging/fmt.hh>

namespace Dune::PDELab {

  /**
   * \addtogroup logging
   * \{
   */

  //! Sink for logging data to a file.
  /**
   * FileSink opens a given file and logs all data to that file. It is typically configured by
   * specifying the sink type "file" in the logging configuration. It supports the following
   * configuration keys:
   *
   * | Key     | Description                                                  |
   * |-------- |--------------------------------------------------------------|
   * | file    | The name of the log file. See below for further information. |
   * | mode    | "truncate" or "append".                                      |
   * | pattern | The pattern for the log message, see PatternFormatSink.      |
   *
   * When running programs in parallel, each rank will write to its own file: the filename will
   * automatically be prepended with the 0-padded rank of each process. Alternatively, if the
   * filename contains a "{}" somewhere, this placeholder will be replaced by the 0-padded rank
   * value. In this case, the replacement also happens when running a sequential program.
   */
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

  /**
   * \}
   */

} // namespace Dune::PDELab

#endif // DUNE_PDELAB_LOGGING_FILESINKS_HH
