#ifndef DUNE_PDELAB_LOGGING_LOGGINGSTREAMBUFFER_HH
#define DUNE_PDELAB_LOGGING_LOGGINGSTREAMBUFFER_HH

#include <sstream>
#include <string>
#include <string_view>

#include <dune/pdelab/logging/logger.hh>

namespace Dune::PDELab {

  /**
   * \addtogroup logging
   * \{
   */

  //! An output-only std::streambuf that forwards to a Logger.
  /**
   * This buffer can be used to make standard C++ `std::ostream`s feed their data into the logging system.
   * It can operate in two different modes:
   *
   * - In line-buffered mode, the buffer will only output complete lines that have been terminated
   *   by a newline character. In particular, this mode will ignore any explicit flushing of the C++
   *   stream. This mode is capable of exactly reproducing the original output layout as designed by
   *   the user of the `std::ostream`, but messages may appear later than expected when the user
   *   explicitly flushes the C++ stream.
   *
   * - In unbuffered mode, the buffer will always forward all pending data everytime the C++ stream
   *   is flushed. As most logging sinks are line-oriented and insert an additional newline after
   *   each log message, this will not correctly reproduce the original layout of the output. As a
   *   lot of people use `std::endl` instead of just `"\n"` for ending their lines, this mode will
   *   not forward empty lines to the logging system to avoid empty lines after every regular line
   *   printed to the C++ stream.
   */
  class LoggingStreamBuffer
    : public std::stringbuf
  {

  public:

    //! Constructs a LoggingStreamBuffer without a working logger.
    LoggingStreamBuffer(bool line_buffered);

    //! Constructs a LoggingStreamBuffer.
    LoggingStreamBuffer(bool line_buffered, Logger stream_logger);

    //! Handles the log message generation.
    int sync() override;

    //! Returns a copy of the Logger used by this buffer.
    Logger logger() const
    {
      return _stream_log;
    }

    //! Sets the logger used by this buffer.
    void setLogger(Logger logger)
    {
      _stream_log = logger;
    }

    //! Returns whether this buffer is line-buffered.
    bool isLineBuffered() const
    {
      return _line_buffered;
    }

    //! Enavbles or disables the line-buffered mode of this buffer.
    void setLineBuffered(bool enabled)
    {
      _line_buffered = enabled;
    }

  private:

    Logger _stream_log;
    bool _line_buffered = true;
    bool _logging = false;

  };

} // namespace Dune::PDELab


#endif // DUNE_PDELAB_LOGGING_LOGGINGSTREAMBUFFER_HH
