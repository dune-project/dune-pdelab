#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string_view>


#include <dune/pdelab/common/exceptions.hh>
#include <dune/pdelab/logging/loggingstreambuffer.hh>

namespace Dune::PDELab {

  LoggingStreamBuffer::LoggingStreamBuffer(bool line_buffered)
    : std::stringbuf(std::ios::out)
    , _line_buffered(line_buffered)
  {}

  LoggingStreamBuffer::LoggingStreamBuffer(bool line_buffered, Logger stream_logger)
    : std::stringbuf(std::ios::out)
    , _stream_log(stream_logger)
    , _line_buffered(line_buffered)
  {}

  int LoggingStreamBuffer::sync()
  {
    // Get a view of the currently available data
    std::string_view buf(pbase(),pptr()-pbase());

    if (_logging)
    {
      std::fprintf(stderr,"\
======================================================================================================\n\
FATAL LOGGING ERROR: Recursive use of redirected stream detected during logging of a redirected stream\n\
======================================================================================================\n"
        );
      std::abort();
    }

    if (not buf.empty())
    {
      _logging = true;
      // Log one message per line of output
      std::string_view::size_type first = 0;
      for (auto last = buf.find_first_of('\n',first) ; last != buf.npos ; last = buf.find_first_of('\n',first))
      {
        // Ignore empty lines in unbuffered mode
        if (_line_buffered or first < last)
          _stream_log("{}"_fmt,buf.substr(first,last-first));

        // skip the current character, which is a newline
        first = last + 1;
      }

      if (_line_buffered)
      {
        if (first == buf.size())
        {
          // We have logged the entire buffer and can clear it
          seekpos(0,std::ios::out);
        }
        else if (first > 0 and buf.size() > 1024)
        {
          // There is still data in the buffer, and the buffer starts to grow a little large
          // Purge logged data from the buffer
          std::memmove(pbase(),pbase()+first,buf.size()-first);
          seekpos(buf.size()-first,std::ios::out);
        }
        // else do nothing
      }
      else
      {
        // Always print out incomplete lines as well
        if (first + 1 < buf.size() or buf[buf.size() - 1] != '\n' )
          _stream_log("{}"_fmt,buf.substr(first));

        // clear buffer
        seekpos(0,std::ios::out);
      }
      _logging = false;
    }
    // Forward to base class
    return std::stringbuf::sync();
  }

} // namespace Dune::PDELab
