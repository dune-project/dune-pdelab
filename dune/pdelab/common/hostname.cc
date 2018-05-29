// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#if HAVE_CONFIG_H
#include "config.h"
#endif

#include <cstddef>
#include <cstring>
#include <string>
#include <vector>

#include <unistd.h>

#include "hostname.hh"

#ifdef __MINGW32__
#include "windows.h"
#endif

namespace Dune {
  namespace PDELab {

    //! C++ friendly wrapper around POSIX' gethostname() or GetComputerName() when compiling for Windows applications
    std::string getHostName() {

      #ifndef __MINGW32__
      std::size_t bufsize = 1024;
      std::vector<char> buf(bufsize);
      while(gethostname(&buf[0], buf.size()),
            buf.back() = '\0',
            std::strlen(&buf[0]) == buf.size()-1)
      {
        buf.clear();
        buf.resize(bufsize*=2);
      }
      #else
      long unsigned int bufsize = 1024;
      std::vector<char> buf(bufsize);
      while(!GetComputerName(&buf[0], &bufsize))
      {
        buf.clear();
        buf.resize(bufsize);
      }
      #endif

      // ignore everything after the first '.', if any
      std::vector<char>::iterator end = buf.begin();
      while(*end != '\0' && *end != '.') ++end;
      std::string str(buf.begin(), end);
      return str;
    }

  } // namespace PDELab
} // namespace Dune
