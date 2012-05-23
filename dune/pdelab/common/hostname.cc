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

namespace Dune {
  namespace PDELab {

    //! C++ friendly wrapper around POSIX' gethostname()
    std::string getHostName() {
      std::size_t bufsize = 1024;
      std::vector<char> buf(bufsize);
      while(gethostname(&buf[0], buf.size()),
            buf.back() = '\0',
            std::strlen(&buf[0]) == buf.size()-1)
      {
        buf.clear();
        buf.resize(bufsize*=2);
      }
      // ignore everything after the first '.', if any
      std::vector<char>::iterator end = buf.begin();
      while(*end != '\0' && *end != '.') ++end;
      std::string str(buf.begin(), end);
      return str;
    }

  } // namespace PDELab
} // namespace Dune
