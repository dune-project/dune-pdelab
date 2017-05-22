// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#if HAVE_CONFIG_H
#include "config.h"
#endif

#include <cstddef>
#include <iomanip>
#include <ios>
#include <ostream>
#include <memory>

#include <sys/types.h>
#include <unistd.h>

#if HAVE_MPI
#include <mpi.h>
#endif // HAVE_MPI

#include <dune/common/ios_state.hh>

#include <dune/pdelab/common/clock.hh>
#include <dune/pdelab/common/hostname.hh>

#include "logtag.hh"

namespace Dune {
  namespace PDELab {

    std::shared_ptr<LogtagFormatterBase>
    makeGeneralLogtagFormatter(std::ostream &(&formatFunc)(std::ostream&)) {
      return std::make_shared<
        GeneralLogtagFormatter<std::ostream &(*)(std::ostream&)>
        >(&formatFunc);
    }

    // anonymous namespace to avoid potential clashed with other compilation
    // units that might provide the name logtagFormatter
    namespace {
      std::shared_ptr<LogtagFormatterBase> &logtagFormatter() {
        static std::shared_ptr<LogtagFormatterBase> formatter =
          makeGeneralLogtagFormatter(hostPidWallUserLogtagFormatFunc);
        return formatter;
      }
    }

    const std::shared_ptr<LogtagFormatterBase> &getLogtagFormatter()
    { return logtagFormatter(); }
    void
    setLogtagFormatter(const std::shared_ptr<LogtagFormatterBase> &formatter) {
      if(!formatter)
        logtagFormatter() =
          makeGeneralLogtagFormatter(hostPidWallUserLogtagFormatFunc);
      else
        logtagFormatter() = formatter;
    }

    std::ostream &logtag(std::ostream &s) {
      Dune::ios_base_all_saver saver(s);
      getLogtagFormatter()->writeTag(s);
      return s;
    }

    WithLogtag::~WithLogtag() { setLogtagFormatter(savedFormatter); }

    namespace {
#if HAVE_MPI
      int &rank() {
        static int rank = -1;
        return rank;
      }
      std::size_t &rankwidth() {
        static std::size_t rankwidth = 1;
        return rankwidth;
      }
#endif // !HAVE_MPI
      std::size_t &hostwidth() {
        static std::size_t hostwidth = 1;
        return hostwidth;
      }
      std::size_t &pidwidth() {
        static std::size_t pidwidth = 1;
        return pidwidth;
      }
    } // anonymous namespace

    void logtagSetupMPI(bool syncWidthes) {
#if HAVE_MPI
      // setup rank
      MPI_Comm_rank(MPI_COMM_WORLD, &rank());
      if(syncWidthes) {
        // setup width of rank
        {
          int size;
          MPI_Comm_size(MPI_COMM_WORLD, &size);
          std::ostringstream s;
          s << size-1;
          rankwidth() = s.str().size();
        }
        // width of the host part
        {
          unsigned width = getHostName().size();
          MPI_Allreduce(MPI_IN_PLACE, &width, 1, MPI_UNSIGNED, MPI_MAX,
                        MPI_COMM_WORLD);
          hostwidth() = width;
        }
        // width of the pid part
        {
          std::ostringstream s;
          s << getpid();
          unsigned width = s.str().size();
          MPI_Allreduce(MPI_IN_PLACE, &width, 1, MPI_UNSIGNED, MPI_MAX,
                        MPI_COMM_WORLD);
          pidwidth() = width;
        }
      }
#endif // HAVE_MPI
    }

    namespace {
      std::ostream &writeHostName(std::ostream &s) {
        Dune::ios_base_all_saver saver(s);
        char fill = s.fill();
        try {
          s << std::setfill(' ') << std::setw(hostwidth()) << getHostName();
        }
        catch(...) { s.fill(fill); }
        s.fill(fill);
        return s;
      }

      std::ostream &writePid(std::ostream &s) {
        Dune::ios_base_all_saver saver(s);
        char fill = s.fill();
        try {
          s << std::setfill(' ') << std::setw(pidwidth())
            << std::setiosflags(std::ios_base::dec | std::ios_base::right)
            << getpid();
        }
        catch(...) { s.fill(fill); }
        s.fill(fill);
        return s;
      }

#if HAVE_MPI
      std::ostream &writeRank(std::ostream &s) {
        if(rank() < 0)
          return s << '?';
        Dune::ios_base_all_saver saver(s);
        char fill = s.fill();
        try {
          s << std::setfill(' ') << std::setw(rankwidth())
            << std::setiosflags(std::ios_base::dec | std::ios_base::right)
            << rank();
        }
        catch(...) { s.fill(fill); }
        s.fill(fill);
        return s;
      }
#endif

      void writeSeconds(std::ostream &s, TimeSpec seconds, std::size_t width) {
        Dune::ios_base_all_saver saver(s);
        char fill = s.fill();
        try {
          s << std::setfill(' ') << std::setw(width) << std::setprecision(6)
            << std::setiosflags(std::ios_base::dec | std::ios_base::fixed |
                                std::ios_base::right |
                                std::ios_base::showpoint)
            << seconds;
          s.fill(fill);
        }
        catch(...) { s.fill(fill); }
      }

    } // anonymous namespace

    std::ostream &hostRankWallUserLogtagFormatFunc(std::ostream &s) {
      s << std::setw(0) << "[h=" << writeHostName;
#if HAVE_MPI
      s << "|r=" << writeRank;
#endif // HAVE_MPI
      s << "|w=";
      writeSeconds(s, getWallTime(), 17);
      s << "|u=";
      writeSeconds(s, getProcessTime(), 12);
      s << "] ";
      return s;
    }

    std::ostream &hostPidWallUserLogtagFormatFunc(std::ostream &s) {
      s << std::setw(0) << "[h:p=" << writeHostName << ":" << writePid
        << "|w=";
      writeSeconds(s, getWallTime(), 17);
      s << "|u=";
      writeSeconds(s, getProcessTime(), 12);
      s << "] ";
      return s;
    }

    std::ostream &nullFormatFunc(std::ostream &s) { return s; }

  } // namespace PDELab
} // namespace Dune
