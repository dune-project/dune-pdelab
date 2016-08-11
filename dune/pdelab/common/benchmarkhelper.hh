#ifndef DUNE_PDELAB_COMMON_BENCHMARKHELPER_HH
#define DUNE_PDELAB_COMMON_BENCHMARKHELPER_HH

#include <ostream>
#include <iomanip>
#include <map>
#include <string>
#include <vector>
#include <limits>
#include <cmath>

#include <dune/common/exceptions.hh>
#include <dune/common/ios_state.hh>

#include <ctime>

namespace Dune {
  namespace PDELab {

    struct CppClockWallTimeSource
    {

      double operator()() const
      {
        return static_cast<double>(std::clock()) / static_cast<double>(CLOCKS_PER_SEC);
      }

    };

#if HAVE_MPI
#include"mpi.h"

    struct MPIWallTimeSource
    {

      double operator()() const
      {
        return MPI_Wtime();
      }

    };

    typedef MPIWallTimeSource DefaultTimeSource;

#else

    typedef CppClockWallTimeSource DefaultTimeSource;

#endif

    struct Timing
    {
      double start;
      double end;

      double elapsed() const
      {
        return end - start;
      }

    };

    struct BenchmarkEntry
    {
      std::vector<Timing> timings;
      double min;
      double max;
      double avg;
      double std_dev;
      std::uint64_t operations;

      BenchmarkEntry()
        : operations(0)
      {}

    };

    template<typename TimeSource = DefaultTimeSource>
    struct BenchmarkHelper
    {

      BenchmarkHelper(std::string name, std::size_t max_runs = 1, TimeSource timeSource = TimeSource())
        : _name(name)
        , _time(timeSource)
        , _run(0)
        , _max_runs(max_runs)
        , _statistics_stale(true)
        , _threads(1)
        , _join_mode(JoinMode::undefined)
      {
        _run_times.timings.resize(max_runs);
      }


      void start_run()
      {
        if (_run >= _max_runs)
          {
            DUNE_THROW(Dune::RangeError,"maximum number of benchmark runs exceeded");
          }
        _statistics_stale = true;
        _run_times.timings[_run].start = _time();
      }

      void start_run(std::ostream& s)
      {
        start_run();
        ios_base_all_saver ios_saver(s);
        s << _name << " (" << std::setw(2) << _run << " of " << std::setw(2) << _max_runs << ") " << std::flush;
      }


      void end_run()
      {
        _run_times.timings[_run].end = _time();
        ++_run;
      }

      void end_run(std::ostream& s)
      {
        end_run();
        ios_base_all_saver ios_saver(s);
        s << " " << std::setw(10) << std::setprecision(3) << _run_times.timings[_run-1].elapsed() << " sec" << std::endl;
      }

      void start(std::string task)
      {
        std::pair<
          std::map<std::string,BenchmarkEntry>::iterator,
          bool
          > res = _tasks.insert(make_pair(task,BenchmarkEntry()));
        if (res.second)
          res.first->second.timings.resize(_max_runs);
        res.first->second.timings[_run].start = _time();
        _statistics_stale = true;
      }

      void start(std::string task, std::ostream& s)
      {
        start(task);
      }

      void end(std::string task)
      {
        _tasks[task].timings[_run].end = _time();
        _statistics_stale = true;
      }

      void setOperations(std::string task, std::uint64_t operations)
      {
        _tasks[task].operations = operations;
      }

      void end(std::string task, std::ostream& s)
      {
        end(task);
        s << "." << std::flush;
      }

      void update_entry(BenchmarkEntry& entry)
      {
        entry.min = std::numeric_limits<double>::max();
        entry.max = 0;
        entry.avg = 0;
        entry.std_dev = 0;

        for (std::vector<Timing>::iterator it = entry.timings.begin(), end = entry.timings.end();
             it != end;
             ++it)
          {
            const double elapsed = it->elapsed();
            entry.min = std::min(entry.min,elapsed);
            entry.max = std::max(entry.max,elapsed);
            entry.avg += elapsed;
            entry.std_dev += elapsed*elapsed;
          }

        entry.avg /= entry.timings.size();
        entry.std_dev /= entry.timings.size();
        entry.std_dev = std::sqrt(entry.std_dev - entry.avg*entry.avg);
      }

      void update_statistics()
      {
        if (!_statistics_stale)
          return;
        _max_name_len = 5; // strlen("total")
        for (std::map<std::string,BenchmarkEntry>::iterator it = _tasks.begin(), end = _tasks.end();
             it != end;
             ++it)
          {
            _max_name_len = std::max(_max_name_len,it->first.size());
            update_entry(it->second);
          }

        update_entry(_run_times);

        _statistics_stale = false;
      }

      double normalize_time(double t) const
      {
        switch (_join_mode)
          {
          case JoinMode::avg:
            return t / _threads;
          default:
            return t;
          }
      }

      void print_entry(std::ostream& s, std::string name, const BenchmarkEntry& entry, bool summary_only = false) const
      {
        s << std::setw(_max_name_len + 1) << std::left << name
          << std::right << std::scientific << std::setw(10) << std::setprecision(2);
        if (_threads == 1 && !summary_only)
          for (std::vector<Timing>::const_iterator it = entry.timings.begin(),
                 end = entry.timings.end();
               it != end;
               ++it)
            {
              s << std::setw(10) << it->elapsed();
            }
        s << std::setw(10) << normalize_time(entry.min)
          << std::setw(10) << normalize_time(entry.max)
          << std::setw(10) << normalize_time(entry.avg);

        if (_threads == 1)
          s << std::setw(10) << entry.std_dev;

        if (entry.operations)
          {
            s << std::setw(10) << std::setprecision(3) << std::fixed << (entry.operations / normalize_time(entry.avg) * 1e-9);
            if (_threads > 1)
              s << std::setw(14) << std::setprecision(3) << std::fixed << (entry.operations / normalize_time(entry.avg) * 1e-9 / _threads);
          }

        s << std::endl;
      }

      void print(std::ostream& s, bool summary_only = false)
      {
        ios_base_all_saver ios_saver(s);

        if (_statistics_stale)
          update_statistics();

        s << _name << " (" << std::setw(2) << _run << " of " << std::setw(2) << _max_runs << ") runs" << std::endl;

        if (_threads > 1)
          s << "joined statistics from " << _threads << " threads" << std::endl;

        s << std::setw(_max_name_len + 1) << "";

        if (_threads == 1 && !summary_only)
          for (std::size_t i = 0; i < _max_runs; ++i)
            s << std::setw(10) << i;

        s << std::setw(10) << "min"
          << std::setw(10) << "max"
          << std::setw(10) << "avg";

        if (_threads == 1)
          s << std::setw(10) << "std_dev";

        s << std::setw(10) << "GFLOPs";

        if (_threads > 1)
          s << std::setw(14) << "GFLOPs/core";

        s << std::endl;

        for (std::map<std::string,BenchmarkEntry>::const_iterator it = _tasks.begin(), end = _tasks.end();
             it != end;
             ++it)
          print_entry(s,it->first,it->second,summary_only);

        print_entry(s,"total",_run_times,summary_only);
      }

      enum class JoinMode {
        undefined,
        max,
        min,
        avg,
        sum,
      };

      void join(BenchmarkHelper& bh, JoinMode jm)
      {
        if (_join_mode == JoinMode::undefined)
          _join_mode = jm;

        if (_join_mode != jm)
          DUNE_THROW(Exception,"Can't change JoinMode");

        update_statistics();
        bh.update_statistics();
        for (auto& t : _tasks)
          {
            auto& t1 = t.second;
            const auto& t2 = bh._tasks[t.first];
            t1.operations += t2.operations;
            switch (jm)
              {
              case JoinMode::max:
                t1.min = std::max(t1.min,t2.min);
                t1.max = std::max(t1.max,t2.max);
                t1.avg = std::max(t1.avg,t2.avg);
                break;
              case JoinMode::min:
                t1.min = std::min(t1.min,t2.min);
                t1.max = std::min(t1.max,t2.max);
                t1.avg = std::min(t1.avg,t2.avg);
                break;
              case JoinMode::avg:
              case JoinMode::sum:
                t1.min += t2.min;
                t1.max += t2.max;
                t1.avg += t2.avg;
                break;
              default:
                DUNE_THROW(Exception,"Invalid JoinMode");
              }
          }
        ++_threads;
      }

    private:
      const std::string _name;
      TimeSource _time;
      std::size_t _run;
      const std::size_t _max_runs;
      std::map<std::string,BenchmarkEntry> _tasks;
      bool _statistics_stale;
      BenchmarkEntry _run_times;
      std::size_t _max_name_len;
      std::size_t _threads;
      JoinMode _join_mode;

    };



  } // namespace PDELAB
} // namespace Dune

#endif // DUNE_PDELAB_COMMON_BENCHMARKHELPER_HH
