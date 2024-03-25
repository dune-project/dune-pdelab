#ifndef DUNE_PDELAB_COMMON_TRACE_HH
#define DUNE_PDELAB_COMMON_TRACE_HH

// This file allows tracing with prefetto when available
// See https://perfetto.dev/docs/instrumentation/tracing-sdk for more info

#if HAVE_PERFETTO
#include <perfetto.h>
#include <filesystem>
#include <thread>
#include <fstream>
#include <sstream>
#include <iostream>

#define PERFETTO_CATEGORIES

PERFETTO_DEFINE_CATEGORIES(
    perfetto::Category("dune")
        .SetDescription("Events from DUNE")
);

namespace Dune::PDELab::inline Experimental {

struct [[nodiscard]] TracingSession {

  TracingSession(std::filesystem::path filename)
  : _filename{filename} {

    perfetto::TracingInitArgs args;
    args.backends |= perfetto::kInProcessBackend;
    perfetto::Tracing::Initialize(args);
    perfetto::TrackEvent::Register();

    perfetto::TraceConfig cfg;
    cfg.add_buffers()->set_size_kb(1024*100);
    auto* ds_cfg = cfg.add_data_sources()->mutable_config();
    ds_cfg->set_name("track_event");
    _trace = perfetto::Tracing::NewTrace();
    _trace->Setup(cfg);
    _trace->StartBlocking();

    // Give a custom name for the traced process.
    perfetto::ProcessTrack process_track = perfetto::ProcessTrack::Current();
    perfetto::protos::gen::TrackDescriptor desc = process_track.Serialize();
    desc.mutable_process()->set_process_name("Main");
    perfetto::TrackEvent::SetTrackDescriptor(process_track, desc);

    TRACE_EVENT_BEGIN("dune", "Main");
  }

  ~TracingSession() {
    TRACE_EVENT_END("dune");

    // Make sure the last event is closed for this example.
    perfetto::TrackEvent::Flush();

    // Stop tracing and read the trace data.
    _trace->StopBlocking();
    std::vector<char> trace_data(_trace->ReadTraceBlocking());

    // Write the result into a file.
    // Note: To save memory with longer traces, you can tell Perfetto to write
    // directly into a file by passing a file descriptor into Setup() above.
    std::ostringstream ss; ss << "." << std::this_thread::get_id();
    _filename.replace_extension("trace" + ss.str());
    std::ofstream output;
    output.open(_filename.c_str(), std::ios::out | std::ios::binary);
    output.write(&trace_data[0], std::streamsize(trace_data.size()));
    output.close();
    std::cout << "Trace written in `"<< _filename.c_str() << "`." << std::endl;
  }

private:
  std::unique_ptr<perfetto::TracingSession> _trace;
  std::filesystem::path _filename;
};

} // namespace Dune::PDELab::inline Experimental

#else

#ifndef TRACE_EVENT
#define TRACE_EVENT(...) {};
#endif // TRACE_EVENT

#ifndef TRACE_EVENT_BEGIN
#define TRACE_EVENT_BEGIN(...) {};
#endif // TRACE_EVENT_BEGIN

#ifndef TRACE_EVENT_END
#define TRACE_EVENT_END(...) {};
#endif // TRACE_EVENT_END


#ifndef TRACE_COUNTER
#define TRACE_COUNTER(...) {};
#endif // TRACE_COUNTER

namespace perfetto::TrackEvent {
  static constexpr inline auto GetTraceTimeNs() { return int{-1}; }
}

#endif // HAVE_PERFETTO

#endif // DUNE_PDELAB_COMMON_TRACE_HH
