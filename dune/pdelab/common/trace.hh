#ifndef DUNE_PDELAB_COMMON_TRACE_HH
#define DUNE_PDELAB_COMMON_TRACE_HH

// This file allows tracing with prefetto when available
// See https://perfetto.dev/docs/instrumentation/tracing-sdk for more info

#if HAVE_PERFETTO
#include <perfetto.h>

#define PERFETTO_CATEGORIES

PERFETTO_DEFINE_CATEGORIES(
    perfetto::Category("dune")
        .SetDescription("Events from DUNE")
);

#else

#ifndef TRACE_EVENT
#define TRACE_EVENT(...) {};
#endif // TRACE_EVENT

#ifndef TRACE_COUNTER
#define TRACE_COUNTER(...) {};
#endif // TRACE_COUNTER

namespace perfetto::TrackEvent {
  static constexpr inline auto GetTraceTimeNs() { return int{-1}; }
}

#endif // HAVE_PERFETTO

#endif // DUNE_PDELAB_COMMON_TRACE_HH
