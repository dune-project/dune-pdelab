#include "config.h"

#include <dune/pdelab/logging/sink.hh>

namespace Dune::PDELab {

  Sink::~Sink()
  {}

  void NullSink::process(const LogMessage &)
  {}

} // end namespace Dune::PDELab
