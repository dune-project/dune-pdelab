// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_VTK_PROBE_HH
#define DUNE_PDELAB_VTK_PROBE_HH

#include <dune/common/smartpointer.hh>

#include <dune/grid/io/file/vtk/vtksequencewriter.hh>

#include "../common/vtkexport.hh"

#include "probe.hh"

//
// VTKProbe
//

template<typename GV>
class ResonatorVTKProbe
  : public Dune::PDELab::DummyProbeDefault<ResonatorVTKProbe<GV> >
{
  Dune::VTKSequenceWriter<GV> writer;
  double timeStretch;

public:
  ResonatorVTKProbe(const GV &gv, const std::string &name, double timeStretch_ = 1)
    : writer(gv, name, ".", "", Dune::VTKOptions::nonconforming)
    , timeStretch(timeStretch_)
  {}

  template<typename GF>
  void measure(const GF &gf, double time = 0) {
    writer.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<GF>
                         (gf,"solution"));
    writer.write(time*timeStretch,Dune::VTKOptions::binaryappended);
    writer.clear();
  }

  template<typename GF, typename EGF>
  void measureExact(const GF &gf, const EGF &egf, double time = 0) {
    writer.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<GF>
                         (gf,"solution"));
    writer.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<EGF>
                         (egf,"reference"));
    writer.write(time*timeStretch,Dune::VTKOptions::binaryappended);
    writer.clear();
  }
};

class ResonatorVTKLevelProbeFactory
{
  const std::string tag;
  const double timeStretch;

public:
  template<typename GV>
  struct Traits {
    typedef ResonatorVTKProbe<GV> Probe;
  };

  ResonatorVTKLevelProbeFactory(const std::string &tag_, const double timeStretch_)
    : tag(tag_), timeStretch(timeStretch_)
  { }

  template<typename GV>
  Dune::SmartPointer<typename Traits<GV>::Probe>
  getProbe(const GV &gv, unsigned level)
  {
    std::ostringstream nameconstruct;
    nameconstruct << tag << ".level" << level;
    return new typename Traits<GV>::Probe(gv, nameconstruct.str(), timeStretch);
  }
};

class ResonatorVTKGridProbeFactory
{
  const std::string fileprefix;
  const double timeStretch;

public:
  template<typename G>
  struct Traits {
    typedef ResonatorVTKLevelProbeFactory LevelProbeFactory;
  };

  ResonatorVTKGridProbeFactory(const std::string &fileprefix_, const double timeStretch_ = 1)
    : fileprefix(fileprefix_), timeStretch(timeStretch_)
  { }

  template<typename G>
  Dune::SmartPointer<typename Traits<G>::LevelProbeFactory>
  levelProbeFactory(const G &grid, const std::string &tag)
  {
    return new typename Traits<G>::LevelProbeFactory(fileprefix + "-" + tag, timeStretch);
  }

};

#endif //DUNE_PDELAB_VTK_PROBE_HH
