#ifndef DUNE_PDELAB_ELECTRICENERGY_PROBE_HH
#define DUNE_PDELAB_ELECTRICENERGY_PROBE_HH

#include "l2norm.hh"
#include "probe.hh"

//
// Electric Energy Probe
//

class ElectricEnergyProbe
  : public Dune::PDELab::DummyProbeDefault<ElectricEnergyProbe>
{
  std::ostream &dat;
  unsigned integrationOrder;

public:
  ElectricEnergyProbe(std::ostream &dat_, unsigned integrationOrder_)
    : dat(dat_), integrationOrder(integrationOrder_)
  { }

  template<typename GF>
  void measure(const GF &gf, double time = 0) {
    dat << std::setprecision(8) << time << "\t" << l2norm(gf, integrationOrder) << std::endl;
  }
};

class ElectricEnergyLevelProbeFactory
{
  GnuplotGraph &graph;
  unsigned integrationOrder;
  std::ostream::pos_type datapos;
  unsigned &index;
  std::string lastplot;
  const std::string tag;

public:
  template<typename GV>
  struct Traits {
    typedef ElectricEnergyProbe Probe;
  };

  ElectricEnergyLevelProbeFactory(GnuplotGraph &graph_, unsigned integrationOrder_,
                             unsigned &index_, const std::string &tag_)
    : graph(graph_), integrationOrder(integrationOrder_), datapos(graph.dat().tellp()),
      index(index_), lastplot(""), tag(tag_)
  { }
  ~ElectricEnergyLevelProbeFactory()
  {
    if(datapos != graph.dat().tellp()) {
      graph.dat() << "\n\n";
      graph.addPlot(lastplot);
      ++index;
    }
  }

  template<typename GV>
  Dune::SmartPointer<typename Traits<GV>::Probe>
  getProbe(const GV &gv, unsigned level)
  {
    if(datapos != graph.dat().tellp()) {
      graph.dat() << "\n\n";
      graph.addPlot(lastplot);
      ++index;
    }
    graph.dat() << "# LEVEL" << level << std::endl;
    datapos = graph.dat().tellp();
    std::ostringstream plotconstruct;
    plotconstruct << "'" << graph.datname() << "'"
                  << " index " << index
                  << " title '" << tag << " level " << level << "'"
                  << " with lines";
    lastplot = plotconstruct.str();
    return new typename Traits<GV>::Probe(graph.dat(), integrationOrder);
  }
};

class ElectricEnergyGridProbeFactory
{
  const unsigned integrationOrder;
  GnuplotGraph graph;
  unsigned index;

public:
  template<typename G>
  struct Traits {
    typedef ElectricEnergyLevelProbeFactory LevelProbeFactory;
  };

  ElectricEnergyGridProbeFactory(const std::string &fileprefix,
                                 const unsigned integrationOrder_)
    : integrationOrder(integrationOrder_), graph(fileprefix), index(0)
  {
    graph.addCommand("set terminal postscript eps color solid");
    graph.addCommand("set output '"+fileprefix+".eps'");
    graph.addCommand("");
    graph.addCommand("set key left top reverse Left");
    graph.addCommand("set title 'Electric Energy Evolution'");
    graph.addCommand("set xlabel 't'");
    graph.addCommand("set ylabel 'Energy'");
    graph.addCommand("");
  }

  template<typename G>
  Dune::SmartPointer<typename Traits<G>::LevelProbeFactory>
  levelProbeFactory(const G &grid, const std::string &tag)
  {
    return new typename Traits<G>::LevelProbeFactory(graph, integrationOrder, index, tag);
  }
};

#endif //DUNE_PDELAB_ELECTRICENERGY_PROBE_HH
