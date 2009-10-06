#ifndef DUNE_PDELAB_L2ERROR_PROBE_HH
#define DUNE_PDELAB_L2ERROR_PROBE_HH

#include "l2difference.hh"
#include "probe.hh"

//
// L2ErrorProbe
//

class L2ErrorProbe
  : public Dune::PDELab::DummyProbeDefault<L2ErrorProbe>
{
  std::ostream &dat;
  unsigned integrationOrder;

  const double mean_h;
  double error;

public:
  template<typename GV>
  L2ErrorProbe(std::ostream &dat_, unsigned integrationOrder_, const GV &gv)
    : dat(dat_), integrationOrder(integrationOrder_)
    , mean_h(std::pow(1.0/gv.size(0), 1.0/GV::dimension)), error(0)
  { }

  double get_error() const
  { return error; }

  template<typename GF, typename EGF>
  void measureFinalExact(const GF &gf, const EGF &egf, double time = 0) {
    error = l2difference(gf, egf, integrationOrder);
    dat << std::setprecision(8) << mean_h << "\t" << error << std::endl;
  }
};

class L2ErrorLevelProbeFactory
{
  std::ostream &dat;
  unsigned integrationOrder;

public:
  template<typename GV>
  struct Traits {
    typedef L2ErrorProbe Probe;
  };

  L2ErrorLevelProbeFactory(std::ostream &dat_, unsigned integrationOrder_)
    : dat(dat_), integrationOrder(integrationOrder_)
  { }

  template<typename GV>
  Dune::SmartPointer<typename Traits<GV>::Probe>
  getProbe(const GV &gv, unsigned level)
  {
    dat << "# level " << level << std::endl;
    return new typename Traits<GV>::Probe(dat, integrationOrder, gv);
  }
};

class L2ErrorGridProbeFactory
{
  GnuplotGraph graph;
  const unsigned integrationOrder;
  unsigned index;
  std::ostream::pos_type lastpos;
  std::string lastplot;

public:
  template<typename G>
  struct Traits {
    typedef L2ErrorLevelProbeFactory LevelProbeFactory;
  };

  L2ErrorGridProbeFactory(const std::string &fileprefix,
                                       const unsigned integrationOrder_)
    : graph(fileprefix), integrationOrder(integrationOrder_)
    , index(0), lastpos(graph.dat().tellp()), lastplot("")
  {
    graph.addCommand("set terminal postscript eps color solid");
    graph.addCommand("set output '"+fileprefix+".eps'");
    graph.addCommand("");
    graph.addCommand("set key left top reverse Left");
    graph.addCommand("set logscale xy");
    graph.addCommand("set title 'L2 Error'");
    graph.addCommand("set xlabel '<h>'");
    graph.addCommand("set ylabel 'Error'");
    graph.addCommand("");
  }

  ~L2ErrorGridProbeFactory() {
    if(lastpos != graph.dat().tellp())
      graph.addPlot(lastplot);
  }

  template<typename G>
  Dune::SmartPointer<typename Traits<G>::LevelProbeFactory>
  levelProbeFactory(const G &grid, const std::string &tag)
  {
    if(lastpos != graph.dat().tellp()) {
      graph.addPlot(lastplot);
      graph.dat() << "\n\n";
      ++index;
    }
    graph.dat() << "# " << tag << std::endl;
    lastpos = graph.dat().tellp();
    std::ostringstream plotconstruct;
    plotconstruct << "'" << graph.datname() << "'"
                  << " index " << index
                  << " title '" << tag << "'"
                  << " with linespoints pt 1";
    lastplot = plotconstruct.str();
    return new typename Traits<G>::LevelProbeFactory(graph.dat(), integrationOrder);
  }

};

//
// L2ErrorEvolutionProbe
//

class L2ErrorEvolutionProbe
  : public Dune::PDELab::DummyProbeDefault<L2ErrorEvolutionProbe>
{
  std::ostream &dat;
  unsigned integrationOrder;

public:
  L2ErrorEvolutionProbe(std::ostream &dat_, unsigned integrationOrder_)
    : dat(dat_), integrationOrder(integrationOrder_)
  { }

  template<typename GF, typename EGF>
  void measureExact(const GF &gf, const EGF &egf, double time = 0) {
    double error = l2difference(gf, egf, integrationOrder);
    dat << std::setprecision(8) << time << "\t" << error << std::endl;
  }
};

class L2ErrorEvolutionLevelProbeFactory
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
    typedef L2ErrorEvolutionProbe Probe;
  };

  L2ErrorEvolutionLevelProbeFactory(GnuplotGraph &graph_, unsigned integrationOrder_,
                                             unsigned &index_, const std::string &tag_)
    : graph(graph_), integrationOrder(integrationOrder_), datapos(graph.dat().tellp()),
      index(index_), lastplot(""), tag(tag_)
  { }
  ~L2ErrorEvolutionLevelProbeFactory()
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

class L2ErrorEvolutionGridProbeFactory
{
  const unsigned integrationOrder;
  GnuplotGraph graph;
  unsigned index;

public:
  template<typename G>
  struct Traits {
    typedef L2ErrorEvolutionLevelProbeFactory LevelProbeFactory;
  };

  L2ErrorEvolutionGridProbeFactory(const std::string &fileprefix,
                                            const unsigned integrationOrder_)
    : integrationOrder(integrationOrder_), graph(fileprefix), index(0)
  {
    graph.addCommand("set terminal postscript eps color solid");
    graph.addCommand("set output '"+fileprefix+".eps'");
    graph.addCommand("");
    graph.addCommand("set key left top reverse Left");
    graph.addCommand("set logscale y");
    graph.addCommand("set title 'L2 Error Evolution'");
    graph.addCommand("set xlabel 't'");
    graph.addCommand("set ylabel 'Error'");
    graph.addCommand("");
  }

  template<typename G>
  Dune::SmartPointer<typename Traits<G>::LevelProbeFactory>
  levelProbeFactory(const G &grid, const std::string &tag)
  {
    return new typename Traits<G>::LevelProbeFactory(graph, integrationOrder, index, tag);
  }
};

#endif //DUNE_PDELAB_L2ERROR_PROBE_HH
