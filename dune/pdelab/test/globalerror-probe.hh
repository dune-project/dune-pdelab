#ifndef DUNE_PDELAB_GLOBALERROR_PROBE_HH
#define DUNE_PDELAB_GLOBALERROR_PROBE_HH

#include <dune/common/smartpointer.hh>

#include "l2difference.hh"
#include "probe.hh"

//
// GlobalErrorProbe
//

class GlobalErrorProbe
  : public Dune::PDELab::DummyProbeDefault<GlobalErrorProbe>
{
  unsigned integrationOrder;
  double mean_h;
  unsigned nsamples;
  double sum;

  std::ostream *dat;

  std::ostream *evodat;

public:
  template<typename GV>
  GlobalErrorProbe(const GV &gv, unsigned integrationOrder_, std::ostream *dat_, std::ostream *evodat_)
    : integrationOrder(integrationOrder_), mean_h(std::pow(1.0/gv.size(0), 1.0/GV::dimension)), nsamples(0), sum(0)
    , dat(dat_), evodat(evodat_)
  { }

  ~GlobalErrorProbe() {
    if(dat)
      *dat << std::setprecision(8) << mean_h << "\t" << get_error() << std::endl;
  }

  double get_error() const
  { return std::sqrt(sum/nsamples); }

  template<typename GF, typename EGF>
  void measureExact(const GF &gf, const EGF &egf, double time = 0) {
    sum += l2difference2(gf, egf, integrationOrder);
    ++nsamples;

    if(evodat)
      *evodat << std::setprecision(8) << time << "\t" << get_error() << std::endl;
  }
};

class GlobalErrorLevelProbeFactory
{
  unsigned integrationOrder;

  std::ostream *dat;

  GnuplotGraph *evograph;
  std::ostream::pos_type evodatapos;
  unsigned &evoindex;
  std::string evolastplot;
  const std::string tag;

public:
  template<typename GV>
  struct Traits {
    typedef GlobalErrorProbe Probe;
  };

  GlobalErrorLevelProbeFactory(unsigned integrationOrder_,
                               std::ostream *dat_,
                               GnuplotGraph *evograph_, unsigned &evoindex_,
                               const std::string &tag_)
    : integrationOrder(integrationOrder_)
    , dat(dat_)
    , evograph(evograph_), evodatapos(0), evoindex(evoindex_), evolastplot("")
    , tag(tag_)
  { 
    if(evograph)
      evodatapos = evograph->dat().tellp();
  }

  ~GlobalErrorLevelProbeFactory() {
    if(evograph)
      if(evodatapos != evograph->dat().tellp()) {
        evograph->dat() << "\n\n";
        evograph->addPlot(evolastplot);
        ++evoindex;
      }
  }

  template<typename GV>
  Dune::SmartPointer<typename Traits<GV>::Probe>
  getProbe(const GV &gv, unsigned level)
  {
    if(dat)
      *dat << "# level " << level << std::endl;

    if(evograph) {
      if(evodatapos != evograph->dat().tellp()) {
        evograph->dat() << "\n\n";
        evograph->addPlot(evolastplot);
        ++evoindex;
      }
      evograph->dat() << "# LEVEL" << level << std::endl;
      evodatapos = evograph->dat().tellp();
      std::ostringstream plotconstruct;
      plotconstruct << "'" << evograph->datname() << "'"
                    << " index " << evoindex
                    << " title '" << tag << " level " << level << "'"
                    << " with lines";
      evolastplot = plotconstruct.str();
    }
      
    return new typename Traits<GV>::Probe(gv, integrationOrder,
                                          dat,
                                          evograph ? &evograph->dat() : 0);
  }
};

class GlobalErrorGridProbeFactory
{
  const unsigned integrationOrder;

  Dune::SmartPointer<GnuplotGraph> graph;
  unsigned index;
  std::ostream::pos_type lastpos;
  std::string lastplot;

  bool do_evolution;
  Dune::SmartPointer<GnuplotGraph> evograph;
  unsigned evoindex;

public:
  template<typename G>
  struct Traits {
    typedef GlobalErrorLevelProbeFactory LevelProbeFactory;
  };

  GlobalErrorGridProbeFactory(const unsigned integrationOrder_,
                              const std::string &convergencePrefix,
                              const std::string &evolutionPrefix = "")
    : integrationOrder(integrationOrder_)
    , graph(0), index(0), lastpos(0), lastplot("")
    , evograph(0), evoindex(0)
  {
    if(convergencePrefix != "") {
      graph = new GnuplotGraph(convergencePrefix);
      lastpos = graph->dat().tellp();

      graph->addCommand("set terminal postscript eps color solid");
      graph->addCommand("set output '"+convergencePrefix+".eps'");
      graph->addCommand("");
      graph->addCommand("set key left top reverse Left");
      graph->addCommand("set logscale xy");
      graph->addCommand("set title 'GlobalError'");
      graph->addCommand("set xlabel '<h>'");
      graph->addCommand("set ylabel 'Error'");
      graph->addCommand("");
    }

    if(evolutionPrefix != "") {
      evograph = new GnuplotGraph(evolutionPrefix);

      evograph->addCommand("set terminal postscript eps color solid");
      evograph->addCommand("set output '"+evolutionPrefix+".eps'");
      evograph->addCommand("");
      evograph->addCommand("set key left top reverse Left");
      evograph->addCommand("set logscale y");
      evograph->addCommand("set title 'Global Error Evolution'");
      evograph->addCommand("set xlabel 't'");
      evograph->addCommand("set ylabel 'Error'");
      evograph->addCommand("");
    }
  }

  ~GlobalErrorGridProbeFactory() {
    if(&*graph)
      if(lastpos != graph->dat().tellp())
        graph->addPlot(lastplot);
  }

  template<typename G>
  Dune::SmartPointer<typename Traits<G>::LevelProbeFactory>
  levelProbeFactory(const G &grid, const std::string &tag)
  {
    if(&*graph) {
      if(lastpos != graph->dat().tellp()) {
        graph->addPlot(lastplot);
        graph->dat() << "\n\n";
        ++index;
      }
      graph->dat() << "# " << tag << std::endl;
      lastpos = graph->dat().tellp();
      std::ostringstream plotconstruct;
      plotconstruct << "'" << graph->datname() << "'"
                    << " index " << index
                    << " title '" << tag << "'"
                    << " with linespoints pt 1";
      lastplot = plotconstruct.str();
    }
    return new typename Traits<G>::LevelProbeFactory(integrationOrder,
                                                     &*graph ? &graph->dat() : 0,
                                                     &*evograph, evoindex, tag);
  }
};


#endif //DUNE_PDELAB_GLOBALERROR_PROBE_HH
