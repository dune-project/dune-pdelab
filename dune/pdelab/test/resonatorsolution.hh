#ifndef DUNE_PDELAB_RESONATORSOLUTION_HH
#define DUNE_PDELAB_RESONATORSOLUTION_HH

#include <dune/common/float_cmp.hh>
#include <dune/common/smartpointer.hh>

#include <dune/grid/io/file/vtk/vtksequencewriter.hh>

#include "../common/function.hh"
#include "../common/vtkexport.hh"

#include "l2difference.hh"
#include "physicalconstants.hh"
#include "probe.hh"

//======================================================================
// Analytic solution
//======================================================================

template<typename GV, typename RF>
class ResonatorSolution
  : public Dune::PDELab::AnalyticGridFunctionBase<
      Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,GV::dimension>,
      ResonatorSolution<GV,RF>
    >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,GV::dimension> Traits;

private:
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,ResonatorSolution<GV,RF> > BaseT;

public:

  ResonatorSolution (const GV& gv,
                     const typename Traits::DomainType &k_,
                     const typename Traits::RangeType &amp_,
                     const typename Traits::DomainType &origin_ = typename Traits::DomainType(0))
    : BaseT(gv)
    , k(k_)
    , amp(amp_)
    , origin(origin_)
  {}

  inline void
  evaluateGlobal (const typename Traits::DomainType& x, 
                  typename Traits::RangeType& y) const
  {
    y = amp;
    for(unsigned i = 0; i < Traits::dimDomain; ++i)
      for(unsigned j = 0; j < Traits::dimDomain; ++j)
        if(i == j)
          y[i] *= std::cos(k[j]*x[j]);
        else
          y[i] *= std::sin(k[j]*x[j]);
  }

private:
  typename Traits::DomainType k;
  typename Traits::RangeType amp;
  typename Traits::DomainType origin;
  
};


template<typename GV, typename RF>
class ResonatorSolutionFactory
{
public:
  typedef ResonatorSolution<GV, RF> Function;

private:
  typedef typename Function::Traits Traits;

  // return a unit vector along cartesian axis dir
  static typename Traits::RangeType
  makeUnitVector(unsigned dir = 0)
  {
    typename Traits::RangeType u(0);
    u[dir] = 1;
    return u;
  }

  // return a modes vector describing the lowest mode when the field points in direction dir
  static Dune::FieldVector<unsigned, Traits::dimDomain>
  makeModesVector(unsigned dir = 0)
  {
    Dune::FieldVector<unsigned, Traits::dimDomain> m(1);
    m[dir] = 0;
    return m;
  }

public:
  ResonatorSolutionFactory
  (const Dune::FieldVector<unsigned, Traits::dimDomain> &modes = makeModesVector(),
   typename Traits::DomainFieldType time0_ = 0,
   const typename Traits::DomainType &size = typename Traits::DomainType(1),
   const typename Traits::RangeType &amp_ = makeUnitVector(),
   const typename Traits::DomainType &origin_ = typename Traits::DomainType(0),
   bool forceAmp = true)
    : time0(time0_)
    , amp(amp_)
    , origin(origin_)
  {
    unsigned zeros = 0; // count zero entries
    for(unsigned i = 0; i < Traits::dimDomain; ++i)
      if(modes[i] == 0)
        ++zeros;
    if(zeros > 1)
      DUNE_THROW(Dune::Exception, "Invalid mode selected: at most one mode-number may be zero.  Mode is: " << modes);

    for(unsigned i = 0; i < Traits::dimDomain; ++i)
      k[i] = pi*modes[i]/size[i];
    typename Traits::DomainFieldType k_norm = k.two_norm();
    typename Traits::DomainType unit_k(k);
    unit_k /= k_norm;

    typename Traits::RangeFieldType amp_longitudinal = amp*unit_k;
    if(!forceAmp &&
       Dune::FloatCmp::eq<typename Traits::RangeFieldType, Dune::FloatCmp::absolute>(amp_longitudinal*amp.two_norm(), 0))
      DUNE_THROW(Dune::Exception, "Amplitude (" << amp << ") must be perpendicular to wave number k (" << k << ")");

    // residual longitudinal part is below our threshold or forceAmp is set
    amp.axpy(-amp_longitudinal, unit_k);
  }

  Dune::SmartPointer<Function> function(const GV &gv, typename Traits::DomainFieldType time) const
  {
    typename Traits::RangeFieldType freq = c0*k.two_norm();
    typename Traits::RangeType amp_with_time(amp);
    amp_with_time *= std::sin(freq*(time-time0));
    return new Function(gv, k, amp_with_time, origin);
  }

private:
  typename Traits::DomainFieldType time0;
  typename Traits::DomainType k;
  typename Traits::RangeType amp;
  typename Traits::DomainType origin;
};

//
// VTKProbe
//

template<typename GV, typename RF>
class ResonatorVTKProbe
  : public Dune::PDELab::DummyProbeDefault<ResonatorVTKProbe<GV, RF> >
{
  Dune::VTKSequenceWriter<GV> writer;
  const ResonatorSolutionFactory<GV, RF> rf;
  double timeStretch;

public:
  ResonatorVTKProbe(const GV &gv, const std::string &name, double timeStretch_ = 1)
    : writer(gv, name, ".", "", Dune::VTKOptions::nonconforming)
    , rf(), timeStretch(timeStretch_)
  {}

  template<typename GF>
  void measure(const GF &gf, double time = 0) {
    typedef typename ResonatorSolutionFactory<GV, RF>::Function Ref;

    Dune::SmartPointer<Ref> reference(rf.function(gf.getGridView(), time));
    writer.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<GF>(gf,"solution"));
    writer.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<Ref>(*reference,"reference"));
    writer.write(time*timeStretch,Dune::VTKOptions::binaryappended);
    writer.clear();
  }

};

template<typename RF>
class ResonatorVTKLevelProbeFactory
{
  const std::string tag;
  const double timeStretch;

public:
  template<typename GV>
  struct Traits {
    typedef ResonatorVTKProbe<GV, RF> Probe;
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

template<typename RF>
class ResonatorVTKGridProbeFactory
{
  const std::string fileprefix;
  const double timeStretch;

public:
  template<typename G>
  struct Traits {
    typedef ResonatorVTKLevelProbeFactory<RF> LevelProbeFactory;
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

//
// GlobalErrorProbe
//

class ResonatorGlobalErrorProbe
  : public Dune::PDELab::DummyProbeDefault<ResonatorGlobalErrorProbe>
{
  unsigned integrationOrder;
  double mean_h;
  unsigned nsamples;
  double sum;

  std::ostream *dat;

  std::ostream *evodat;

public:
  template<typename GV>
  ResonatorGlobalErrorProbe(const GV &gv, unsigned integrationOrder_, std::ostream *dat_, std::ostream *evodat_)
    : integrationOrder(integrationOrder_), mean_h(std::pow(1.0/gv.size(0), 1.0/GV::dimension)), nsamples(0), sum(0)
    , dat(dat_), evodat(evodat_)
  { }

  ~ResonatorGlobalErrorProbe() {
    if(dat)
      *dat << std::setprecision(8) << mean_h << "\t" << get_error() << std::endl;
  }

  double get_error() const
  { return std::sqrt(sum/nsamples); }

  template<typename GF>
  void measure(const GF &gf, double time = 0) {
    typedef ResonatorSolutionFactory<
      typename GF::Traits::GridViewType,
      typename GF::Traits::RangeFieldType> RF;
    sum += l2difference2(gf, *RF().function(gf.getGridView(), time),
                         integrationOrder);
    ++nsamples;

    if(evodat)
      *evodat << std::setprecision(8) << time << "\t" << get_error() << std::endl;
  }
};

class ResonatorGlobalErrorLevelProbeFactory
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
    typedef ResonatorGlobalErrorProbe Probe;
  };

  ResonatorGlobalErrorLevelProbeFactory(unsigned integrationOrder_,
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

  ~ResonatorGlobalErrorLevelProbeFactory() {
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

class ResonatorGlobalErrorGridProbeFactory
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
    typedef ResonatorGlobalErrorLevelProbeFactory LevelProbeFactory;
  };

  ResonatorGlobalErrorGridProbeFactory(const unsigned integrationOrder_,
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

  ~ResonatorGlobalErrorGridProbeFactory() {
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

//
// L2ErrorProbe
//

class ResonatorL2ErrorProbe
  : public Dune::PDELab::DummyProbeDefault<ResonatorL2ErrorProbe>
{
  std::ostream &dat;
  unsigned integrationOrder;

  const double mean_h;
  double error;

public:
  template<typename GV>
  ResonatorL2ErrorProbe(std::ostream &dat_, unsigned integrationOrder_, const GV &gv)
    : dat(dat_), integrationOrder(integrationOrder_)
    , mean_h(std::pow(1.0/gv.size(0), 1.0/GV::dimension)), error(0)
  { }

  double get_error() const
  { return error; }

  template<typename GF>
  void measureFinal(const GF &gf, double time = 0) {
    typedef ResonatorSolutionFactory<
      typename GF::Traits::GridViewType,
      typename GF::Traits::RangeFieldType> RSF;

    error = l2difference(gf, *RSF().function(gf.getGridView(), time),
                         integrationOrder);
    dat << std::setprecision(8) << mean_h << "\t" << error << std::endl;
  }
};

class ResonatorL2ErrorLevelProbeFactory
{
  std::ostream &dat;
  unsigned integrationOrder;

public:
  template<typename GV>
  struct Traits {
    typedef ResonatorL2ErrorProbe Probe;
  };

  ResonatorL2ErrorLevelProbeFactory(std::ostream &dat_, unsigned integrationOrder_)
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

class ResonatorL2ErrorGridProbeFactory
{
  GnuplotGraph graph;
  const unsigned integrationOrder;
  unsigned index;
  std::ostream::pos_type lastpos;
  std::string lastplot;

public:
  template<typename G>
  struct Traits {
    typedef ResonatorL2ErrorLevelProbeFactory LevelProbeFactory;
  };

  ResonatorL2ErrorGridProbeFactory(const std::string &fileprefix,
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

  ~ResonatorL2ErrorGridProbeFactory() {
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

class ResonatorL2ErrorEvolutionProbe
  : public Dune::PDELab::DummyProbeDefault<ResonatorL2ErrorEvolutionProbe>
{
  std::ostream &dat;
  unsigned integrationOrder;

public:
  ResonatorL2ErrorEvolutionProbe(std::ostream &dat_, unsigned integrationOrder_)
    : dat(dat_), integrationOrder(integrationOrder_)
  { }

  template<typename GF>
  void measure(const GF &gf, double time = 0) {
    typedef ResonatorSolutionFactory<
      typename GF::Traits::GridViewType,
      typename GF::Traits::RangeFieldType> RSF;

    double error = l2difference(gf, *RSF().function(gf.getGridView(), time),
                                integrationOrder);
    dat << std::setprecision(8) << time << "\t" << error << std::endl;
  }
};

class ResonatorL2ErrorEvolutionLevelProbeFactory
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
    typedef ResonatorL2ErrorEvolutionProbe Probe;
  };

  ResonatorL2ErrorEvolutionLevelProbeFactory(GnuplotGraph &graph_, unsigned integrationOrder_,
                                             unsigned &index_, const std::string &tag_)
    : graph(graph_), integrationOrder(integrationOrder_), datapos(graph.dat().tellp()),
      index(index_), lastplot(""), tag(tag_)
  { }
  ~ResonatorL2ErrorEvolutionLevelProbeFactory()
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

class ResonatorL2ErrorEvolutionGridProbeFactory
{
  const unsigned integrationOrder;
  GnuplotGraph graph;
  unsigned index;

public:
  template<typename G>
  struct Traits {
    typedef ResonatorL2ErrorEvolutionLevelProbeFactory LevelProbeFactory;
  };

  ResonatorL2ErrorEvolutionGridProbeFactory(const std::string &fileprefix,
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

#endif //DUNE_PDELAB_RESONATORSOLUTION_HH
