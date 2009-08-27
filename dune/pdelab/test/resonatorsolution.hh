#ifndef DUNE_PDELAB_RESONATORSOLUTION_HH
#define DUNE_PDELAB_RESONATORSOLUTION_HH

#include <dune/common/smartpointer.hh>

#include <dune/grid/io/file/vtk/vtksequencewriter.hh>

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
  : public Dune::PDELab::DummyProbe
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

template<typename GV, typename RF>
class ResonatorGlobalErrorProbe
  : public Dune::PDELab::DummyProbe
{
  const ResonatorSolutionFactory<GV, RF> rf;
  std::ostream &dat;
  unsigned integrationOrder;

  unsigned nsamples;
  double sum;
  double mean_h;

public:
  ResonatorGlobalErrorProbe(std::ostream &dat_, unsigned integrationOrder_, const GV &gv)
    : rf(), dat(dat_), integrationOrder(integrationOrder_), nsamples(0), sum(0)
    , mean_h(std::pow(1.0/gv.size(0), 1.0/GV::dimension))
  { }

  ~ResonatorGlobalErrorProbe() {
    dat << std::setprecision(8) << mean_h << "\t" << std::sqrt(sum/nsamples) << std::endl;
  }

  template<typename GF>
  void measure(const GF &gf, double time = 0) {
    sum += l2difference2(gf.getGridView(),
                         gf, *rf.function(gf.getGridView(), time),
                         integrationOrder);
    ++nsamples;
  }
};

template<typename RF>
class ResonatorGlobalErrorLevelProbeFactory
{
  std::ostream &dat;
  unsigned integrationOrder;

public:
  template<typename GV>
  struct Traits {
    typedef ResonatorGlobalErrorProbe<GV, RF> Probe;
  };

  ResonatorGlobalErrorLevelProbeFactory(std::ostream &dat_, unsigned integrationOrder_)
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

template<typename RF>
class ResonatorGlobalErrorGridProbeFactory
{
  GnuplotGraph graph;
  const unsigned integrationOrder;
  unsigned index;
  std::ostream::pos_type lastpos;
  std::string lastplot;

public:
  template<typename G>
  struct Traits {
    typedef ResonatorGlobalErrorLevelProbeFactory<RF> LevelProbeFactory;
  };

  ResonatorGlobalErrorGridProbeFactory(const std::string &fileprefix,
                                       const unsigned integrationOrder_)
    : graph(fileprefix), integrationOrder(integrationOrder_)
    , index(0), lastpos(graph.dat().tellp()), lastplot("")
  {
    graph.addCommand("set terminal postscript eps color solid");
    graph.addCommand("set output '"+fileprefix+".eps'");
    graph.addCommand("");
    graph.addCommand("set key left top reverse Left");
    graph.addCommand("set logscale xy");
    graph.addCommand("set title 'GlobalError'");
    graph.addCommand("set xlabel '<h>'");
    graph.addCommand("set ylabel 'Error'");
    graph.addCommand("");
  }

  ~ResonatorGlobalErrorGridProbeFactory() {
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

#endif //DUNE_PDELAB_RESONATORSOLUTION_HH
