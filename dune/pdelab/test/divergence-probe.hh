// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_DIVERGENCE_PROBE_HH
#define DUNE_PDELAB_DIVERGENCE_PROBE_HH

#include <cmath>

#include <dune/common/smartpointer.hh>

#include <dune/grid/common/quadraturerules.hh>

#include "probe.hh"

template<typename T>
inline T sqr(const T& v) {
  return v*v;
}

// Calculate the squared L2 norm of the divergence of a function
// it is assumed the function does not have divergence on the elements itself
// (only between the elements) and that its value is 0 outside the domain.
template<typename U>
typename U::Traits::RangeFieldType
divergencel2norm2 (const U& u, int qorder=2)
{
  // constants and types
  typedef typename U::Traits::GridViewType GV;
  typedef typename GV::IndexSet IS;
  static const int dimDomain = U::Traits::dimDomain;
  typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
  typedef typename GV::IntersectionIterator IntersectionIterator;
  typedef typename U::Traits::DomainFieldType ct;
  typedef Dune::QuadratureRule<ct,dimDomain-1> QR;
  typedef Dune::QuadratureRules<ct,dimDomain-1> QRs;

  const GV& gv = u.getGridView();
  const IS& is = gv.indexSet();
  typename U::Traits::RangeType u_val;
  typename U::Traits::RangeType u_val_other;
  
  // loop over grid view
  typename U::Traits::RangeFieldType sum = 0.0;
  for (ElementIterator eit = gv.template begin<0>();
       eit != gv.template end<0>(); ++eit)
    for (IntersectionIterator iit = gv.ibegin(*eit);
         iit != gv.iend(*eit); ++iit)
    {
      if(iit->neighbor() && is.index(*eit) > is.index(*iit->outside()))
        continue;

      Dune::GeometryType gt = iit->type();
      const QR& rule = QRs::rule(gt,qorder);

      for (typename QR::const_iterator qit=rule.begin(); qit!=rule.end();
           ++qit)
      {
        // evaluate the given grid functions at integration point
        u.evaluate(*eit,
                   iit->geometryInInside().global(qit->position()),
                   u_val);
        if(iit->neighbor()) {
          u.evaluate(*iit->outside(),
                     iit->geometryInOutside().global(qit->position()),
                     u_val_other);
          u_val -= u_val_other;
        }

        // accumulate error
        sum += sqr(u_val*iit->unitOuterNormal(qit->position()))*qit->weight()*
          iit->geometry().integrationElement(qit->position());
      }
    }
  return sum;
}

//
// Divergence Probe
//

class DivergenceProbe
  : public Dune::PDELab::DummyProbe
{
  std::ostream &dat;
  unsigned integrationOrder;

public:
  DivergenceProbe(std::ostream &dat_, unsigned integrationOrder_)
    : dat(dat_), integrationOrder(integrationOrder_)
  { }

  template<typename GF>
  void measure(const GF &gf, double time = 0) {
    dat << std::setprecision(8) << time << "\t"
        << std::sqrt(divergencel2norm2(gf, integrationOrder)) << std::endl;
  }
};

class DivergenceLevelProbeFactory
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
    typedef DivergenceProbe Probe;
  };

  DivergenceLevelProbeFactory(GnuplotGraph &graph_, unsigned integrationOrder_,
                             unsigned &index_, const std::string &tag_)
    : graph(graph_), integrationOrder(integrationOrder_),
      datapos(graph.dat().tellp()),
      index(index_), lastplot(""), tag(tag_)
  { }
  ~DivergenceLevelProbeFactory()
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

class DivergenceGridProbeFactory
{
  const unsigned integrationOrder;
  GnuplotGraph graph;
  unsigned index;

public:
  template<typename G>
  struct Traits {
    typedef DivergenceLevelProbeFactory LevelProbeFactory;
  };

  DivergenceGridProbeFactory(const std::string &fileprefix,
                                 const unsigned integrationOrder_)
    : integrationOrder(integrationOrder_), graph(fileprefix), index(0)
  {
    graph.addCommand("set terminal postscript eps color solid");
    graph.addCommand("set output '"+fileprefix+".eps'");
    graph.addCommand("");
    graph.addCommand("set key left top reverse Left");
    graph.addCommand("set title 'Divergence Evolution'");
    graph.addCommand("set xlabel 't'");
    graph.addCommand("set ylabel 'Divergence'");
    graph.addCommand("");
  }

  template<typename G>
  Dune::SmartPointer<typename Traits<G>::LevelProbeFactory>
  levelProbeFactory(const G &grid, const std::string &tag)
  {
    return new typename Traits<G>::LevelProbeFactory(graph, integrationOrder,
                                                     index, tag);
  }
};

#endif //DUNE_PDELAB_DIVERGENCE_PROBE_HH
