#ifndef DUNE_PDELAB_RESONATORSOLUTION_HH
#define DUNE_PDELAB_RESONATORSOLUTION_HH

#include <dune/common/float_cmp.hh>
#include <dune/common/smartpointer.hh>

#include "../common/function.hh"

#include "physicalconstants.hh"

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


#endif //DUNE_PDELAB_RESONATORSOLUTION_HH
