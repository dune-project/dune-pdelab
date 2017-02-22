//-*- tab-width: 4; c-basic-offset: 2; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_FUNCTION_BINDTIME_HH
#define DUNE_PDELAB_FUNCTION_BINDTIME_HH

#include <dune/common/indices.hh>
#include <dune/functions/common/signature.hh>

namespace Dune {
namespace PDELab {

template<typename F, typename Placeholder>
class SetTimeWrapper
{
public:
  SetTimeWrapper(F&& f, Placeholder p)
    : _f(f), _p(p)
  {}
  template<typename Domain>
  auto operator()(const Domain & d) const
  {
    return swap_arguments(d,_t,_p);
  }
  void setTime(double t)
  {
    _t = t;
  }
private:
  double _t;
  F _f;
  Placeholder _p;
  template<typename Domain>
  auto swap_arguments(const Domain& d, double t, const index_constant<1>&) const
  {
    return _f(t,d);
  }
  template<typename Domain>
  auto swap_arguments(const Domain& d, double t, const index_constant<2>&) const
  {
    return _f(d,t);
  }
};

template<typename F, typename Placeholder>
SetTimeWrapper<F,Placeholder>
bindTime(F&& f, Placeholder p)
{
  return SetTimeWrapper<F,Placeholder>(f,p);
}

} // end namespace PDELab
} // end namespace Dune

#endif // DUNE_PDELAB_FUNCTION_BINDTIME_HH
