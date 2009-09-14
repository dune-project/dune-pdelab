// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_RT0QFEM_HH
#define DUNE_PDELAB_RT0QFEM_HH

#include"rt0q2dfem.hh"
#include"rt0q3dfem.hh"

namespace Dune {
  namespace PDELab {

    template<typename GV, typename D, typename R, int d>
    class RT0QLocalFiniteElementMap;

    template<typename GV, typename D, typename R>
    class RT0QLocalFiniteElementMap<GV,D,R,2> 
      : public RT0Q2DLocalFiniteElementMap<GV,D,R>
    {
    public:
      RT0QLocalFiniteElementMap (const GV& gv_) 
        : RT0Q2DLocalFiniteElementMap<GV,D,R>(gv_)
      {
      }
    };

    template<typename GV, typename D, typename R>
    class RT0QLocalFiniteElementMap<GV,D,R,3> 
      : public RT0Q3DLocalFiniteElementMap<GV,D,R>
    {
    public:
      RT0QLocalFiniteElementMap (const GV& gv_) 
        : RT0Q3DLocalFiniteElementMap<GV,D,R>(gv_)
      {
      }
    };
  }
}

#endif
