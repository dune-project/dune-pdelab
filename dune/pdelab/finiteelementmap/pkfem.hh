// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_PkFEM_HH
#define DUNE_PDELAB_PkFEM_HH

#include<dune/common/exceptions.hh>

#include<dune/localfunctions/lagrange/pk2d.hh>

#include"finiteelementmap.hh"
#include"pk1dbasis.hh"
#include"pk2dfem.hh"
#include"pk3dfem.hh"

namespace Dune {
  namespace PDELab {

    //! wrap up element from local functions
    //! \ingroup FiniteElementMap

    // dummy template; use only implementations for d=2, d=3
    template<typename GV, typename D, typename R, unsigned int k, unsigned int d>
    class PkLocalFiniteElementMap
    {
    };

    // 1D version
    template<typename GV, typename D, typename R, unsigned int k>
    class PkLocalFiniteElementMap<GV,D,R,k,1> : public Pk1dLocalFiniteElementMap<D,R>
    {
    public:
      PkLocalFiniteElementMap (const GV& gv) : Pk1dLocalFiniteElementMap<D,R>(k)
      {}
    };

    // 2D version
    template<typename GV, typename D, typename R, unsigned int k>
    class PkLocalFiniteElementMap<GV,D,R,k,2> : public Pk2DLocalFiniteElementMap<GV,D,R,k>
    {
    public:
      PkLocalFiniteElementMap (const GV& gv) : Pk2DLocalFiniteElementMap<GV,D,R,k>(gv)
      {}
    };

    // 3D version
    template<typename GV, typename D, typename R, unsigned int k>
    class PkLocalFiniteElementMap<GV,D,R,k,3> : public Pk3DLocalFiniteElementMap<GV,D,R,k>
    {
    public:
      PkLocalFiniteElementMap (const GV& gv) : Pk3DLocalFiniteElementMap<GV,D,R,k>(gv)
      {}
    };

  }
}

#endif
