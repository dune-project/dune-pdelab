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


    template<unsigned int d>
    struct DUNE_DEPRECATED_MSG("The fifth template parameter for PkLocalFiniteElementMap (the dimension) is deprecated") warn_deprecated_pk_interface
    {
      warn_deprecated_pk_interface() DUNE_DEPRECATED_MSG("The fifth template parameter for PkLocalFiniteElementMap (the dimension) is deprecated")
        {}
    };

    template<>
    struct warn_deprecated_pk_interface<0>
    {
      warn_deprecated_pk_interface()
      {}
    };

    //! wrap up element from local functions
    //! \ingroup FiniteElementMap


    template<typename GV, typename D, typename R, unsigned int k, unsigned int d>
    class PkFEMDimensionDispatch
    {};

    // 1D version
    template<typename GV, typename D, typename R, unsigned int k>
    class PkFEMDimensionDispatch<GV,D,R,k,1> : public Pk1dLocalFiniteElementMap<D,R>
    {
    public:
      PkFEMDimensionDispatch (const GV& gv) : Pk1dLocalFiniteElementMap<D,R>(k)
      {}
    };

    // 2D version
    template<typename GV, typename D, typename R, unsigned int k>
    class PkFEMDimensionDispatch<GV,D,R,k,2> : public Pk2DLocalFiniteElementMap<GV,D,R,k>
    {
    public:
      PkFEMDimensionDispatch (const GV& gv) : Pk2DLocalFiniteElementMap<GV,D,R,k>(gv)
      {}
    };

    // 3D version
    template<typename GV, typename D, typename R, unsigned int k>
    class PkFEMDimensionDispatch<GV,D,R,k,3> : public Pk3DLocalFiniteElementMap<GV,D,R,k>
    {
    public:
      PkFEMDimensionDispatch (const GV& gv) : Pk3DLocalFiniteElementMap<GV,D,R,k>(gv)
      {}
    };



    // dummy template; use only implementations for d=2, d=3
    template<typename GV, typename D, typename R, unsigned int k, unsigned int d = 0>
    class PkLocalFiniteElementMap
      : public PkFEMDimensionDispatch<GV,D,R,k,GV::dimension>
      , public warn_deprecated_pk_interface<d>

    {

    public:

      PkLocalFiniteElementMap(const GV& gv)
        : PkFEMDimensionDispatch<GV,D,R,k,GV::dimension>(gv)
        , warn_deprecated_pk_interface<d>()
      {}

    };


  }
}

#endif
