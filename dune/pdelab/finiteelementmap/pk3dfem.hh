// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_FINITELEMENTMAP_HH_PK3DFEM_HH
#define DUNE_PDELAB_FINITELEMENTMAP_HH_PK3DFEM_HH

#warning dune/pdelab/finiteelementmap/pk3dfem.hh and Pk3DLocalFiniteElementMap are deprecated, please use dune/pdelab/finiteelementmap/pkfem.hh and PkLocalFiniteElementMap instead

#include <dune/common/deprecated.hh>

#include <dune/pdelab/finiteelementmap/pkfem.hh>

namespace Dune
{
  namespace PDELab
  {

    template<typename GV, typename D, typename R, unsigned int k>
    class DUNE_DEPRECATED_MSG("Please use PkLocalFiniteElementMap instead") Pk3DLocalFiniteElementMap
      : public fem::PkLocalFiniteElementMapBase<GV,D,R,k,3>
    {

    public:

      Pk3DLocalFiniteElementMap(const GV& gv)
        : fem::PkLocalFiniteElementMapBase<GV,D,R,k,3>(gv)
      {}

  }
}

#endif // DUNE_PDELAB_FINITELEMENTMAP_HH_PK3DFEM_HH
