// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_FINITEELEMENTMAP_QKFEM_HH
#define DUNE_PDELAB_FINITEELEMENTMAP_QKFEM_HH

#include <cstddef>

#include <dune/localfunctions/lagrange/qk.hh>
#include <dune/pdelab/finiteelementmap/finiteelementmap.hh>

namespace Dune {
  namespace PDELab {

    //! wrap up element from local functions
    //! \ingroup FiniteElementMap
    template<typename GV, typename D, typename R, std::size_t k>
    class QkLocalFiniteElementMap
      : public SimpleLocalFiniteElementMap<Dune::QkLocalFiniteElement<D,R,GV::dimension,k>,GV::dimension>
    {

    public:

      static_assert(
        (k == 1) or (GV::dimension == 2 or GV::dimension == 3),
        "QkLocalFiniteElementMap with k = 2 is only implemented for d = 2,3"
        );

      static_assert(
        k == 1 or k == 2,
        "QkLocalFiniteElementMap is only implemented for k = 1,2"
        );

      QkLocalFiniteElementMap(const GV& gv)
      {}

      static constexpr bool fixedSize()
      {
        return true;
      }

      static constexpr bool hasDOFs(int codim)
      {
        switch(k)
          {
          case 1:
            return codim == GV::dimension;
          case 2:
            return true;
          default:
            return false;
          }
      }

      static constexpr std::size_t size(GeometryType gt)
      {
        switch (k)
          {
          case 1:
            return gt == GeometryTypes::vertex ? 1 : 0;
          case 2:
            // Q1 simply attaches a single DOF to each subentity
            return 1;
          default:
            return 0;
          }
      }

      static constexpr std::size_t maxLocalSize()
      {
        std::size_t r = 1;
        for (std::size_t i = 0; i < GV::dimension; ++i)
          r *= (k + 1);
        return r;
      }

    };

  }
}

#endif // DUNE_PDELAB_FINITEELEMENTMAP_QKFEM_HH
