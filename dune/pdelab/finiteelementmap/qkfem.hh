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

      QkLocalFiniteElementMap(const GV& gv)
      {}

      bool fixedSize() const
      {
        return true;
      }

      bool hasDOFs(int codim) const
      {
        switch(k)
          {
          case 1:
            return codim == GV::dimension;
          case 2:
            if (GV::dimension != 2 && GV::dimension != 3)
              DUNE_THROW(NotImplemented,"QkLocalFiniteElementMap with k = 2 is only implemented for d = 2,3");
            return 1;
          default:
            DUNE_THROW(NotImplemented,"QkLocalFiniteElementMap is only implemented for k <= 2");
          }
      }

      std::size_t size(GeometryType gt) const
      {
        switch (k)
          {
          case 1:
            return gt.isVertex() ? 1 : 0;
          case 2:
            {
              if (GV::dimension != 2 && GV::dimension != 3)
                DUNE_THROW(NotImplemented,"QkLocalFiniteElementMap with k = 2 is only implemented for d = 2,3");
              // Q1 simply attaches a single DOF to each subentity
              return 1;
            }
          default:
            DUNE_THROW(NotImplemented,"QkLocalFiniteElementMap is only implemented for k <= 2");
          }
      }

      std::size_t maxLocalSize() const
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
