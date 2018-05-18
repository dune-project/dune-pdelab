// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_MYPOISSON_QKFEM_REORDER_HH
#define DUNE_MYPOISSON_QKFEM_REORDER_HH

#include <dune/pdelab/finiteelementmap/finiteelementmap.hh>
#include <dune/pdelab/finiteelementmap/blockstructured/qk.hh>

namespace Dune {
  namespace Blockstructured {

    //! wrap up element from local functions
    //! \ingroup FiniteElementMap
    template<typename GV, typename D, typename R, std::size_t k, std::size_t blocks>
    class QkLocalFiniteElementMap
        : public Dune::PDELab::SimpleLocalFiniteElementMap
            <QkLocalFiniteElement<D, R, GV::dimension, k, blocks>, GV::dimension> {

      static constexpr std::size_t DOFs1d = k * blocks + 1;

    public:

      QkLocalFiniteElementMap(const GV &gv) {}

      bool fixedSize() const {
        return true;
      }

      bool hasDOFs(int codim) const {
        if (DOFs1d == 1)
          return codim == 0;
        else if (DOFs1d == 2)
          return codim == GV::dimension;
        else
          return 1;
      }

      std::size_t size(GeometryType gt) const {
        std::size_t acc = 1;
        for (std::size_t i = 0; i < gt.dim(); ++i)
          acc *= DOFs1d - 2;
        return acc;
      }

      std::size_t maxLocalSize() const {
        return Dune::StaticPower<DOFs1d, GV::dimension>::power;
      }

    };

  }
}
#endif //DUNE_MYPOISSON_QKFEM_HH
