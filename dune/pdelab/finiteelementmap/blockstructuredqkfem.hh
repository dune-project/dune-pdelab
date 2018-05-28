// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_MYPOISSON_QKFEM_REORDER_HH
#define DUNE_MYPOISSON_QKFEM_REORDER_HH

#include <dune/pdelab/finiteelementmap/finiteelementmap.hh>
#include <dune/localfunctions/blockstructuredqk.hh>
#include <dune/localfunctions/lagrange/qk.hh>

namespace Dune {
  namespace Blockstructured {

    template<typename D, typename R, int d, int k_, int blocks_>
    struct QklocalFiniteElementWrapper
        : public Dune::QkLocalFiniteElement<D, R, d, k_ * blocks_> {
      static constexpr int k = k_;
      static constexpr int blocks = blocks_;
    };
  }

  namespace PDELab {

    //! wrap up element from local functions
    //! \ingroup FiniteElementMap
    template<typename GV, typename D, typename R, std::size_t k, std::size_t blocks>
    class BlockstructuredQkLocalFiniteElementMap
        : public SimpleLocalFiniteElementMap
            <Blockstructured::QklocalFiniteElementWrapper<D, R, GV::dimension, k, blocks>, GV::dimension> {

      static constexpr std::size_t DOFs1d = k * blocks + 1;

    public:

      BlockstructuredQkLocalFiniteElementMap(const GV &gv) {}

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