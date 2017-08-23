// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// Pk in one dimension with k as runtime variable

#ifndef DUNE_PDELAB_FINITEELEMENTMAP_PK1D_HH
#define DUNE_PDELAB_FINITEELEMENTMAP_PK1D_HH

#include<dune/pdelab/finiteelementmap/finiteelementmap.hh>
#include<dune/pdelab/finiteelement/pk1d.hh>

namespace Dune {

  namespace PDELab {

    /** \brief FiniteElementMap for the Pk basis in 1d
     *
     *  \note k is a runtime varbiable.
     *
     *  \tparam D Type to represent domain.
     *  \tparam R Type to represent range.
     */
    template<class D, class R>
    class Pk1dLocalFiniteElementMap
      : public Dune::PDELab::SimpleLocalFiniteElementMap<Pk1dLocalFiniteElement<D,R>,1>
    {
    public:

      Pk1dLocalFiniteElementMap (std::size_t k)
        : Dune::PDELab::SimpleLocalFiniteElementMap<Pk1dLocalFiniteElement<D,R>,1>(Pk1dLocalFiniteElement<D,R>(k))
        , _k(k)
      {}

      static constexpr bool fixedSize()
      {
        return true;
      }

      bool hasDOFs(int codim) const
      {
        switch (codim)
          {
          case 0:
            return _k != 1;
          case 1:
            return _k > 0;
          }
        return false;
      }

      std::size_t size(GeometryType gt) const
      {
        if (gt.isVertex())
          return _k > 0 ? 1 : 0;
        if (gt.isLine())
          return _k > 0 ? _k - 1 : 1;
        return 0;
      }

      std::size_t maxLocalSize() const
      {
        return _k + 1;
      }

    private:
      const std::size_t _k;
    };
  }
}
#endif // DUNE_PDELAB_FINITEELEMENTMAP_PK1D_HH
