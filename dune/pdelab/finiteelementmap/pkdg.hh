// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_FINITEELEMENTMAP_PKDG_HH
#define DUNE_PDELAB_FINITEELEMENTMAP_PKDG_HH

#include <dune/common/deprecated.hh>
#include <dune/common/exceptions.hh>
#include <dune/geometry/type.hh>
#include <dune/localfunctions/lagrange/pk.hh>

#include <dune/pdelab/finiteelementmap/finiteelementmap.hh>
#include <dune/pdelab/finiteelement/pkdglagrange/pkdglagrange.hh>


namespace Dune {
  namespace PDELab {

      template<typename GV, typename D, typename R, unsigned int k, unsigned int d>
      class PkDGLocalFiniteElementMap
        : public SimpleLocalFiniteElementMap<Dune::PkDGLocalFiniteElement<D,R,d,k>,d>
      {

      public:

        PkDGLocalFiniteElementMap()
        {}

        static constexpr bool fixedSize()
        {
          return true;
        }

        static constexpr bool hasDOFs(int codim)
        {
        return codim == 0;
        }

        static constexpr std::size_t maxLocalSize()
        {
          return Dune::PkStuff::PkDGSize<k, d>::value;
        }

        static constexpr std::size_t size(GeometryType gt)
        {
          if (gt == GeometryTypes::simplex(d))
            return Dune::PkStuff::PkDGSize<k, d>::value;
          else
            return 0;
        }
              //! return order of polynomial basis
        static constexpr std::size_t order()
        {
          return k;
        }
      };
  }
}

#endif // DUNE_PDELAB_FINITEELEMENTMAP_PKDG_HH
