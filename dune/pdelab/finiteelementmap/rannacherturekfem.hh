// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_FINITEELEMENTMAP_RANNACHERTUREKFEM_HH
#define DUNE_PDELAB_FINITEELEMENTMAP_RANNACHERTUREKFEM_HH

#include<dune/localfunctions/rannacherturek.hh>
#include"finiteelementmap.hh"

namespace Dune {
  namespace PDELab {

    //! wrap up element from local functions
    //! \ingroup FiniteElementMap
    template<class D, class R, std::size_t d>
    class  RannacherTurekLocalFiniteElementMap
      : public SimpleLocalFiniteElementMap<RannacherTurekLocalFiniteElement<D,R,d>,d>
    {
    public:

      static constexpr bool fixedSize()
      {
        return true;
      }

      static constexpr bool hasDOFs(int codim)
      {
        return codim == 1;
      }

      static constexpr std::size_t size(GeometryType gt)
      {
        return gt == GeometryTypes::line ? 1 : 0;
      }

      static constexpr std::size_t maxLocalSize()
      {
        return 4;
      }

    };
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_FINITEELEMENTMAP_RANNACHERTUREKFEM_HH
