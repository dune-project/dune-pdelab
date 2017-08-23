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
      bool fixedSize() const
      {
        return true;
      }

      bool hasDOFs(int codim) const
      {
        return codim == 1;
      }

      std::size_t size(GeometryType gt) const
      {
        return gt.isLine() ? 1 : 0;
      }

      std::size_t maxLocalSize() const
      {
        return 4;
      }
    };
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_FINITEELEMENTMAP_RANNACHERTUREKFEM_HH
