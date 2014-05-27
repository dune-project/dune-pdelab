// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_RANNACHER_TUREK2DFEM_HH
#define DUNE_PDELAB_RANNACHER_TUREK2DFEM_HH

#warning dune/pdelab/finiteelementmap/rannacher_turek2dfem.hh and RannacherTurek2DLocalFiniteElementMap are deprecated, please use dune/pdelab/finiteelementmap/rannacherturekfem.hh and RannacherTurekLocalFiniteElementMap instead

#include<dune/localfunctions/rannacherturek.hh>
#include"finiteelementmap.hh"

namespace Dune {
  namespace PDELab {

    //! wrap up element from local functions
    //! \ingroup FiniteElementMap
    template<class D, class R>
    class DUNE_DEPRECATED_MSG("Please use RannacherTurekLocalFiniteElementMap instead") RannacherTurek2DLocalFiniteElementMap
      : public SimpleLocalFiniteElementMap<RannacherTurekLocalFiniteElement<D,R,2> >
    {
    public:
      bool fixedSize() const
      {
        return true;
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

#endif // DUNE_PDELAB_RANNACHER_TUREK2DFEM_HH
