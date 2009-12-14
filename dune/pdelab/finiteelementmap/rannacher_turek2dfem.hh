// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_RANNACHER_TUREK2DFEM_HH
#define DUNE_PDELAB_RANNACHER_TUREK2DFEM_HH

#include<dune/localfunctions/rannacher_turek2d.hh>
#include"finiteelementmap.hh"

namespace Dune {
  namespace PDELab {

    //! wrap up element from local functions
    //! \ingroup FiniteElementMap
    template<class D, class R>
    class RannacherTurek2DLocalFiniteElementMap
      : public SimpleLocalFiniteElementMap<RannacherTurek2DLocalFiniteElement<D,R> >
    {};
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_RANNACHER_TUREK2DFEM_HH
