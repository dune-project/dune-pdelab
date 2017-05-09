// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// Qk DG basis with Gauss Lobatto points

#ifndef DUNE_PDELAB_FINITEELEMENTMAP_QKDGGL_HH
#define DUNE_PDELAB_FINITEELEMENTMAP_QKDGGL_HH

#warning This file is deprecated and will be removed after the Dune-PDELab 2.5 release! Switch to the class QkDGLocalFiniteElement map with QkDGBasisPolynomial::lobatto in the header dune/pdelab/finiteelementmap/qkdg.hh.

#include <dune/pdelab/finiteelement/qkdglobatto.hh>
#include <dune/pdelab/finiteelementmap/finiteelementmap.hh>

namespace Dune {
  namespace PDELab {

    template<class D, class R, int k, int d>
    using QkDGGLLocalFiniteElementMap DUNE_DEPRECATED_MSG("This class is deprecated and will be removed after the Dune-PDELab 2.5 release! Switch to the class QkDGLocalFiniteElement map with QkDGBasisPolynomial::lobatto in the header dune/pdelab/finiteelementmap/qkdg.hh.") = QkDGLocalFiniteElementMap<D,R,k,d,QkDGBasisPolynomial::lobatto>;

  }
}

#endif // DUNE_PDELAB_FINITEELEMENTMAP_QKDGGL_HH
