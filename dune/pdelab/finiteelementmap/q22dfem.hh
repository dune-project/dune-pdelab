// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_Q22DFEM_HH
#define DUNE_PDELAB_Q22DFEM_HH

#include<dune/localfunctions/lagrange/q22d.hh>

#include"finiteelementmap.hh"
#include <dune/pdelab/finiteelementmap/global.hh>
namespace Dune {
  namespace PDELab {

	//! wrap up element from local functions
    //! \ingroup FiniteElementMap
	template<class D, class R>
	class Q22DLocalFiniteElementMap
	  : public SimpleLocalFiniteElementMap< Dune::Q22DLocalFiniteElement<D,R> >
	{};

    //! Global-valued finite element map for Q22D elements
    /**
     * \ingroup FiniteElementMap
     *
     * \tparam Geometry Type of the geometry od the elements.
     * \tparam RF       Range field type.
     */
    template<class Geometry, class RF>
    class Q22DFiniteElementMap
      : public GeometryFiniteElementMap<
          Q22DFiniteElementFactory<Geometry, RF>
          >
    {
      typedef Q22DFiniteElementFactory<Geometry, RF> FEFactory;
      typedef GeometryFiniteElementMap<FEFactory> Base;

      static FEFactory feFactory;

    public:
      Q22DFiniteElementMap() : Base(feFactory) { }
    };

    template<class GV, class RF>
    typename Q22DFiniteElementMap<GV, RF>::FEFactory
    Q22DFiniteElementMap<GV, RF>::feFactory;
  }
}

#endif
