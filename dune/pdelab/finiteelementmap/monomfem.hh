// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_FINITEELEMENTMAP_MONOMFEM_HH
#define DUNE_PDELAB_FINITEELEMENTMAP_MONOMFEM_HH

#include <cstddef>

#include <dune/geometry/type.hh>

#include<dune/localfunctions/monomial.hh>

#include"finiteelementmap.hh"
#include <dune/pdelab/finiteelementmap/global.hh>

namespace Dune {
  namespace PDELab {

	//! wrap up element from local functions
    //! \ingroup FiniteElementMap
	template<class D, class R, int d, int p>
	class MonomLocalFiniteElementMap
      : public SimpleLocalFiniteElementMap< Dune::MonomialLocalFiniteElement<D,R,d,p>,d>
	{
    public:

      MonomLocalFiniteElementMap (const Dune::GeometryType& type)
        : SimpleLocalFiniteElementMap< Dune::MonomialLocalFiniteElement<D,R,d,p>,d>(Dune::MonomialLocalFiniteElement<D,R,d,p>(type)), _gt(type)
      {
      }

      static constexpr bool fixedSize()
      {
        return true;
      }

      static constexpr bool hasDOFs(int codim)
      {
        return codim == 0;
      }

      std::size_t size(GeometryType gt) const
      {
        return gt == _gt ? MonomialLocalBasis<D,R,d,p>::size() : 0;
      }

      static constexpr std::size_t maxLocalSize()
      {
        return MonomialLocalBasis<D,R,d,p>::size();
      }

    private:
      const GeometryType _gt;

    };

    //! Global-valued finite element map for Monom elements
    /**
     * \ingroup FiniteElementMap
     *
     * \tparam Geometry Type of the geometry od the elements.
     * \tparam RF       Range field type.
     * \tparam p        Maximum polynomial order of the elements.
     */
    template<class Geometry, class RF, std::size_t p>
    class MonomFiniteElementMap
      : public GeometryFiniteElementMap<
          MonomialFiniteElementFactory<Geometry, RF, p>
          >
    {
      typedef MonomialFiniteElementFactory<Geometry, RF, p> FEFactory;
      typedef GeometryFiniteElementMap<FEFactory> Base;

      static FEFactory feFactory;

    public:

      //! The dimension of the finite elements returned by this map.
      static constexpr int dimension = Geometry::mydimension;

      MonomFiniteElementMap() : Base(feFactory) { }
    };

    template<class GV, class RF, std::size_t p>
    typename MonomFiniteElementMap<GV, RF, p>::FEFactory
    MonomFiniteElementMap<GV, RF, p>::feFactory;

  }
}

#endif // DUNE_PDELAB_FINITEELEMENTMAP_MONOMFEM_HH
