// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_MONOMFEM_HH
#define DUNE_PDELAB_MONOMFEM_HH

#include <cstddef>

#include <dune/geometry/type.hh>

#include<dune/localfunctions/monom.hh>

#include"finiteelementmap.hh"
#include <dune/pdelab/finiteelementmap/global.hh>

namespace Dune {
  namespace PDELab {

	//! wrap up element from local functions
    //! \ingroup FiniteElementMap
	template<class D, class R, int d, int p>
	class MonomLocalFiniteElementMap
	  : public SimpleLocalFiniteElementMap< Dune::MonomLocalFiniteElement<D,R,d,p> >
	{
    public:
      MonomLocalFiniteElementMap (Dune::GeometryType::BasicType basicType) DUNE_DEPRECATED
        : SimpleLocalFiniteElementMap< Dune::MonomLocalFiniteElement<D,R,d,p> >(Dune::MonomLocalFiniteElement<D,R,d,p>(GeometryType(basicType,d)))
        , _gt(basicType,d)
      {
      }
      MonomLocalFiniteElementMap (const Dune::GeometryType& type)
        : SimpleLocalFiniteElementMap< Dune::MonomLocalFiniteElement<D,R,d,p> >(Dune::MonomLocalFiniteElement<D,R,d,p>(type)), _gt(type)
      {
      }

      bool fixedSize() const
      {
        return true;
      }

      std::size_t size(GeometryType gt) const
      {
        return gt == _gt ? Dune::MonomImp::Size<d,p>::val : 0;
      }

      std::size_t maxLocalSize() const
      {
        return MonomImp::Size<d,p>::val;
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
          MonomFiniteElementFactory<Geometry, RF, p>
          >
    {
      typedef MonomFiniteElementFactory<Geometry, RF, p> FEFactory;
      typedef GeometryFiniteElementMap<FEFactory> Base;

      static FEFactory feFactory;

    public:
      MonomFiniteElementMap() : Base(feFactory) { }
    };

    template<class GV, class RF, std::size_t p>
    typename MonomFiniteElementMap<GV, RF, p>::FEFactory
    MonomFiniteElementMap<GV, RF, p>::feFactory;

  }
}

#endif //DUNE_PDELAB_MONOMFEM_HH
