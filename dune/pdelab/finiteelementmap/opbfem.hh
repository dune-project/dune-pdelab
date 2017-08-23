// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_FINITEELEMENTMAP_OPBFEM_HH
#define DUNE_PDELAB_FINITEELEMENTMAP_OPBFEM_HH

#include<dune/pdelab/finiteelementmap/finiteelementmap.hh>
#include<dune/pdelab/finiteelement/l2orthonormal.hh>

namespace Dune {
  namespace PDELab {

    //! wrap up element from local functions
    //! \ingroup FiniteElementMap
    template<class D, class R, int k, int d, Dune::GeometryType::BasicType bt, typename ComputationFieldType=R, PB::BasisType basisType = PB::BasisType::Pk>
    class OPBLocalFiniteElementMap
      : public Dune::PDELab::SimpleLocalFiniteElementMap<Dune::OPBLocalFiniteElement<D,R,k,d,bt,ComputationFieldType,basisType>,d>
    {
      typedef PB::BasisTraits<basisType> BasisTraits;
    public:

      static constexpr bool fixedSize()
      {
        return true;
      }

      static constexpr bool hasDOFs(int codim)
      {
        return codim == 0;
      }

      static constexpr std::size_t size(GeometryType gt)
      {
        if (gt == GeometryType(bt,d))
          return BasisTraits::template Size<k,d>::value;
        else
          return 0;
      }

      static constexpr std::size_t maxLocalSize()
      {
        return BasisTraits::template Size<k,d>::value;
      }

    };

  }
}

#endif // DUNE_PDELAB_FINITEELEMENTMAP_OPBFEM_HH
