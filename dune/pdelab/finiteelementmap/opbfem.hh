// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_OPBFEM_HH
#define DUNE_PDELAB_OPBFEM_HH

#include<dune/pdelab/finiteelementmap/finiteelementmap.hh>
#include"l2orthonormal.hh"

namespace Dune {
  namespace PDELab {

    //! wrap up element from local functions
    //! \ingroup FiniteElementMap
    template<class D, class R, int k, int d, Dune::GeometryType::BasicType bt, typename ComputationFieldType=R, PB::BasisType basisType = PB::BasisType::Pk>
    class OPBLocalFiniteElementMap
      : public Dune::PDELab::SimpleLocalFiniteElementMap< Dune::OPBLocalFiniteElement<D,R,k,d,bt,ComputationFieldType,basisType> >
    {
      typedef PB::BasisTraits<basisType> BasisTraits;
    public:

      bool fixedSize() const
      {
        return true;
      }

      std::size_t size(GeometryType gt) const
      {
        if (gt == GeometryType(bt,d))
          return BasisTraits::template Size<k,d>::value;
        else
          return 0;
      }

      std::size_t maxLocalSize() const
      {
        return BasisTraits::template Size<k,d>::value;
      }

    };

  }
}

#endif
