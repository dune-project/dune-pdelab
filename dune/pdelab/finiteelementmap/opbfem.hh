// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_OPBFEM_HH
#define DUNE_PDELAB_OPBFEM_HH

#include<dune/pdelab/finiteelementmap/finiteelementmap.hh>
#include"l2orthonormal.hh"

namespace Dune {
  namespace PDELab {

    //! wrap up element from local functions
    //! \ingroup FiniteElementMap
    template<class D, class R, int k, int d, Dune::GeometryType::BasicType bt, typename ComputationFieldType=R>
    class OPBLocalFiniteElementMap
      : public Dune::PDELab::SimpleLocalFiniteElementMap< Dune::OPBLocalFiniteElement<D,R,k,d,bt,ComputationFieldType> >
    {

      static const std::size_t _per_cell_size =
        (bt == GeometryType::cube) ? Dune::PB::QkSize<k,d>::value
        : (bt == GeometryType::simplex) ? Dune::PB::PkSize<k,d>::value
        : 0;

    public:

      bool fixedSize() const
      {
        return true;
      }

      std::size_t size(GeometryType gt) const
      {
        if (gt == GeometryType(bt,d))
          return _per_cell_size;
        else
          return 0;
      }

      std::size_t maxLocalSize() const
      {
        return _per_cell_size;
      }

    };

  }
}

#endif
