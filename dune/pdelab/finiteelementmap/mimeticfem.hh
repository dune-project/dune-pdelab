// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_FINITEELEMENTMAP_MIMETICFEM_HH
#define DUNE_PDELAB_FINITEELEMENTMAP_MIMETICFEM_HH

#include<vector>
#include<dune/geometry/type.hh>
#include<dune/localfunctions/mimetic.hh>
#include"finiteelementmap.hh"

namespace Dune {
  namespace PDELab {

    template<typename IIS, typename D, typename R, int dim>
    class MimeticLocalFiniteElementMap :
      public Dune::PDELab::LocalFiniteElementMapInterface<Dune::PDELab::LocalFiniteElementMapTraits< MimeticLocalFiniteElement<D,R,dim> >,
                                                          MimeticLocalFiniteElementMap<IIS,D,R,dim> >
    {
      typedef MimeticLocalFiniteElement<D,R,dim> FE;

    public:
      //! \brief export type of the signature
      typedef Dune::PDELab::LocalFiniteElementMapTraits<FE> Traits;

      //! \brief Use when Imp has a standard constructor
      MimeticLocalFiniteElementMap (const IIS& iis_, Dune::GeometryType::BasicType basicType)
        : iis(iis_), bt(basicType)

      {
        // create a standard number of variants
        variant.resize(20);
        for (int i=0; i<20; i++) variant[i] = FE(bt,i);
      }

      //! \brief get local basis functions for entity
      template<class EntityType>
      const typename Traits::FiniteElementType& find (const EntityType& e) const
      {
        size_t n = static_cast<size_t>(iis.size(e));
        if (n<variant.size())
          return variant[n];
        else
          {
            size_t old_n = variant.size();
            variant.resize(n+1);
            for (size_t i=old_n; i<n+1; i++) variant[i] = FE(bt,i);
            return variant[n];
          }
      }

    private:
      const IIS& iis;
      Dune::GeometryType::BasicType bt;
      mutable std::vector<FE> variant;
    };

  }
}

#endif // DUNE_PDELAB_FINITEELEMENTMAP_MIMETICFEM_HH
