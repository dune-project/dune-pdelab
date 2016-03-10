// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_FINITEELEMENT_LOCALBASISCACHE_HH
#define DUNE_PDELAB_FINITEELEMENT_LOCALBASISCACHE_HH

#include<vector>
#include<map>

#include<dune/common/exceptions.hh>

namespace Dune {
  namespace PDELab {

    //! \brief store values of basis functions and gradients in a cache
    template<class LocalBasisType>
    class LocalBasisCache
    {
      typedef typename LocalBasisType::Traits::DomainFieldType DomainFieldType;
      typedef typename LocalBasisType::Traits::DomainType DomainType;
      typedef typename LocalBasisType::Traits::RangeType RangeType;
      typedef typename LocalBasisType::Traits::JacobianType JacobianType;

      struct less_than
      {
        bool operator() (const DomainType& v1, const DomainType& v2) const
        {
          for (typename DomainType::size_type i=0; i<DomainType::dimension; i++)
            {
              if ( v1[i] < v2[i]-1e-5 ) return true;   // is less than
              if ( v1[i] > v2[i]+1e-5 ) return false;  // is greater than
            }
          return false; // is equal
        }
      };

      typedef std::map<DomainType,std::vector<RangeType>,less_than> FunctionCache;
      typedef std::map<DomainType,std::vector<JacobianType>,less_than> JacobianCache;

    public:

      //! \brief constructor
      LocalBasisCache () {}

      //! evaluate basis functions at a point
      const std::vector<RangeType>&
      evaluateFunction (const DomainType& position, const LocalBasisType& localbasis) const
      {
        typename FunctionCache::iterator it = functioncache.find(position);
        if (it!=functioncache.end()) return it->second;
        std::vector<RangeType> values;
        localbasis.evaluateFunction(position,values);
        it = functioncache.insert(functioncache.begin(),std::pair<DomainType,std::vector<RangeType> >(position,values));
        return it->second;
      }

      //! evaluate Jacobians at a point
      const std::vector<JacobianType>&
      evaluateJacobian (const DomainType& position, const LocalBasisType& localbasis) const
      {
        typename JacobianCache::iterator it = jacobiancache.find(position);
        if (it!=jacobiancache.end()) return it->second;
        std::vector<JacobianType> values;
        localbasis.evaluateJacobian(position,values);
        it = jacobiancache.insert(jacobiancache.begin(),std::pair<DomainType,std::vector<JacobianType> >(position,values));
        return it->second;
      }

    private:
      mutable FunctionCache functioncache;
      mutable JacobianCache jacobiancache;
    };

  }
}

#endif // DUNE_PDELAB_FINITEELEMENT_LOCALBASISCACHE_HH
