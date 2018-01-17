#ifndef DUNE_GENEO_ONEDIMSUBDOMAINBASIS_HH
#define DUNE_GENEO_ONEDIMSUBDOMAINBASIS_HH

#include "subdomainbasis.hh"

namespace Dune {
  namespace PDELab {

    /*!
    * \brief A local basis containing only one basis function.
    *
    * This is primarily intended for constructing a basis
    * exclusively containing a partition of unity.
    */
    template<class X>
    class OneDimensionalSubdomainBasis : public SubdomainBasis<X>
    {

    public:
      OneDimensionalSubdomainBasis(X& part_unity) {
        this->local_basis.resize(1);
        this->local_basis[0] = std::make_shared<X>(part_unity);
      }

    };

  }
}
#endif //DUNE_GENEO_ONEDIMSUBDOMAINBASIS_HH
