
#ifndef DUNE_GENEO_SUBDOMAINBASIS_HH
#define DUNE_GENEO_SUBDOMAINBASIS_HH

namespace Dune {
  namespace PDELab {
    /*!
     * \brief This represents a general per-subdomain basis.
     * \tparam X Vector type.
     */
    template <class X>
    class SubdomainBasis {

    public:
      SubdomainBasis() {}

      /*!
       * \brief Constructor creating a basis with only one basis function per subdomain.
       *
       */
      SubdomainBasis(X& basis_function) {
        this->local_basis.resize(1);
        this->local_basis[0] = std::make_shared<X>(basis_function);
      }

      std::vector<std::shared_ptr<X> > local_basis;

    };
  }
}



#endif //DUNE_GENEO_SUBDOMAINBASIS_HH
