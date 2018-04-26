
#ifndef DUNE_PDELAB_BACKEND_ISTL_GENEO_SUBDOMAINBASIS_HH
#define DUNE_PDELAB_BACKEND_ISTL_GENEO_SUBDOMAINBASIS_HH

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

      typedef X VectorType;

      /*!
       * \brief Constructor creating a basis with only one basis function per subdomain.
       *
       */
      SubdomainBasis(X& basis_function) {
        this->local_basis.resize(1);
        this->local_basis[0] = std::make_shared<X>(basis_function);
      }

      /*!
       * \brief Returns basis vector i
       */
      std::shared_ptr<X> get_basis_vector(int i) {
        return local_basis[i];
      }

      /*!
       * \brief Size of the local basis
       */
      int basis_size() {
        return local_basis.size();
      }

    protected:
      std::vector<std::shared_ptr<X> > local_basis;

    };
  }
}



#endif //DUNE_PDELAB_BACKEND_ISTL_GENEO_SUBDOMAINBASIS_HH
