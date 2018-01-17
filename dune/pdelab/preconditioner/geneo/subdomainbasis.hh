
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

      std::vector<std::shared_ptr<X> > local_basis;

    };
  }
}



#endif //DUNE_GENEO_SUBDOMAINBASIS_HH
