#ifndef DUNE_PDELAB_BACKEND_ISTL_GENEO_NonoverlappingGenEOBasisFromFiles_HH
#define DUNE_PDELAB_BACKEND_ISTL_GENEO_NonoverlappingGenEOBasisFromFiles_HH

#include <dune/pdelab/backend/istl/geneo/subdomainbasis.hh>
#include <string>

namespace Dune {
  namespace PDELab {

    /*!
     */
    template<class GridView, class M, class X>
    class GenEOBasisFromFiles : public SubdomainBasis<X>
    {

    public:

      /*!
       * \brief Constructor of subdomainbasis from files where the
       * local EV from a pristine model are saved
       */
      GenEOBasisFromFiles(std::string& path_to_storage, int basis_size, int subdomain_number, int verbose = 0) {

        if (verbose > 1) std::cout << "Getting EV basis for subdomain: " << subdomain_number << " from offline." << std::endl;

        this->local_basis.resize(basis_size);

        for (int basis_index = 0; basis_index < basis_size; basis_index++) {
          std::shared_ptr<X> ev = std::make_shared<X>();
          std::string filename_EV = path_to_storage + std::to_string(subdomain_number) + "_EV_" + std::to_string(basis_index) + ".mm";
          std::ifstream file_EV;
          file_EV.open(filename_EV.c_str(), std::ios::in);
          Dune::readMatrixMarket(*ev,file_EV);
          file_EV.close();

          this->local_basis[basis_index] = ev;
        }
      }
    };
  }
}

#endif
