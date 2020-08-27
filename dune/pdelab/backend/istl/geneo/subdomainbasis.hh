
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
       * \brief Append a vector into the basis
       */
      void append(X& basis_function) {
        this->local_basis.push_back(std::make_shared<X>(basis_function));
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

      /*!
       * \brief Write EV basis in .mm format
       */
      void to_file(std::string basename, int rank) {
        // For each processor, printing the size of the associated subdomainbasis.
        // As this size is variable for each processor, it is the way to keep
        // this information after consecutive runs
        std::ofstream output_basis_size;
        std::string filename = basename + "_" + std::to_string(rank) + "_size.txt";
        output_basis_size.open(filename, std::ios::out);
        output_basis_size << local_basis.size();
        output_basis_size.close();

        // Writing subdomainbasis using matrixmarket format
        for (int basis_index = 0; basis_index < local_basis.size(); basis_index++) {
          std::string filename = basename + "_" + std::to_string(basis_index) + "_" + std::to_string(rank) + ".mm";
          Dune::storeMatrixMarket(*local_basis[basis_index], filename, 15);
        }
      }

    protected:
      std::vector<std::shared_ptr<X> > local_basis;

    };
  }
}



#endif //DUNE_PDELAB_BACKEND_ISTL_GENEO_SUBDOMAINBASIS_HH
