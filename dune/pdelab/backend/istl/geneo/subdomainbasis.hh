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

      /*!
       * \brief Write EV basis to a .txt file
       */
       void to_file(std::string basename) {
         std::ofstream output_basis_size;
         std::string filename_basis_size = basename + "_basis_size.txt";
         output_basis_size.open(filename_basis_size, std::ios::out);
         output_basis_size << local_basis.size();
         output_basis_size.close();

          for (int basis_index = 0; basis_index < local_basis.size(); basis_index++) {
            std::ostringstream os;
            os << basis_index;
            std::string filename = basename + "_EV_" + os.str() + ".txt";
            std::ofstream output;
            output.open(filename, std::ios::out);
            //if (!output.is_open()) {std::cout << "Error: Cannot open file" << std::endl;}
            //else {std::cout << "Writting " << filename << std::endl;}
            Dune::printvector(output, *local_basis[basis_index], "", "", 3, 10, 16);
            output.close();
          }
       }

    protected:
      std::vector<std::shared_ptr<X> > local_basis;

    };





  }
}



#endif //DUNE_PDELAB_BACKEND_ISTL_GENEO_SUBDOMAINBASIS_HH
