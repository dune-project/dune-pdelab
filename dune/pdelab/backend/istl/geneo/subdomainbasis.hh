
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
       * \brief Write EV basis in .mm format
       */
       void to_file(std::string basename, int rank) {
         // For each processor, printing the size of the associated subdomainbasis.
         // As this size is variable for each processor, it is the way to keep
         // this information after consecutive runs
         std::ofstream output_basis_size;
         std::ostringstream osrank;
         osrank << rank;
         std::string filename = basename + "_r" + osrank.str() + "_size.txt";
         output_basis_size.open(filename, std::ios::out);
         output_basis_size << local_basis.size();
         output_basis_size.close();

         // Writing subdomainbasis using matrixmarket format
         for (int basis_index = 0; basis_index < local_basis.size(); basis_index++) {
           std::ostringstream rfilename;
           rfilename << basename << "_r" << rank << "_" << basis_index << ".mm";
           std::ofstream file(rfilename.str().c_str());
           file.setf(std::ios::scientific,std::ios::floatfield);
           Dune::writeMatrixMarket(*local_basis[basis_index], file);
           file.close();
         }
       }

    protected:
      std::vector<std::shared_ptr<X> > local_basis;

    };





  }
}



#endif //DUNE_PDELAB_BACKEND_ISTL_GENEO_SUBDOMAINBASIS_HH
