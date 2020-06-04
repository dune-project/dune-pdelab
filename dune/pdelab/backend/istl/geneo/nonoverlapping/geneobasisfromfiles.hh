#ifndef DUNE_PDELAB_BACKEND_ISTL_GENEO_NonoverlappingGenEOBasisFromFiles_HH
#define DUNE_PDELAB_BACKEND_ISTL_GENEO_NonoverlappingGenEOBasisFromFiles_HH

#include <dune/pdelab/backend/istl/geneo/subdomainbasis.hh>
#include <string>

namespace Dune {
  namespace PDELab {

    /*!
     */
    template<class GridView, class M, class X>
    class NonoverlappingGenEOBasisFromFiles : public SubdomainBasis<X>
    {

    public:

      /*!
       * \brief Constructor of subdomainbasis from files of the pristine solution
       */
      NonoverlappingGenEOBasisFromFiles(NonoverlappingOverlapAdapter<GridView, X, M>& adapter) {
        std::cout << "Getting EV basis at proc " << adapter.gridView().comm().rank() << " from offline." << std::endl;
        std::ostringstream os1;
        os1 << adapter.gridView().comm().rank();
        std::string basename = "Offline/Proc_"+os1.str();

        // Get the basis size from file
        std::string filename_basis_size = basename+"_basis_size.txt";
        std::ifstream input_basis_size;
        input_basis_size.open(filename_basis_size, std::ios::in);
        //if (!input_basis_size.is_open()) {std::cout << "Error: Cannot open file" << std::endl;}
        //else {std::cout << "Reading " << filename_basis_size << std::endl;}
        int basis_size;
        input_basis_size >> basis_size;
        input_basis_size.close();

        this->local_basis.resize(basis_size);

        for (int basis_index = 0; basis_index < this->local_basis.size(); basis_index++) {

          std::shared_ptr<X> ev = std::make_shared<X>();

          std::ifstream input;
          std::ostringstream os;
          os << basis_index;
          std::string filename = basename + "_EV_" + os.str() + ".txt";
          input.open(filename, std::ios::in);

          std::string line;
          std::getline(input, line);
          auto start = line.find("blocks=")+7U;
          auto end = line.find(",", start);
          int EV_size = std::stoi(line.substr(start, end - start));
          ev->resize(EV_size);
          for(int i=0; i<ev->size();i++){
            int tmp;
            input >> tmp >> (*ev)[i][0] >> (*ev)[i][1] >> (*ev)[i][2];
          }

          input.close();

          this->local_basis[basis_index] = ev;
        }


      }

    };


  }
}

#endif
