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
       * \brief Constructor of subdomainbasis from files where the
       * local EV from a pristine model are saved
       */
      NonoverlappingGenEOBasisFromFiles(NonoverlappingOverlapAdapter<GridView, X, M>& adapter, std::string& basename) {

        std::cout << "Getting EV basis at proc " << adapter.gridView().comm().rank() << " from offline." << std::endl;

        std::ostringstream osrank;
        osrank << adapter.gridView().comm().rank();
        int basis_size;

        // Get the basis size from file
        std::string filename_basis_size = basename+ "_r" + osrank.str() + "_size.txt";
        std::ifstream input_basis_size;
        input_basis_size.open(filename_basis_size, std::ios::in);
        if (!input_basis_size.is_open()) {std::cout << "Error: Cannot open file " << filename_basis_size << std::endl;}
        input_basis_size >> basis_size;
        input_basis_size.close();
        this->local_basis.resize(basis_size);

        for (int basis_index = 0; basis_index < this->local_basis.size(); basis_index++) {

          std::shared_ptr<X> ev = std::make_shared<X>();

          std::ostringstream rfilename;
          rfilename<< basename << "_r" << adapter.gridView().comm().rank() << "_" << basis_index << ".mm";
          std::ifstream file;
          file.open(rfilename.str().c_str(), std::ios::in);
          if(!file)
            DUNE_THROW(IOError, "Could not open file: " << rfilename.str().c_str());
          Dune::readMatrixMarket(*ev,file);
          file.close();
          /*std::ifstream input;
          std::ostringstream os;
          os << basis_index;
          std::string filename = basename + "_" + os.str() + ".txt";
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

          input.close();*/

          this->local_basis[basis_index] = ev;
        }


      }

    };


  }
}

#endif
