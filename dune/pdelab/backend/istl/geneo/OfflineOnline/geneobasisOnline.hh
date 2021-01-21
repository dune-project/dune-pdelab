#ifndef DUNE_PDELAB_BACKEND_ISTL_GENEO_Nonoverlapping_GenEOBasisOnline_HH
#define DUNE_PDELAB_BACKEND_ISTL_GENEO_Nonoverlapping_GenEOBasisOnline_HH

#include <algorithm>
#include <functional>

#include <dune/pdelab/backend/istl/geneo/subdomainbasis.hh>
#include <dune/pdelab/backend/istl/geneo/arpackpp_geneo.hh>

#if HAVE_ARPACKPP

namespace Dune {
  namespace PDELab {

    /*!
     * \brief Implementation of the GenEO coarse basis
     * See Spillane et al., 2014, 'Abstract robust coarse spaces for systems of PDEs via generalized eigenproblems in the overlaps'
     * This coarse space is based on generalized eigenpoblems defined on the full stiffness matrix of a subdomain and one assembled only
     * on the area where this subdomain overlaps with others.
     */
    template<class GO, class Matrix, class Vector>
    class GenEOBasisOnline : public SubdomainBasis<Vector>
    {
      using GFS = typename GO::Traits::TrialGridFunctionSpace;
      using GridView = typename GFS::Traits::GridView;
      using GV = GridView;
      using V = typename GO::Domain;
      using M = typename GO::Jacobian;

    public:


      /**
       * @brief Constructor
       *
       * @param adapter Adapter for extending matrices by virtual overlap
       * @param A_extended Neumann matrix extended by virtual overlap
       * @param A_ovlp_extended Neumann matrix restricted to overlap region and extended by virtual overlap
       * @param part_unity Partition of unity extended by virtual overlap
       * @param eigenvalue_threshold Threshold up to which eigenvalues should be included in the GenEO space; no threshold is applied if value negative
       * @param nev Number of eigenvectors to be included in coarse space; returns chosen value in case of thresholding
       * @param nev_arpack Number of eigenvectors to be computed by ARPACK, possibly a larger value for stability
       * @param shift Shift value for shift invert mode solution of eigenproblem
       * @param add_part_unity Whether to add the partition of unity to the coarse space
       */
      GenEOBasisOnline(Matrix& A, Vector& part_unity, const double eigenvalue_threshold,
                    int& nev, int nev_arpack = -1, const double shift = 0.001, const bool add_part_unity = false,
                    const int verbose = 0) {

        const int v_size = A.M();

        if (nev_arpack == -1)
          nev_arpack = std::max(nev, 2);
        if (nev_arpack < nev)
          DUNE_THROW(Dune::Exception,"nev_arpack is less then nev!");

        // X * A_0 * X
        Matrix ovlp_mat(A);
        using Dune::PDELab::Backend::native;
        const int block_size = Vector::block_type::dimension;
        for (auto row_iter = native(ovlp_mat).begin(); row_iter != native(ovlp_mat).end(); row_iter++) {
          for (auto col_iter = row_iter->begin(); col_iter != row_iter->end(); col_iter++) {
            for (int i = 0; i < block_size; i++) {
              for (int j = 0; j < block_size; j++) {
                (*col_iter)[i][j] *= native(part_unity)[row_iter.index()][i] * native(part_unity)[col_iter.index()][j];
              }
            }
          }
        }

        // Setup Arpack for solving generalized eigenproblem
        std::cout << "ARPACK setup...";
        ArpackGeneo::ArPackPlusPlus_Algorithms<Matrix, Vector> arpack(A);
        std::cout << " done" << std::endl;
        double eps = 0.0;

        std::vector<double> eigenvalues(nev_arpack,0.0);
        std::vector<Vector> eigenvectors(nev_arpack,Vector(v_size));

        std::cout << "ARPACK solve...";
        arpack.computeGenNonSymMinMagnitude(ovlp_mat, eps, eigenvectors, eigenvalues, shift);
        std::cout << " done" << std::endl;

        // Count eigenvectors below threshold
        int cnt = -1;
        if (eigenvalue_threshold >= 0) {
          for (int i = 0; i < nev; i++) {
            if (eigenvalues[i] > eigenvalue_threshold) {
              cnt = i;
              break;
            }
          }
          if (verbose > 0)
            std::cout << "Subdomain eigen basis recomputed picked " << cnt << " eigenvectors" << std::endl;
          if (cnt == -1)
            DUNE_THROW(Dune::Exception,"No eigenvalue above threshold - not enough eigenvalues computed!");
        } else {
          cnt = nev;
        }

        // Write results to basis
        this->local_basis.resize(cnt);
        for (int base_id = 0; base_id < cnt; base_id++) {
          this->local_basis[base_id] = std::make_shared<Vector>(v_size);
          *this->local_basis[base_id] = eigenvectors[base_id];
          for (int i = 0; i < part_unity.N(); i++) {
            for (int j = 0; j < block_size; j++) {
              (*(this->local_basis[base_id]))[i][j] *= part_unity[i][j];
            }
          }
        }

        // Normalize basis vectors
        for (auto& v : this->local_basis) {
          *v *= 1.0 / v->two_norm2();
        }

        // Optionally add partition of unity to eigenvectors
        // Only if there is no near-zero eigenvalue (that usually already corresponds to a partition of unity!)
        if (add_part_unity && eigenvalues[0] > 1E-10) {
          this->local_basis.insert (this->local_basis.begin(), std::make_shared<Vector>(part_unity));
          this->local_basis.pop_back();
        }
      }
    };

    template<class GridView, class M, class Vector, typename vector1i>
    class GenEOBasisFromFiles : public SubdomainBasis<Vector>
    { // For testing

    public:

      GenEOBasisFromFiles(std::string& path_to_storage, int basis_size, int subdomain_number, vector1i& indiceChange, std::string defectID="", int verbose = 0) {

        if (verbose > 1) std::cout << "Getting EV basis for subdomain: " << subdomain_number << " from offline." << std::endl;

        this->local_basis.resize(basis_size);

        for (int basis_index = 0; basis_index < basis_size; basis_index++) {
          std::shared_ptr<Vector> ev = std::make_shared<Vector>();
          std::string filename_EV = path_to_storage + std::to_string(subdomain_number) + "_" + defectID + "EV_" + std::to_string(basis_index) + ".mm";
          std::ifstream file_EV;
          file_EV.open(filename_EV.c_str(), std::ios::in);
          Vector tmp;
          Dune::readMatrixMarket(tmp,file_EV);
          file_EV.close();

          ev->resize(tmp.N());
          for (int j=0; j<tmp.N(); j++){
            (*ev)[indiceChange[j]] = tmp[j];
          }

          this->local_basis[basis_index] = ev;
        }
      }
    };

    template<class GridView, class M, class Vector, class vector1i>
    class NeighbourBasis : public SubdomainBasis<Vector>
    {

    public:

      NeighbourBasis(std::string& path_to_storage, int basis_size, int subdomain_number, vector1i& offlineDoF2GI, vector1i& offlineNeighbourDoF2GI, std::string defectID="", int verbose = 0) {

        if (verbose > 1) std::cout << "Getting EV basis for neighbour subdomain: " << subdomain_number << " from offline." << std::endl;

        this->local_basis.resize(basis_size);

        std::vector<std::pair<int, int>> N2T;
        for(int i=0; i<offlineNeighbourDoF2GI.size();i++){
          auto it = std::find(offlineDoF2GI.begin(), offlineDoF2GI.end(), offlineNeighbourDoF2GI[i]);
          if (it != offlineDoF2GI.end()){
            N2T.push_back(std::make_pair(i,std::distance(offlineDoF2GI.begin(),it)));
            // std::cout << i << " : " << std::distance(offlineDoF2GI.begin(),it) << std::endl;
          }
        }

        for (int basis_index = 0; basis_index < basis_size; basis_index++) {
          Vector ev;
          std::string filename_EV = path_to_storage + std::to_string(subdomain_number) + "_" + defectID + "EV_" + std::to_string(basis_index) + ".mm";
          std::ifstream file_EV;
          file_EV.open(filename_EV.c_str(), std::ios::in);
          Dune::readMatrixMarket(ev,file_EV);
          file_EV.close();

          std::shared_ptr<Vector> ev_moved = std::make_shared<Vector>();

          ev_moved->resize(offlineDoF2GI.size());
          for(int i=0; i<N2T.size();i++){
            (*ev_moved)[N2T[i].second] = ev[N2T[i].first];
          }

          this->local_basis[basis_index] = ev_moved;
        }
      }
    };
  }
}

#endif

#endif
