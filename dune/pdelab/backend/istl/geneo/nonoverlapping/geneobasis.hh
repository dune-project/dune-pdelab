#ifndef DUNE_PDELAB_BACKEND_ISTL_GENEO_NonoverlappingGenEOBasis_HH
#define DUNE_PDELAB_BACKEND_ISTL_GENEO_NonoverlappingGenEOBasis_HH

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
    class NonoverlappingGenEOBasis : public SubdomainBasis<Vector>
    {
      typedef Dune::PDELab::Backend::Native<Matrix> ISTLM;
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
      NonoverlappingGenEOBasis(NonoverlappingOverlapAdapter<GridView, Vector, Matrix>& adapter, std::shared_ptr<Matrix> A_extended, std::shared_ptr<Matrix> A_ovlp_extended, std::shared_ptr<Vector> part_unity, const double eigenvalue_threshold,
                    int& nev, int nev_arpack = -1, const double shift = 0.001, const bool add_part_unity = false, const int verbose = 0) {

        if (nev_arpack == -1)
          nev_arpack = std::max(nev, 2);
        if (nev_arpack < nev)
          DUNE_THROW(Dune::Exception,"nev_arpack is less then nev!");

        // X * A_0 * X
        Matrix ovlp_mat(*A_ovlp_extended);
        using Dune::PDELab::Backend::native;
        for (auto row_iter = native(ovlp_mat).begin(); row_iter != native(ovlp_mat).end(); row_iter++) {
          for (auto col_iter = row_iter->begin(); col_iter != row_iter->end(); col_iter++) {
            *col_iter *= native(*part_unity)[row_iter.index()] * native(*part_unity)[col_iter.index()];
          }
        }

        // Setup Arpack for solving generalized eigenproblem
        std::cout << "ARPACK setup...";
        ArpackGeneo::ArPackPlusPlus_Algorithms<ISTLM, Vector> arpack(*A_extended);
        std::cout << " done" << std::endl;
        double eps = 0.0;

        std::vector<double> eigenvalues(nev_arpack,0.0);
        std::vector<Vector> eigenvectors(nev_arpack,Vector(adapter.getExtendedSize()));

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
            std::cout << "Process " << adapter.gridView().comm().rank() << " picked " << cnt << " eigenvectors" << std::endl;
          if (cnt == -1)
            DUNE_THROW(Dune::Exception,"No eigenvalue above threshold - not enough eigenvalues computed!");
        } else {
          cnt = nev;
        }

        // Write results to basis
        this->local_basis.resize(cnt);
        for (int base_id = 0; base_id < cnt; base_id++) {
          this->local_basis[base_id] = std::make_shared<Vector>(*part_unity);
          // scale partition of unity with eigenvector
          std::transform(
            this->local_basis[base_id]->begin(),this->local_basis[base_id]->end(),
            eigenvectors[base_id].begin(),
            this->local_basis[base_id]->begin(),
            std::multiplies<>()
            );
        }

        // Normalize basis vectors
        for (auto& v : this->local_basis) {
          *v *= 1.0 / v->two_norm2();
        }

        // Optionally add partition of unity to eigenvectors
        // Only if there is no near-zero eigenvalue (that usually already corresponds to a partition of unity!)
        if (add_part_unity && eigenvalues[0] > 1E-10) {
          this->local_basis.insert (this->local_basis.begin(), std::make_shared<Vector>(*part_unity));
          this->local_basis.pop_back();
        }
      }
    };


  }
}

#endif

#endif
