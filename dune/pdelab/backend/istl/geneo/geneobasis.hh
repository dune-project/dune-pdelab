#ifndef DUNE_PDELAB_BACKEND_ISTL_GENEO_GENEOBASIS_HH
#define DUNE_PDELAB_BACKEND_ISTL_GENEO_GENEOBASIS_HH

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
    template<class GFS, class M, class X, int dim>
    class GenEOBasis : public SubdomainBasis<X>
    {
      typedef Dune::PDELab::Backend::Native<M> ISTLM;
      typedef Dune::PDELab::Backend::Native<X> ISTLX;

    public:

      /*!
       * \brief Constructor.
       * \param gfs Grid function space.
       * \param AF_exterior Stiffness matrix with boundary conditions from problem definition and Neumann on processor boundaries.
       * \param AF_ovlp The same matrix as AF_exterior, but only assembled on overlap region (where more than 1 subdomain exists).
       *
       */
      GenEOBasis(const GFS& gfs, const M& AF_exterior, const M& AF_ovlp, const double eigenvalue_threshold, X& part_unity,
                int& nev, int nev_arpack = -1, double shift = 0.001, bool add_part_unity = false, int verbose = 0) {
        using Dune::PDELab::Backend::native;

        if (nev_arpack == -1)
          nev_arpack = std::max(nev, 2);
        if (nev_arpack < nev)
          DUNE_THROW(Dune::Exception,"eigenvectors_compute is less then eigenvectors or not specified!");

        // X * A_0 * X
        M ovlp_mat(AF_ovlp);
        for (auto row_iter = native(ovlp_mat).begin(); row_iter != native(ovlp_mat).end(); row_iter++) {
          for (auto col_iter = row_iter->begin(); col_iter != row_iter->end(); col_iter++) {
            *col_iter *= native(part_unity)[row_iter.index()] * native(part_unity)[col_iter.index()];
          }
        }

        // Setup Arpack for solving generalized eigenproblem
        ArpackGeneo::ArPackPlusPlus_Algorithms<ISTLM, ISTLX> arpack(native(AF_exterior));
        double eps = 0.0;

        std::vector<double> eigenvalues;
        std::vector<ISTLX> eigenvectors;
        eigenvectors.resize(nev_arpack);
        for (int i = 0; i < nev_arpack; i++) {
          eigenvectors[i] = native(X(gfs,0.0));
        }

        arpack.computeGenNonSymMinMagnitude(native(ovlp_mat), eps, nev_arpack, eigenvectors, eigenvalues, shift);

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
            std::cout << "Process " << gfs.gridView().comm().rank() << " picked " << cnt << " eigenvectors" << std::endl;
          if (cnt == -1)
            DUNE_THROW(Dune::Exception,"No eigenvalue above threshold - not enough eigenvalues computed!");
        } else {
          cnt = nev;
        }

        // Write results to basis
        this->local_basis.resize(cnt);
        for (int base_id = 0; base_id < cnt; base_id++) {
          this->local_basis[base_id] = std::make_shared<X>(part_unity);
          for (int it = 0; it < native(eigenvectors[base_id]).N(); it++) {
            for(int j = 0; j < dim; j++)
              native(*this->local_basis[base_id])[it][j] *= eigenvectors[base_id][it][j];
          }
        }

        // Normalize basis vectors
        for (auto& v : this->local_basis) {
          *v *= 1.0 / v->two_norm2();
        }

        // Optionally add partition of unity to eigenvectors
        // Only if there is no near-zero eigenvalue (that usually already corresponds to a partition of unity!)
        if (add_part_unity && eigenvalues[0] > 1E-10) {
          this->local_basis.insert (this->local_basis.begin(), std::make_shared<X>(part_unity));
          this->local_basis.pop_back();
        }
      }
    };


  }
}

#endif

#endif //DUNE_PDELAB_BACKEND_ISTL_GENEO_GENEOBASIS_HH
