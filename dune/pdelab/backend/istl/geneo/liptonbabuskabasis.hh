#ifndef DUNE_PDELAB_BACKEND_ISTL_GENEO_LIPTONBABUSKABASIS_HH
#define DUNE_PDELAB_BACKEND_ISTL_GENEO_LIPTONBABUSKABASIS_HH

#include <algorithm>
#include <functional>

#include <dune/pdelab/backend/istl/geneo/subdomainbasis.hh>
#include <dune/pdelab/backend/istl/geneo/arpackpp_geneo.hh>

#if HAVE_ARPACKPP

namespace Dune {
  namespace PDELab {

    /*!
     * \brief Experimental implementation of a basis according to Babuska and Lipton
     * 'Optimal Local Approximation Spaces for Generalized Finite Elmente Methods with Application to Multiscale Problems' (2011)
     */
    template<class GFS, class M, class X, class Y, int dim>
    class LiptonBabuskaBasis : public SubdomainBasis<X>
    {
      typedef Dune::PDELab::Backend::Native<M> ISTLM;
      typedef Dune::PDELab::Backend::Native<X> ISTLX;

    public:
      LiptonBabuskaBasis(const GFS& gfs, const M& AF_exterior, const M& AF_ovlp, const double eigenvalue_threshold, X& part_unity,
              int& nev, int nev_arpack, double shift = 0.001, bool add_part_unity = false, int verbose = 0) {
        using Dune::PDELab::Backend::native;

        if (nev_arpack == -1)
          nev_arpack = std::max(nev, 2);
        if (nev_arpack < nev)
          DUNE_THROW(Dune::Exception,"eigenvectors_compute is less then eigenvectors or not specified!");

        auto AF_interior = AF_exterior;
        native(AF_interior) -= native(AF_ovlp);

        ArpackGeneo::ArPackPlusPlus_Algorithms<ISTLM, X> arpack(native(AF_exterior));
        double eps = .0001;

        std::vector<double> eigenvalues(nev_arpack,0.0);
        std::vector<X> eigenvectors(nev_arpack,X(gfs,0.0));

        arpack.computeGenNonSymMinMagnitude(native(AF_interior), eps, eigenvectors, eigenvalues, shift);

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

        this->local_basis.resize(cnt);

        for (int base_id = 0; base_id < cnt; base_id++) {
          this->local_basis[base_id] = std::make_shared<X>(part_unity);
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

        if (add_part_unity && eigenvalues[0] > 1E-10) {
          this->local_basis.insert (this->local_basis.begin(), std::make_shared<X>(part_unity));
          this->local_basis.pop_back();
        }
      }

    };

  }
}

#endif

#endif //DUNE_PDELAB_BACKEND_ISTL_GENEO_LIPTONBABUSKABASIS_HH
