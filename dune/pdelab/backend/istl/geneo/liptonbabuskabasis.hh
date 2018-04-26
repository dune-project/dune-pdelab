#ifndef DUNE_PDELAB_BACKEND_ISTL_GENEO_LIPTONBABUSKABASIS_HH
#define DUNE_PDELAB_BACKEND_ISTL_GENEO_LIPTONBABUSKABASIS_HH

#include <dune/pdelab/backend/istl/geneo/subdomainbasis.hh>
#include <dune/pdelab/backend/istl/geneo/arpackpp_geneo.hh>

#if HAVE_ARPACKPP

namespace Dune {
  namespace PDELab {

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

        ArpackGeneo::ArPackPlusPlus_Algorithms<ISTLM, ISTLX> arpack(native(AF_exterior));
        double eps = .0001;

        std::vector<double> eigenvalues;
        std::vector<ISTLX> eigenvectors;
        eigenvectors.resize(nev_arpack);
        for (int i = 0; i < nev_arpack; i++) {
          eigenvectors[i] = native(X(gfs,0.0));
        }

        arpack.computeGenNonSymMinMagnitude(native(AF_interior), eps, nev_arpack, eigenvectors, eigenvalues, shift);


        /*for(int i = 0; i < nev_arpack; i++){
            auto check = eigenvectors[i];
            native(AF_exterior).mv(eigenvectors[i],check);
            auto check2 = eigenvectors[i];
            native(ovlp_mat).mv(eigenvectors[i],check2);
            check2 *= eigenvalues[i];
            check -= check2;
            if(native(check).infinity_norm() > 1e-6) std::cout << "Rank " << gfs.gridView().comm().rank() << " Error in EV calculation " << native(check).infinity_norm() << std::endl;
        }*/


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
          for (int it = 0; it < native(eigenvectors[base_id]).N(); it++) {
            for(int j = 0; j < dim; j++)
              native(*this->local_basis[base_id])[it][j] *= eigenvectors[base_id][it][j];
          }
        }

        // Normalize basis vectors
        for (int i = 0; i < this->local_basis.size(); i++) {
          native(*(this->local_basis[i])) *= 1.0 / (native(*(this->local_basis[i])) * native(*(this->local_basis[i])));
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
