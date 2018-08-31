#ifndef DUNE_PDELAB_BACKEND_ISTL_GENEO_RNDGENEOBASIS_HH
#define DUNE_PDELAB_BACKEND_ISTL_GENEO_RNDGENEOBASIS_HH

#include <algorithm>
#include <functional>

#include <dune/pdelab/backend/istl/geneo/subdomainbasis.hh>
#include <dune/pdelab/backend/istl/geneo/arpackpp_geneo.hh>
#include <dune/istl/matrixmarket.hh>

#include <stdio.h>
#include <stdlib.h>

namespace Dune {
  namespace PDELab {

    /*!
     * \brief Implementation of the GenEO coarse basis
     * See Spillane et al., 2014, 'Abstract robust coarse spaces for systems of PDEs via generalized eigenproblems in the overlaps'
     * This coarse space is based on generalized eigenpoblems defined on the full stiffness matrix of a subdomain and one assembled only
     * on the area where this subdomain overlaps with others.
     */
    template<class GFS, class M, class X, int dim>
    class RndGenEOBasis : public SubdomainBasis<X>
    {
      typedef Dune::PDELab::Backend::Native<M> ISTLM;
      typedef Dune::PDELab::Backend::Native<X> ISTLX;

    public:

      /*!
       * \brief Constructor.
       * \param gfs Grid function space.
       * \param AF_exterior Stiffness matrix with boundary conditions from problem definition and Neumann on processor boundaries.
       * \param AF_ovlp The same matrix as AF_exterior, but only assembled on overlap region (where more than 1 subdomain exists).
       * \param eigenvalue_threshold Threshold up to which eigenvalue an eigenpair should be included in the basis. If negative, no thresholding.
       * \param part_unity Partition of unity to construct the basis with.
       * \param nev With thresholding, returns number of eigenvectors below threshold. Else, prescribes how many to use.
       * \param nev_arpack How many eigenpairs ARPACK is supposed to compute. Larger numbers may increase its stability. -1 for default.
       * \param shift The shift to be used in ARPACK's shift invert solver mode. May need to be adjusted for extreme eigenvalue distributions.
       * \param add_part_unity Whether to explicitly add the partition of unity itself in the coarse basis.
       * \param verbose Verbosity value.
       */
      RndGenEOBasis(const GFS& gfs, const M& AF_exterior, const M& AF_ovlp, const double eigenvalue_threshold, X& part_unity,
                int& nev, int nev_arpack = -1, double shift = 0.001, bool add_part_unity = false, int verbose = 0) {
        using Dune::PDELab::Backend::native;

        if (nev_arpack == -1)
          nev_arpack = std::max(nev, 2);
        if (nev_arpack < nev)
          DUNE_THROW(Dune::Exception,"nev_arpack is less then nev!");

        // Defer solving PDE to python script via MatrixMarket
        std::string suffix = std::to_string(gfs.gridView().comm().rank()) + ".mm";

        std::string mm_AF_exterior = "AF_exterior" + suffix;
        std::string mm_AF_ovlp = "AF_ovlp" + suffix;
        std::string mm_part_unity = "part_unity" + suffix;
        std::string mm_basis = "basis" + suffix;
        std::string mm_dims = "dimensions" + std::to_string(gfs.gridView().comm().rank());
        Dune::storeMatrixMarket (native(AF_exterior), mm_AF_exterior);
        Dune::storeMatrixMarket (native(AF_ovlp), "AF_ovlp" + suffix);
        Dune::storeMatrixMarket (native(part_unity), "part_unity" + suffix);

        system(("python gen_basis.py " + mm_part_unity + " " + mm_AF_ovlp + " " + mm_AF_exterior + " " + mm_basis + " " + mm_dims + " rnd").c_str());

        // Read number of resulting basis vectors etc.
        std::ifstream input(mm_dims);
        int dofs, cnt;
        input >> dofs >> cnt;

        // Import results
        this->local_basis.resize(cnt);
        X input_vector (gfs, 0.0);
        for (int base_id = 0; base_id < cnt; base_id++) {
          Dune::loadMatrixMarket (native(input_vector), std::to_string(base_id) + "_" + mm_basis);
          this->local_basis[base_id] = std::make_shared<X>(input_vector);
        }

        // Normalize basis vectors
        for (auto& v : this->local_basis) {
          *v *= 1.0 / v->two_norm2();
        }

      }
    };


  }
}

#endif //DUNE_PDELAB_BACKEND_ISTL_GENEO_RNDGENEOBASIS_HH
