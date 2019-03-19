#ifndef DUNE_PDELAB_BACKEND_ISTL_GENEO_RNDGENEOBASIS_HH
#define DUNE_PDELAB_BACKEND_ISTL_GENEO_RNDGENEOBASIS_HH

#include <algorithm>
#include <functional>

#include <dune/pdelab/backend/istl/geneo/subdomainbasis.hh>
#include <dune/pdelab/backend/istl/geneo/arpackpp_geneo.hh>
#include <dune/istl/cholmod.hh>
#include <random>

#include <dune/istl/matrixmarket.hh>

#if HAVE_SUITESPARSE_UMFPACK

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

    private:

      typedef Dune::PDELab::Backend::Native<M> ISTLM;
      typedef Dune::PDELab::Backend::Native<X> ISTLX;

      class IgnoreBlock
      {
      public:
        IgnoreBlock (bool val) : _val(val) {}
        bool operator[](std::size_t) const {return _val;}
      private:
        bool _val;
      };

      // dummy for empty ignore block vector
      //template <class M>
      class Ignore
      {
      public:
        Ignore (const M& mat)
         : _mat(mat), cnt(0) {
           using Dune::PDELab::Backend::native;
           for (int i = 0; i < native(_mat).N(); i++) { // TODO: instead go with constraints container, set constrained dofs to 1 via pdelab?
             if (native(_mat)[i][i] == 1.0)
               cnt ++;
           }
         }
        IgnoreBlock operator[](std::size_t i) const {
          using Dune::PDELab::Backend::native;
          return IgnoreBlock(native(_mat)[i][i] == 1.0);
        }
        std::size_t count() const { return cnt; }
      private:
        const M& _mat;
        int cnt;
      };

    public:

      /*!
       * \brief Constructor.
       * \param gfs Grid function space.
       * \param AF_exterior Stiffness matrix with boundary conditions from problem definition and Neumann on processor boundaries.
       * \param AF_ovlp The same matrix as AF_exterior, but only assembled on overlap region (where more than 1 subdomain exists).
       * \param eigenvalue_threshold Threshold up to which eigenvalue an eigenpair should be included in the basis. If negative, no thresholding.
       * \param part_unity Partition of unity to construct the basis with.
       * \param nev With thresholding, returns number of eigenvectors below threshold. Else, prescribes how many to use.
       * \param shift The shift to be used in ARPACK's shift invert solver mode. May need to be adjusted for extreme eigenvalue distributions.
       * \param add_part_unity Whether to explicitly add the partition of unity itself in the coarse basis.
       * \param verbose Verbosity value.
       */
      RndGenEOBasis(const GFS& gfs, const M& AF_exterior, const M& AF_ovlp, const double eigenvalue_threshold, X& part_unity,
                int nev, int verbose = 0) {
        using Dune::PDELab::Backend::native;

        Dune::Timer timer;
        Dune::Timer timer_rnd(false);

        // X * A_0 * X
        M ovlp_mat(AF_ovlp);
        for (auto row_iter = native(ovlp_mat).begin(); row_iter != native(ovlp_mat).end(); row_iter++) {
          for (auto col_iter = row_iter->begin(); col_iter != row_iter->end(); col_iter++) {
            *col_iter *= native(part_unity)[row_iter.index()] * native(part_unity)[col_iter.index()];
          }
        }
        MPI_Barrier (MPI_COMM_WORLD);
        if (verbose > 0) std::cout << "X * A0 * X: " << timer.elapsed() << std::endl; timer.reset();

        /*if (verbose > 0) {
          Dune::storeMatrixMarket(native(AF_ovlp), "AF_ovlp.mm");
          Dune::storeMatrixMarket(native(AF_exterior), "AF_exterior.mm");
          Dune::storeMatrixMarket(native(part_unity), "part_unity.mm");
        }*/

        Ignore ignore(AF_exterior);

        timer.reset();

        // Transposed cholesky factor inverse
        Cholmod<ISTLM> cholmod;
        //cholmod.setMatrix(native(AF_exterior), &ignore);
        {
          // TODO: Needed?
          // We do a shifted version to regularise the problem
          double shift = 0.01;
          M AF_ext_shift (AF_exterior);
          native(AF_ext_shift).axpy(-shift, native(AF_ovlp));
          cholmod.setMatrix(native(AF_ext_shift), &ignore);
        }
        if (verbose > 0) std::cout << "cholesky factor: " << timer.elapsed() << std::endl; timer.reset();
        ISTLM factor = cholmod.getFactor();
        if (verbose > 0) std::cout << "cholesky extract: " << timer.elapsed() << std::endl; timer.reset();
        Dune::UMFPack<ISTLM> factor_inverse(factor);
        if (verbose > 0) std::cout << "cholesky inverse: " << timer.elapsed() << std::endl; timer.reset();


        Dune::UMFPack<ISTLM> source_inverse(native(AF_exterior));
        MPI_Barrier (MPI_COMM_WORLD);
        if (verbose > 0) std::cout << "source_inverse: " << timer.elapsed() << std::endl; timer.reset();

        std::default_random_engine generator;
        std::normal_distribution<double> distribution(0.0,1.0);

        Dune::InverseOperatorResult result;

        this->local_basis.resize(nev);
        for (int i = 0; i < nev; i++) {

          ISTLX ext(AF_exterior.N());


          ISTLX rnd (factor.N());
          for (int j = 0; j < rnd.N(); j++)
            rnd[j] = distribution(generator);

          ISTLX restr(factor.N());
          factor_inverse.apply(restr, rnd, result);

          // Extend to full domain size including Dirichlet DOFs
          auto restr_iter = restr.begin();
          for (int j = 0; j < ext.N(); j++) {
            if (ignore[j][0]) {
              ext[j] = 0.0;
            } else {
              ext[j] = *restr_iter;
            }
            restr_iter++;
          }

          std::transform(
            ext.begin(),ext.end(),
                         part_unity.begin(),
                         ext.begin(),
                         std::multiplies<>()
          );

          for (int reiter = 0; reiter < 1; reiter++) {

            ISTLX temp(AF_exterior.N());
            native(ovlp_mat).mv(ext, temp);
            ext = temp;

            std::transform(
              ext.begin(),ext.end(),
                           part_unity.begin(),
                           ext.begin(),
                           std::multiplies<>()
            );

            source_inverse.apply(temp, ext, result);
            ext = temp;

            std::transform(
              ext.begin(),ext.end(),
                           part_unity.begin(),
                           ext.begin(),
                           std::multiplies<>()
            );

          }

          this->local_basis[i] = std::make_shared<X>(gfs, 0.0);
          std::copy(
            ext.begin(),ext.end(),
                         native(*this->local_basis[i]).begin()
          );
        }
        MPI_Barrier (MPI_COMM_WORLD);
        if (verbose > 0) std::cout << "rnd basis: " << timer.elapsed() << std::endl; timer.reset();


        // Gram-Schmidt

        std::vector<double> squarednorms(0);
        X scaled(gfs, 0.0);
        for (int i = 0; i < nev; i++) { // TODO: Orthonormalisierung bzgl range product?
          for (int j = 0; j < i; j++) {

            // ovlp_mat als sk.prod?
            double f = dotProduct(AF_exterior, *this->local_basis[i], *this->local_basis[j], scaled) / squarednorms[j];//(*this->local_basis[j] * *this->local_basis[j]);
            scaled = *this->local_basis[j]; // axpy instead?
            scaled *= f;
            *this->local_basis[i] -= scaled;
            //if (verbose > 1)
            //  std::cout << "scale " << f << "\t to " << *this->local_basis[i] * *this->local_basis[j] << std::endl;
          }
          squarednorms.push_back(dotProduct(AF_exterior, *this->local_basis[i], *this->local_basis[i], scaled));
        }

        /*X scaled(gfs, 0.0);
        for (int i = 0; i < nev; i++) {
          *this->local_basis[i] *= 1.0 / std::sqrt(dotProduct(AF_exterior, *this->local_basis[i], *this->local_basis[i], scaled));
        }*/

        MPI_Barrier (MPI_COMM_WORLD);
        if (verbose > 0) std::cout << "random vector gen: " << timer_rnd.elapsed() << std::endl;
        if (verbose > 0) std::cout << "Gram-Schmidt: " << timer.elapsed() << std::endl; timer.reset();
      }

    private:
      double matDotProduct (const M& mat, const X& x, const X& y, X& temp) {
        using Dune::PDELab::Backend::native;
        native(mat).mv(native(y), native(temp));
        return x * temp;
      }
      double dotProduct (const M& mat, const X& x, const X& y, X& temp) {
        //using Dune::PDELab::Backend::native;
        //native(mat).mv(native(y), native(temp));
        return x * y;
      }
    };


  }
}

#endif

#endif //DUNE_PDELAB_BACKEND_ISTL_GENEO_RNDGENEOBASIS_HH
