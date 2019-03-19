#ifndef DUNE_PDELAB_BACKEND_ISTL_GENEO_FASTRNDGENEOBASIS_HH
#define DUNE_PDELAB_BACKEND_ISTL_GENEO_FASTRNDGENEOBASIS_HH

#include <algorithm>
#include <functional>

#include <dune/pdelab/backend/istl/geneo/subdomainbasis.hh>
#include <dune/pdelab/backend/istl/geneo/arpackpp_geneo.hh>
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
    class FastRndGenEOBasis : public SubdomainBasis<X>
    {

    private:

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
       * \param shift The shift to be used in ARPACK's shift invert solver mode. May need to be adjusted for extreme eigenvalue distributions.
       * \param add_part_unity Whether to explicitly add the partition of unity itself in the coarse basis.
       * \param verbose Verbosity value.
       */
      FastRndGenEOBasis(const GFS& gfs, const M& AF_exterior, const M& AF_ovlp, const double eigenvalue_threshold, X& part_unity,
                int nev, int verbose = 0, int reiterations = 1) {
        using Dune::PDELab::Backend::native;

        Dune::Timer timer;
        Dune::Timer timer_rnd(false);

        /*if (nev_arpack == -1)
          nev_arpack = std::max(nev, 2);
        if (nev_arpack < nev)
          DUNE_THROW(Dune::Exception,"nev_arpack is less then nev!");*/

        std::default_random_engine generator;
        std::normal_distribution<double> distribution;

        int testspace_size = nev * 2;
        std::vector<std::shared_ptr<X> > testspace(testspace_size);
        for (int i = 0; i < testspace_size; i++) {
          testspace[i] = std::make_shared<X>(gfs, 0.0);
          X& testvec = *(testspace[i]);

          for (int j = 0; j < AF_exterior.N(); j++) {
            for(int j_block = 0; j_block < ISTLM::block_type::rows; j_block++){
              if (native(AF_exterior)[j][j][j_block][j_block] != 1.0)
                native(testvec)[j][j_block] = distribution(generator) / std::sqrt(native(AF_exterior)[j][j][j_block][j_block]);
              else
                native(testvec)[j][j_block] = 0.0;
            }
          }
          for (int j = 0; j < AF_exterior.N(); j++) {
            for(int j_block = 0; j_block < ISTLM::block_type::rows; j_block++){
              native(testvec)[j][j_block] *= native(part_unity)[j][j_block];
            }
          }
          //if (verbose > 0) std::cout << "testevec init norm: " << (testspace[i])->two_norm() << std::endl;
        }

        {
          X proj(gfs, 0.0);
          double maxnorm = 0.0;
          for (int ti = 0; ti < testspace_size; ti++) {
            X& testvec = *(testspace[ti]);
            /*native(AF_ovlp).mv(native(testvec), native(proj));
             *              double factor = proj * testvec;*/
            double factor = matDotProduct(AF_ovlp, testvec, testvec, proj);
            maxnorm = std::max(maxnorm, factor);
          }
          maxnorm = std::sqrt(maxnorm);
          if (verbose > 0) std::cout << "maxnorm: " << maxnorm << std::endl;
        }

        // X * A_0 * X
        M ovlp_mat(AF_ovlp);
        for (auto row_iter = native(ovlp_mat).begin(); row_iter != native(ovlp_mat).end(); row_iter++) {
          for (auto col_iter = row_iter->begin(); col_iter != row_iter->end(); col_iter++) {
            *col_iter *= native(part_unity)[row_iter.index()] * native(part_unity)[col_iter.index()];
          }
        }
        MPI_Barrier (MPI_COMM_WORLD);
        if (verbose > 0) std::cout << "XA0X: " << timer.elapsed() << std::endl; timer.reset();

        /*if (verbose > 0) {
          Dune::storeMatrixMarket(native(AF_ovlp), "AF_ovlp.mm");
          Dune::storeMatrixMarket(native(AF_exterior), "AF_exterior.mm");
          Dune::storeMatrixMarket(native(part_unity), "part_unity.mm");
        }*/

        timer.reset();


        Dune::UMFPack<ISTLM> source_inverse(native(AF_exterior));
        /*Dune::MatrixAdapter<ISTLM, ISTLX, ISTLX> matop(native(AF_exterior));
        Dune::SeqGS<ISTLM, ISTLX, ISTLX> seqgs(native(AF_exterior), 10, 1.0);
        Dune::CGSolver<ISTLX> cgsolver(matop, seqgs, 1e-6, 1000, 1);*/

        MPI_Barrier (MPI_COMM_WORLD);
        if (verbose > 0) std::cout << "source_inverse: " << timer.elapsed() << std::endl; timer.reset();


        Dune::InverseOperatorResult result;


        this->local_basis.resize(nev);
        for (int i = 0; i < nev; i++) {

          ISTLX ext(AF_exterior.N());

          timer_rnd.start();

          // Apply D^-1/2
          for (int j = 0; j < ext.N(); j++) {
            for(int j_block = 0; j_block < ISTLM::block_type::rows; j_block++){
              if (native(AF_exterior)[j][j][j_block][j_block] != 1.0)
                ext[j][j_block] = distribution(generator) / std::sqrt(native(AF_exterior)[j][j][j_block][j_block]);
              else
                ext[j][j_block] = 0.0;
            }
          }
          timer_rnd.stop();


          for (int j = 0; j < ext.N(); j++) {
            for(int j_block = 0; j_block < ISTLM::block_type::rows; j_block++){
              ext[j][j_block] *= native(part_unity)[j][j_block];
            }
          }
          /*std::transform(
            ext.begin(),ext.end(),
                         part_unity.begin(),
                         ext.begin(),
                         std::multiplies<>()
          );*/

          for (int reiter = 0; reiter < reiterations; reiter++) {

            ISTLX temp(AF_exterior.N());
            native(ovlp_mat).mv(ext, temp);
            ext = temp;

            for (int j = 0; j < ext.N(); j++) {
              for(int j_block = 0; j_block < ISTLM::block_type::rows; j_block++){
                ext[j][j_block] *= native(part_unity)[j][j_block];
              }
            }
            /*std::transform(
              ext.begin(),ext.end(),
                           part_unity.begin(),
                           ext.begin(),
                           std::multiplies<>()
            );*/

            //cgsolver.apply(temp, ext, result);
            source_inverse.apply(temp, ext, result);
            ext = temp;

            for (int j = 0; j < ext.N(); j++) {
              for(int j_block = 0; j_block < ISTLM::block_type::rows; j_block++){
                ext[j][j_block] *= native(part_unity)[j][j_block];
              }
            }
            /*std::transform(
              ext.begin(),ext.end(),
                           part_unity.begin(),
                           ext.begin(),
                           std::multiplies<>()
            );*/

          }

          this->local_basis[i] = std::make_shared<X>(gfs, 0.0);
          std::copy(
            ext.begin(),ext.end(),
                         native(*this->local_basis[i]).begin()
          );


          // Gram-Schmidt

          std::vector<double> squarednorms(0);
          X scaled(gfs, 0.0);
          for (int gi = 0; gi <= i; gi++) { // TODO: Orthonormalisierung bzgl range product?
            for (int gj = 0; gj < gi; gj++) {

              // ovlp_mat als sk.prod?
              double f = matDotProduct(AF_ovlp, *this->local_basis[gi], *this->local_basis[gj], scaled) / squarednorms[gj];//(*this->local_basis[j] * *this->local_basis[j]);
              scaled = *this->local_basis[gj]; // axpy instead?
              scaled *= f;
              *this->local_basis[gi] -= scaled;
              //if (verbose > 1)
              //  std::cout << "scale " << f << "\t to " << *this->local_basis[i] * *this->local_basis[j] << std::endl;
            }
            squarednorms.push_back(matDotProduct(AF_ovlp, *this->local_basis[gi], *this->local_basis[gi], scaled));
          }


          // Project out of test space

          X proj(gfs, 0.0);
          for (int bi = i; bi <= i; bi++) {
            for (int ti = 0; ti < testspace_size; ti++) {
              X& testvec = *(testspace[ti]);
              //if (verbose > 0) std::cout << "testevec new norm: " << (testspace[i])->two_norm() << std::endl;

              /*native(AF_ovlp).mv(native(testvec), native(proj));
              double factor = native(proj).dot(native(*this->local_basis[bi]));*/
              double factor = matDotProduct(AF_ovlp, *this->local_basis[bi], testvec, proj)
                            / matDotProduct(AF_ovlp, *this->local_basis[bi], *this->local_basis[bi], proj);
              //if (verbose > 0) std::cout << native(*(testspace[ti])).two_norm();
              testvec.axpy(-factor, *this->local_basis[bi]);
              //if (verbose > 0) std::cout << "\t" << native(testvec).two_norm() << "\t by factor " << -factor << std::endl;
            }
          }
          {
            double maxnorm = 0.0;
            for (int ti = 0; ti < testspace_size; ti++) {
              X& testvec = *(testspace[ti]);
              /*native(AF_ovlp).mv(native(testvec), native(proj));
              double factor = proj * testvec;*/
              double factor = matDotProduct(AF_ovlp, testvec, testvec, proj);
              maxnorm = std::max(maxnorm, factor);
            }
            maxnorm = std::sqrt(maxnorm);
            if (verbose > 0) std::cout << "maxnorm: " << maxnorm << std::endl;
          }

        }
        MPI_Barrier (MPI_COMM_WORLD);
        if (verbose > 0) std::cout << "rnd basis: " << timer.elapsed() << std::endl; timer.reset();


        // Gram-Schmidt

        std::vector<double> squarednorms(0);
        X scaled(gfs, 0.0);
        for (int i = 0; i < nev; i++) { // TODO: Orthonormalisierung bzgl range product?
          for (int j = 0; j < i; j++) {

            // ovlp_mat als sk.prod?
            double f = matDotProduct(AF_ovlp, *this->local_basis[i], *this->local_basis[j], scaled) / squarednorms[j];//(*this->local_basis[j] * *this->local_basis[j]);
            scaled = *this->local_basis[j]; // axpy instead?
            scaled *= f;
            *this->local_basis[i] -= scaled;
            //if (verbose > 1)
            //  std::cout << "scale " << f << "\t to " << *this->local_basis[i] * *this->local_basis[j] << std::endl;
          }
          squarednorms.push_back(matDotProduct(AF_ovlp, *this->local_basis[i], *this->local_basis[i], scaled));
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

#endif //DUNE_PDELAB_BACKEND_ISTL_GENEO_FASTRNDGENEOBASIS_HH
