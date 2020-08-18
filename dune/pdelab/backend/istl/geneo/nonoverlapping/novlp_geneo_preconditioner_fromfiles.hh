#ifndef DUNE_PDELAB_BACKEND_ISTL_GENEO_NOVLP_GENEO_PRECONDITIONER_FROMFILES_HH
#define DUNE_PDELAB_BACKEND_ISTL_GENEO_NOVLP_GENEO_PRECONDITIONER_FROMFILES_HH

#include <dune/pdelab/backend/istl/geneo/nonoverlapping/geneobasis.hh>
#include <dune/pdelab/backend/istl/geneo/nonoverlapping/geneobasisfromfiles.hh>

namespace Dune {
  namespace PDELab {

    template <typename GO, typename Matrix, typename Vector>
    class NonoverlappingGenEOPreconditionerFromFiles : public Dune::Preconditioner<Vector,Vector> {

      using V = typename GO::Domain;
      using M = typename GO::Jacobian;
      using GFS = typename GO::Traits::TrialGridFunctionSpace;
      using GV = typename GFS::Traits::GridView;


    public:

      // define the category
      virtual Dune::SolverCategory::Category category() const override
      {
        return prec->category();
      }


      /*!
      *          \brief Prepare the preconditioner.
      *
      *          \copydoc Preconditioner::pre(X&,Y&)
      */
      virtual void pre (Vector& x, Vector& b) override
      {
        prec->pre(x, b);
      }

      /*!
      *          \brief Apply the precondioner.
      *
      *          \copydoc Preconditioner::apply(X&,const Y&)
      */
      virtual void apply (Vector& v, const Vector& d) override
      {
        prec->apply(v, d);
      }
      /*!
      *          \brief Clean up.
      *
      *          \copydoc Preconditioner::post(X&)
      */
      virtual void post (Vector& x) override
      {
        prec->post(x);
      }

      NonoverlappingGenEOPreconditionerFromFiles(const GO& go, const typename GO::Jacobian& A, int algebraic_overlap, int avg_nonzeros, const double eigenvalue_threshold, int& nev,
                          int nev_arpack = -1, const double shift = 0.001, int verbose = 0, int multiscale = 0, std::vector<int> proc_to_be_solved = {0}, std::string path_to_storage = "")
      {

        using Dune::PDELab::Backend::native;

        const GV& gv = go.trialGridFunctionSpace().gridView();

        Dune::NonoverlappingOverlapAdapter<GV, Vector, Matrix> adapter(gv, native(A), avg_nonzeros, algebraic_overlap);

        std::shared_ptr<Matrix> A_extended = adapter.extendMatrix(native(A));

        std::shared_ptr<Vector> part_unity = Dune::makePartitionOfUnity<GV, Matrix, Vector>(adapter, *A_extended);


        std::string basename = path_to_storage+"EV";

        if (multiscale==1) { // first step of the multiscale FRAMEWORK: solving the pristine model & saving it to file
          subdomainbasis = std::make_shared<Dune::PDELab::NonoverlappingGenEOBasis<GO, Matrix, Vector>>(go, adapter, A, A_extended, part_unity, eigenvalue_threshold, nev, nev_arpack, shift, false, 2);
          // Save the EV basis
          int rank = adapter.gridView().comm().rank();
          subdomainbasis->to_file(basename, rank);
        } else if (multiscale==2) { // other step of the multiscale FRAMEWORK: loading the subdomain basis from files
          std::vector<int>::iterator it = std::find(std::begin(proc_to_be_solved), std::end(proc_to_be_solved), adapter.gridView().comm().rank());
          if (it != proc_to_be_solved.end()) {
            subdomainbasis = std::make_shared<Dune::PDELab::NonoverlappingGenEOBasis<GO, Matrix, Vector>>(go, adapter, A, A_extended, part_unity, eigenvalue_threshold, nev, nev_arpack, shift, false, 2);
          } else {
            subdomainbasis = std::make_shared<Dune::PDELab::NonoverlappingGenEOBasisFromFiles<GV, Matrix, Vector>>(adapter, basename);
          }
        } else if (multiscale==3) {
           // Test case for write/read database
           // Check numbers of digit initially and after the saving/reading procedure
           subdomainbasis = std::make_shared<Dune::PDELab::NonoverlappingGenEOBasis<GO, Matrix, Vector>>(go, adapter, A, A_extended, part_unity, eigenvalue_threshold, nev, nev_arpack, shift, false, 2);

           int rank = adapter.gridView().comm().rank();
           subdomainbasis->to_file(basename, rank);

           //auto fromfile_subdomainbasis;
           auto fromfile_subdomainbasis = std::make_shared<Dune::PDELab::NonoverlappingGenEOBasisFromFiles<GV, Matrix, Vector>>(adapter, basename);
           fromfile_subdomainbasis->to_file(basename+"_rewritten", rank);
           // auto tmp = (*subdomainbasis->get_basis_vector(0));
           // tmp -= (*fromfile_subdomainbasis->get_basis_vector(0));
           // std::cout << tmp.two_norm() << std::endl;

        } else { // Classic case: no need to use multiscale FRAMEWORK
          // subdomainbasis = std::make_shared<Dune::PDELab::NonoverlappingGenEOBasis<GV, Matrix, Vector>>(adapter, *A_extended, *A_ovlp_extended, *part_unity, eigenvalue_threshold, nev, nev_arpack, shift,false,2);

          subdomainbasis = std::make_shared<Dune::PDELab::NonoverlappingGenEOBasis<GO, Matrix, Vector>>(go, adapter, A, A_extended, part_unity, eigenvalue_threshold, nev, nev_arpack, shift, false, 2);
        }

        auto coarse_space = std::make_shared<Dune::PDELab::NewSubdomainProjectedCoarseSpace<GV, Matrix, Vector>>(adapter, gv, *A_extended, subdomainbasis, verbose);

        prec = std::make_shared<Dune::PDELab::ISTL::NonoverlappingTwoLevelOverlappingAdditiveSchwarz<GV, Matrix, Vector>>(adapter, A_extended, *part_unity, coarse_space, true, verbose);
      }

      std::shared_ptr<Dune::PDELab::SubdomainBasis<Vector>> getSubdomainbasis() {
        return subdomainbasis;
      }

    private:

      std::shared_ptr<Dune::PDELab::ISTL::NonoverlappingTwoLevelOverlappingAdditiveSchwarz<GV, Matrix, Vector>> prec = nullptr;

      std::shared_ptr<Dune::PDELab::SubdomainBasis<Vector>> subdomainbasis;

      // M A_ovlp;
      // std::shared_ptr<Matrix> A_ovlp_extended;
      // std::shared_ptr<Matrix> A_extended;
      // std::shared_ptr<Vector> part_unity;

      // std::shared_ptr<Dune::PDELab::ISTL::NonoverlappingTwoLevelOverlappingAdditiveSchwarz<typename GO::Traits::TrialGridFunctionSpace::Traits::GridView, Matrix, Vector>> prec = nullptr;

    };

  }
}

#endif
