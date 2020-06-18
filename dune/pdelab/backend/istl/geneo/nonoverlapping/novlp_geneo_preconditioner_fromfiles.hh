#ifndef DUNE_PDELAB_BACKEND_ISTL_GENEO_NOVLP_GENEO_PRECONDITIONER_FROMFILES_HH
#define DUNE_PDELAB_BACKEND_ISTL_GENEO_NOVLP_GENEO_PRECONDITIONER_FROMFILES_HH

#include <dune/pdelab/backend/istl/geneo/nonoverlapping/geneobasis.hh>
#include <dune/pdelab/backend/istl/geneo/nonoverlapping/geneobasisfromfiles.hh>

namespace Dune {
  namespace PDELab {

    template <typename GO, typename ScalarMatrix, typename Matrix, typename ScalarVector, typename Vector>
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
                          int nev_arpack = -1, const double shift = 0.001, int verbose = 0, int multiscale = 0, std::vector<int> proc_to_be_solved = {0})
      : A_ovlp(go){

        V x(go.trialGridFunctionSpace(),0.0); // NOTE: We assume linear problems, so simply set x to zero here!

        const GFS& gfs = go.trialGridFunctionSpace();
        const GV& gv = gfs.gridView();


        using Dune::PDELab::Backend::native;


        Dune::NonoverlappingOverlapAdapter<GV, Vector, Matrix> adapter(gv, native(A), avg_nonzeros, algebraic_overlap);

        A_extended = adapter.extendMatrix(native(A));
        part_unity = Dune::makePartitionOfUnity<GV, Matrix, Vector>(adapter, *A_extended);


        using Attribute = Dune::EPISAttribute;
        Dune::AllSet<Attribute> allAttribute;
        auto allinterface = std::shared_ptr<Dune::Interface>(new Dune::Interface());
        allinterface->build(*adapter.getRemoteIndices(),allAttribute,allAttribute); // all to all communication
        auto communicator = std::shared_ptr<DuneWithRank::BufferedCommunicator>(new DuneWithRank::BufferedCommunicator());
        communicator->build<Vector>(*allinterface);

        Dune::PDELab::MultiVectorBundle<GV, Vector, Matrix> remotePartUnities(adapter);
        remotePartUnities.localVector_ = part_unity;
        communicator->forward<Dune::PDELab::MultiGatherScatter<Dune::PDELab::MultiVectorBundle<GV, Vector, Matrix>>>(remotePartUnities,remotePartUnities); // make function known in other subdomains


        // Assemble fine grid matrix defined only on overlap region
        auto part_unity_restricted = std::make_shared<Vector>(native(A).N());
        adapter.restrictVector(*part_unity, *part_unity_restricted);



        auto es_pou_excluder = std::make_shared<EntitySetPartUnityExcluder<Vector, GV, GFS>> (gfs, part_unity_restricted);
        gfs.entitySet().setExcluder(es_pou_excluder);



        go.jacobian(x,A_ovlp);



        M newmat(go);
        // Provide neighbors with matrices assembled exclusively on respective overlap area
        auto extended_matrices = adapter.lambdaMultiExtendMatrix(native(A_ovlp), native(A), [&](int i){
          std::shared_ptr<Vector> neighbor_part_unity = remotePartUnities.getVectorForRank(i);

          adapter.restrictVector(*neighbor_part_unity, *part_unity_restricted);

          std::shared_ptr<EntitySetExcluder<Vector, GV>> es_local_pou_excluder = std::make_shared<EntitySetPartUnityExcluder<Vector, GV, GFS>> (gfs, part_unity_restricted);
          gfs.entitySet().setExcluder(es_local_pou_excluder);

          for (auto rIt=native(newmat).begin(); rIt!=native(newmat).end(); ++rIt)
            for (auto cIt=rIt->begin(); cIt!=rIt->end(); ++cIt) {
              *cIt = 0.0;
            }

          go.jacobian(x,newmat);

          return stackobject_to_shared_ptr(native(newmat));
        });

        A_ovlp_extended = extended_matrices.first;
        A_extended = extended_matrices.second;


        // Enforce problem's Dirichlet condition on PoU
        const int block_size = Vector::block_type::dimension;
        for (auto rIt=A_extended->begin(); rIt!=A_extended->end(); ++rIt) {
          for(int block_i = 0; block_i < block_size; block_i++){ //loop over block, TODO block_size
            bool isDirichlet = true;
            for (auto cIt=rIt->begin(); cIt!=rIt->end(); ++cIt)
            {
              for(int block_j = 0; block_j<block_size; block_j++){
                if ((rIt.index() != cIt.index() || block_i!=block_j) && (*cIt)[block_i][block_j] != 0.0) {
                  isDirichlet = false;
                  break;
                }
              }
              if(!isDirichlet) break;
            }
            if (isDirichlet) {
              (*part_unity)[rIt.index()][block_i] = .0;
            }
          }
        }

        std::string basename = "Offline/EV";

        if (multiscale==0) { // Classic case: no need to use multiscale FRAMEWORK
          subdomainbasis = std::make_shared<Dune::PDELab::NonoverlappingGenEOBasis<GV, Matrix, Vector>>(adapter, *A_extended, *A_ovlp_extended, *part_unity, eigenvalue_threshold, nev, nev_arpack, shift,false,2);
        } else if (multiscale==1) { // first step of the multiscale FRAMEWORK: solving the pristine model & saving it to file
          subdomainbasis = std::make_shared<Dune::PDELab::NonoverlappingGenEOBasis<GV, Matrix, Vector>>(adapter, *A_extended, *A_ovlp_extended, *part_unity, eigenvalue_threshold, nev, nev_arpack, shift,false,2);
          // Save the EV basis
          int rank = adapter.gridView().comm().rank();
          subdomainbasis->to_file(basename, rank);
        } else if (multiscale==2) { // other step of the multiscale FRAMEWORK: loading the subdomain basis from files
          std::vector<int>::iterator it = std::find(std::begin(proc_to_be_solved), std::end(proc_to_be_solved), adapter.gridView().comm().rank());
          if (it != proc_to_be_solved.end()) {
            subdomainbasis = std::make_shared<Dune::PDELab::NonoverlappingGenEOBasis<GV, Matrix, Vector>>(adapter, *A_extended, *A_ovlp_extended, *part_unity, eigenvalue_threshold, nev, nev_arpack, shift,false,2);
          } else {
            subdomainbasis = std::make_shared<Dune::PDELab::NonoverlappingGenEOBasisFromFiles<GV, Matrix, Vector>>(adapter, basename);
          }
        }

        else if (multiscale==3) {
           // Test case for write/read database
           // Check numbers of digit initially and after the saving/reading procedure
           subdomainbasis = std::make_shared<Dune::PDELab::NonoverlappingGenEOBasis<GV, Matrix, Vector>>(adapter, *A_extended, *A_ovlp_extended, *part_unity, eigenvalue_threshold, nev, nev_arpack, shift,false,2);

           int rank = adapter.gridView().comm().rank();
           subdomainbasis->to_file(basename, rank);

           //auto fromfile_subdomainbasis;
           auto fromfile_subdomainbasis = std::make_shared<Dune::PDELab::NonoverlappingGenEOBasisFromFiles<GV, Matrix, Vector>>(adapter, basename);

           std::cout << (*subdomainbasis->get_basis_vector(0))[0] << std::endl;
           std::cout << (*fromfile_subdomainbasis->get_basis_vector(0))[0] << std::endl;
           std::cout << (*subdomainbasis->get_basis_vector(0))[0] - (*fromfile_subdomainbasis->get_basis_vector(0))[0] << std::endl;
        }

        auto coarse_space = std::make_shared<Dune::PDELab::NewSubdomainProjectedCoarseSpace<GV, Matrix, Vector>>(adapter, gv, *A_extended, subdomainbasis, verbose);

        // Apply Dirichlet conditions on processor boundaries, needed for Schwarz method
        for (auto rIt=A_extended->begin(); rIt!=A_extended->end(); ++rIt){
          for(int block_i = 0; block_i < block_size; block_i++){
            if ((*part_unity)[rIt.index()][block_i] == .0) {
              for (auto cIt=rIt->begin(); cIt!=rIt->end(); ++cIt)
              {
                for(int block_j = 0; block_j < block_size; block_j++){
                  (*cIt)[block_i][block_j] = (rIt.index() == cIt.index() && block_i == block_j) ? 1.0 : 0.0;
                }
              }
            }
          }
        }

        prec = std::make_shared<Dune::PDELab::ISTL::NonoverlappingTwoLevelOverlappingAdditiveSchwarz<GV, ScalarMatrix, Matrix, ScalarVector, Vector>>(adapter, *A_extended, coarse_space, true, verbose);

      }

      std::shared_ptr<Dune::PDELab::SubdomainBasis<Vector>> getSubdomainbasis() {
        return subdomainbasis;
      }

    private:

      std::shared_ptr<Dune::PDELab::SubdomainBasis<Vector>> subdomainbasis;

      M A_ovlp;
      std::shared_ptr<Matrix> A_ovlp_extended;
      std::shared_ptr<Matrix> A_extended;
      std::shared_ptr<Vector> part_unity;

      std::shared_ptr<Dune::PDELab::ISTL::NonoverlappingTwoLevelOverlappingAdditiveSchwarz<typename GO::Traits::TrialGridFunctionSpace::Traits::GridView, ScalarMatrix, Matrix, ScalarVector, Vector>> prec = nullptr;

    };

  }
}

#endif
