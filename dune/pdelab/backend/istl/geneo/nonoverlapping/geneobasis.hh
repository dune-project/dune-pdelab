#ifndef DUNE_PDELAB_BACKEND_ISTL_GENEO_NonoverlappingGenEOBasis_HH
#define DUNE_PDELAB_BACKEND_ISTL_GENEO_NonoverlappingGenEOBasis_HH

#include <algorithm>
#include <functional>

#include <dune/pdelab/backend/istl/geneo/subdomainbasis.hh>
#include <dune/pdelab/backend/istl/geneo/arpackpp_geneo.hh>

#if HAVE_ARPACKPP

namespace Dune {
  namespace PDELab {

    template <typename Vector, typename GV>
    class EntitySetExcluder {
    public:

      typedef typename GV::template Codim<0>::Entity Entity;
      typedef typename GV::template Codim<1>::Entity Entity1;
      typedef typename GV::template Codim<2>::Entity Entity2;

      bool includeEntity(const Entity1& entity) const {
        return true;
      }
      bool includeEntity(const Entity2& entity) const {
        return true;
      }

      virtual bool includeEntity(const Entity& entity) const {
        return true;
      }

    };

    template <typename Vector, typename GV, typename GFS>
    class EntitySetPartUnityExcluder : public EntitySetExcluder<Vector, GV> {
    public:

      typedef typename GV::template Codim<0>::Entity Entity;

      EntitySetPartUnityExcluder(const GFS& gfs, std::shared_ptr<Vector> partUnity) : lfs(gfs), partUnity_(partUnity) {}

      bool includeEntity(const Entity& entity) const override {

        if (entity.partitionType() == Dune::PartitionType::GhostEntity)
          return false;

        typedef Dune::PDELab::LFSIndexCache<LFS,Dune::PDELab::EmptyTransformation> LFSCache;
        LFSCache lfs_cache(lfs);
        lfs.bind( entity );
        lfs_cache.update();

        for (std::size_t i = 0; i < lfs_cache.size(); i++)
        {
          if ((*partUnity_)[lfs_cache.containerIndex(i).back()][0] > 0.0 &&
            (*partUnity_)[lfs_cache.containerIndex(i).back()][0] < 1.0) {
            return true;
            }
        }
        return false;
      }

    private:
      typedef Dune::PDELab::LocalFunctionSpace<GFS, Dune::PDELab::TrialSpaceTag> LFS;
      mutable LFS lfs;
      std::shared_ptr<Vector> partUnity_ = nullptr;
    };

    template <typename Vector, typename GV>
    class EntitySetGhostExcluder : public EntitySetExcluder<Vector, GV> {
    public:

      typedef typename GV::template Codim<0>::Entity Entity;

      bool includeEntity(const Entity& entity) const override {
        return entity.partitionType() != Dune::PartitionType::GhostEntity;
      }
    };

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

      /*!
       * \brief Constructor.
       * \param gfs Grid function space.
       * \param A Left-hand side matrix of the GenEO eigenproblem; with problem-dependent boundary conditions on the domain boundary
       *        and Neumann on processor boundaries
       * \param AF_ovlp The same matrix as AF_exterior, but only assembled on overlap region (where more than 1 subdomain exists).
       * \param eigenvalue_threshold Threshold up to which eigenvalue an eigenpair should be included in the basis. If negative, no thresholding.
       * \param part_unity Partition of unity to construct the basis with.
       * \param nev With thresholding, returns number of eigenvectors below threshold. Else, prescribes how many to use.
       * \param nev_arpack How many eigenpairs ARPACK is supposed to compute. Larger numbers may increase its stability. -1 for default.
       * \param shift The shift to be used in ARPACK's shift invert solver mode. May need to be adjusted for extreme eigenvalue distributions.
       * \param add_part_unity Whether to explicitly add the partition of unity itself in the coarse basis.
       * \param verbose Verbosity value.
       */
      NonoverlappingGenEOBasis(const GO& go, NonoverlappingOverlapAdapter<GridView, Vector, Matrix>& adapter, const M& A, std::shared_ptr<Matrix> A_extended, std::shared_ptr<Vector> part_unity, const double eigenvalue_threshold,
                    int& nev, int nev_arpack = -1, const double shift = 0.001, const bool add_part_unity = false, const int verbose = 0) {

        const GFS& gfs = go.trialGridFunctionSpace();
        const GV& gv = gfs.gridView();

        V x(go.trialGridFunctionSpace(),0.0); // NOTE: We assume linear problems, so simply set x to zero here!


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
        using Dune::PDELab::Backend::native;
        auto part_unity_restricted = std::make_shared<Vector>(native(A).N());
        adapter.restrictVector(*part_unity, *part_unity_restricted);



        auto es_pou_excluder = std::make_shared<EntitySetPartUnityExcluder<Vector, GV, GFS>> (gfs, part_unity_restricted);
        gfs.entitySet().setExcluder(es_pou_excluder);

        M A_ovlp(go);
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

        std::shared_ptr<Matrix> A_ovlp_extended = extended_matrices.first;
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
          //if (eigenvalues[base_id] < .0) // Normalization for debugging; Remove this for performance
          //  *(this->local_basis[base_id]) *= -1.0;
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
