
#ifndef DUNE_PDELAB_BACKEND_ISTL_GENEO_SUBDOMAINPROJECTEDCOARSESPACE_HH
#define DUNE_PDELAB_BACKEND_ISTL_GENEO_SUBDOMAINPROJECTEDCOARSESPACE_HH


#include <dune/pdelab/boilerplate/pdelab.hh>

#include <dune/common/timer.hh>

#include <dune/pdelab/backend/istl/geneo/multicommdatahandle.hh>

#include <dune/pdelab/backend/istl/geneo/coarsespace.hh>
#include <dune/pdelab/backend/istl/geneo/subdomainbasis.hh>

namespace Dune {
  namespace PDELab {

    /*!
    * \brief This constructs a coarse space from a per-subdomain local basis.
    *
    * The per-subdomain coarse basis is communicated to each subdomain's neighbors,
    * a global coarse system based on those is constructed and distributed to
    * across all processes. In the process, the per-subdomain basis functions are
    * extended by zeros, resulting in a sparse system.
    */
    template<class GFS, class M, class X, class PIH>
    class SubdomainProjectedCoarseSpace : public CoarseSpace<X>
    {

      typedef int rank_type;
    public:
      typedef typename CoarseSpace<X>::COARSE_V COARSE_V;
      typedef typename CoarseSpace<X>::COARSE_M COARSE_M;
      typedef typename COARSE_M::field_type field_type;

      /*! \brief Constructor.
      * \param gfs Grid function space.
      * \param AF_exterior_ Stiffness matrix of the problem to be solved.
      * \param subdomainbasis Per-subdomain coarse basis.
      * \param verbosity Verbosity.
      */
      SubdomainProjectedCoarseSpace (const GFS& gfs, const M& AF_exterior_, std::shared_ptr<SubdomainBasis<X> > subdomainbasis, const PIH& parallelhelper, int verbosity = 1)
       : gfs_(gfs), AF_exterior_(AF_exterior_),
         verbosity_(verbosity),
         ranks_(gfs.gridView().comm().size()),
         my_rank_(gfs.gridView().comm().rank()),
         subdomainbasis_(subdomainbasis)
      {
        neighbor_ranks_ = parallelhelper.getNeighborRanks();

        setup_coarse_system();
      }

    private:
      void setup_coarse_system() {
        using Dune::PDELab::Backend::native;

        // Barrier for proper time measurement
        gfs_.gridView().comm().barrier();
        if (my_rank_ == 0 && verbosity_ > 0) std::cout << "Matrix setup" << std::endl;
        Dune::Timer timer_setup;

        // Communicate local coarse space dimensions
        local_basis_sizes_.resize(ranks_);
        rank_type local_size = subdomainbasis_->basis_size();
        gfs_.gridView().comm().allgather(&local_size, 1, local_basis_sizes_.data());

        // Count coarse space dimensions
        global_basis_size_ = 0;
        for (rank_type n : local_basis_sizes_) {
          global_basis_size_ += n;
        }
        my_basis_array_offset_ = basis_array_offset(my_rank_);
        rank_type max_local_basis_size = *std::max_element(local_basis_sizes_.begin(),local_basis_sizes_.end());

        if (my_rank_ == 0 && verbosity_ > 0) std::cout << "Global basis size B=" << global_basis_size_ << std::endl;


        coarse_system_ = std::make_shared<COARSE_M>(global_basis_size_, global_basis_size_, COARSE_M::row_wise);

        // Set up container for storing rows of coarse matrix associated with current rank
        // Hierarchy: Own basis functions -> Current other basis function from each neighbor -> Actual entries
        std::vector<std::vector<std::vector<field_type> > > local_rows(local_basis_sizes_[my_rank_]);
        for (rank_type basis_index = 0; basis_index < local_basis_sizes_[my_rank_]; basis_index++) {
          local_rows[basis_index].resize(neighbor_ranks_.size()+1);
        }

        // Container for neighbors' basis functions
        std::vector<std::shared_ptr<X> > neighbor_basis(neighbor_ranks_.size());
        for (auto& basis : neighbor_basis) {
          basis = std::make_shared<X>(gfs_, 0.0);
        }

        // Assemble local section of coarse matrix
        for (rank_type basis_index_remote = 0; basis_index_remote < max_local_basis_size; basis_index_remote++) {

          // Communicate one basis vectors of every subdomain to all of its neighbors in one go
          // If the current rank has already communicated all its basis vectors, just pass zeros
          if (basis_index_remote < local_basis_sizes_[my_rank_]) {
            Dune::PDELab::MultiCommDataHandle<GFS,X,rank_type> commdh(gfs_, *subdomainbasis_->get_basis_vector(basis_index_remote), neighbor_basis, neighbor_ranks_);
            gfs_.gridView().communicate(commdh,Dune::All_All_Interface,Dune::ForwardCommunication);
          } else {
            X dummy(gfs_, 0.0);
            Dune::PDELab::MultiCommDataHandle<GFS,X,rank_type> commdh(gfs_, dummy, neighbor_basis, neighbor_ranks_);
            gfs_.gridView().communicate(commdh,Dune::All_All_Interface,Dune::ForwardCommunication);
          }

          // Compute local products of basis functions with discretization matrix
          if (basis_index_remote < local_basis_sizes_[my_rank_]) {
            auto basis_vector = *subdomainbasis_->get_basis_vector(basis_index_remote);
            X Atimesv(gfs_,0.0);
            native(AF_exterior_).mv(native(basis_vector), native(Atimesv));
            for (rank_type basis_index = 0; basis_index < local_basis_sizes_[my_rank_]; basis_index++) {
              field_type entry = *subdomainbasis_->get_basis_vector(basis_index)*Atimesv;
              local_rows[basis_index][neighbor_ranks_.size()].push_back(entry);
            }
          }

          // Compute products of discretization matrix with local and remote vectors
          for (std::size_t neighbor_id = 0; neighbor_id < neighbor_ranks_.size(); neighbor_id++) {
            if (basis_index_remote >= local_basis_sizes_[neighbor_ranks_[neighbor_id]])
              continue;

            auto basis_vector = *neighbor_basis[neighbor_id];
            X Atimesv(gfs_,0.0);
            native(AF_exterior_).mv(native(basis_vector), native(Atimesv));

            for (rank_type basis_index = 0; basis_index < local_basis_sizes_[my_rank_]; basis_index++) {

              field_type entry = *subdomainbasis_->get_basis_vector(basis_index)*Atimesv;
              local_rows[basis_index][neighbor_id].push_back(entry);
            }
          }

        }

        // Construct coarse matrix from local sections
        auto setup_row = coarse_system_->createbegin();
        rank_type row_id = 0;
        for (rank_type rank = 0; rank < ranks_; rank++) {
          for (rank_type basis_index = 0; basis_index < local_basis_sizes_[rank]; basis_index++) {

            // Communicate number of entries in this row
            rank_type couplings = 0;
            if (rank == my_rank_) {
              couplings = local_basis_sizes_[my_rank_];
              for (rank_type neighbor_id : neighbor_ranks_) {
                couplings += local_basis_sizes_[neighbor_id];
              }
            }
            gfs_.gridView().comm().broadcast(&couplings, 1, rank);

            // Communicate row's pattern
            rank_type entries_pos[couplings];
            if (rank == my_rank_) {
              rank_type cnt = 0;
              for (rank_type basis_index2 = 0; basis_index2 < local_basis_sizes_[my_rank_]; basis_index2++) {
                entries_pos[cnt] = my_basis_array_offset_ + basis_index2;
                cnt++;
              }
              for (std::size_t neighbor_id = 0; neighbor_id < neighbor_ranks_.size(); neighbor_id++) {
                rank_type neighbor_offset = basis_array_offset (neighbor_ranks_[neighbor_id]);
                for (rank_type basis_index2 = 0; basis_index2 < local_basis_sizes_[neighbor_ranks_[neighbor_id]]; basis_index2++) {
                  entries_pos[cnt] = neighbor_offset + basis_index2;
                  cnt++;
                }
              }
            }
            gfs_.gridView().comm().broadcast(entries_pos, couplings, rank);

            // Communicate actual entries
            field_type entries[couplings];
            if (rank == my_rank_) {
              rank_type cnt = 0;
              for (rank_type basis_index2 = 0; basis_index2 < local_basis_sizes_[my_rank_]; basis_index2++) {
                entries[cnt] = local_rows[basis_index][neighbor_ranks_.size()][basis_index2];
                cnt++;
              }
              for (std::size_t neighbor_id = 0; neighbor_id < neighbor_ranks_.size(); neighbor_id++) {
                for (rank_type basis_index2 = 0; basis_index2 < local_basis_sizes_[neighbor_ranks_[neighbor_id]]; basis_index2++) {
                  entries[cnt] = local_rows[basis_index][neighbor_id][basis_index2];
                  cnt++;
                }
              }
            }
            gfs_.gridView().comm().broadcast(entries, couplings, rank);

            // Build matrix row based on pattern
            for (rank_type i = 0; i < couplings; i++)
              setup_row.insert(entries_pos[i]);
            ++setup_row;

            // Set matrix entries
            for (rank_type i = 0; i < couplings; i++) {
              (*coarse_system_)[row_id][entries_pos[i]] = entries[i];
            }

            row_id++;
          }
        }


        if (my_rank_ == 0 && verbosity_ > 0) std::cout << "Matrix setup finished: M=" << timer_setup.elapsed() << std::endl;
      }

      /*! \brief Returns the offset of the block of local coarse basis functions w.r.t. global ordering
       */
      rank_type basis_array_offset (rank_type rank) {
        rank_type offset = 0;
        for (rank_type i = 0; i < rank; i++) {
          offset += local_basis_sizes_[i];
        }
        return offset;
      }

    public:

      void restrict (const X& fine, COARSE_V& restricted) const override {

        using Dune::PDELab::Backend::native;

        rank_type recvcounts[ranks_];
        rank_type displs[ranks_];
        for (rank_type rank = 0; rank < ranks_; rank++) {
          displs[rank] = 0;
        }
        for (rank_type rank = 0; rank < ranks_; rank++) {
          recvcounts[rank] = local_basis_sizes_[rank];
          for (rank_type i = rank+1; i < ranks_; i++)
            displs[i] += local_basis_sizes_[rank];
        }

        field_type buf_defect[global_basis_size_];
        field_type buf_defect_local[local_basis_sizes_[my_rank_]];

        for (rank_type basis_index = 0; basis_index < local_basis_sizes_[my_rank_]; basis_index++) {
          buf_defect_local[basis_index] = 0.0;
          for (std::size_t i = 0; i < native(fine).N(); i++)
            buf_defect_local[basis_index] += native(*subdomainbasis_->get_basis_vector(basis_index))[i] * native(fine)[i];
        }

        MPI_Allgatherv(&buf_defect_local, local_basis_sizes_[my_rank_], MPITraits<field_type>::getType(), &buf_defect, recvcounts, displs, MPITraits<field_type>::getType(), gfs_.gridView().comm());
        for (rank_type basis_index = 0; basis_index < global_basis_size_; basis_index++) {
          restricted[basis_index] = buf_defect[basis_index];
        }
      }

      void prolongate (const COARSE_V& coarse, X& prolongated) const override {
        using Dune::PDELab::Backend::native;

        prolongated = 0.0;

        // Prolongate result
        for (rank_type basis_index = 0; basis_index < local_basis_sizes_[my_rank_]; basis_index++) {
          X local_result(*subdomainbasis_->get_basis_vector(basis_index));
          native(local_result) *= coarse[my_basis_array_offset_ + basis_index];
          prolongated += local_result;
        }
      }

      std::shared_ptr<COARSE_M> get_coarse_system () override {
        return coarse_system_;
      }

      rank_type basis_size() override {
        return global_basis_size_;
      }

    private:

      const GFS& gfs_;
      const M& AF_exterior_;

      int verbosity_;

      std::vector<rank_type> neighbor_ranks_;

      rank_type ranks_, my_rank_;

      std::shared_ptr<SubdomainBasis<X> > subdomainbasis_;

      std::vector<rank_type> local_basis_sizes_; // Dimensions of local coarse space per subdomain
      rank_type my_basis_array_offset_; // Start of local basis functions in a consecutive global ordering
      rank_type global_basis_size_; // Dimension of entire coarse space

      std::shared_ptr<COARSE_M> coarse_system_; // Coarse space matrix
    };
  }
}

#endif //DUNE_PDELAB_BACKEND_ISTL_GENEO_SUBDOMAINPROJECTEDCOARSESPACE_HH
