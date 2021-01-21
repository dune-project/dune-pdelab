
#ifndef DUNE_PDELAB_BACKEND_ISTL_GENEO_NEWSUBDOMAINPROJECTEDCOARSESPACE_HH
#define DUNE_PDELAB_BACKEND_ISTL_GENEO_NEWSUBDOMAINPROJECTEDCOARSESPACE_HH


#include <dune/pdelab/boilerplate/pdelab.hh>

#include <dune/common/timer.hh>

#include <dune/pdelab/backend/istl/geneo/multicommdatahandle.hh>

#include <dune/pdelab/backend/istl/geneo/coarsespace.hh>
#include <dune/pdelab/backend/istl/geneo/subdomainbasis.hh>

namespace Dune {
  namespace PDELab {


    template<typename GridView, typename Vector, typename Matrix>
    class MultiVectorBundle {
    public:
      typedef typename Vector::value_type value_type;

      MultiVectorBundle(NonoverlappingOverlapAdapter<GridView, Vector, Matrix>& adapter)
      : neighboringRanks_(adapter.findNeighboringRanks()),
        neighbor_basis(0)
      {
        for (auto& rank : neighboringRanks_) {
          neighbor_basis.push_back(std::make_shared<Vector>(adapter.getExtendedSize()));
        }
      }

      std::shared_ptr<Vector> getVectorForRank(int rank) {
        for (int i = 0; i < neighboringRanks_.size(); i++) {
          if (neighboringRanks_[i] == rank) {
            return neighbor_basis[i];
          }
        }
        DUNE_THROW(Dune::Exception, "Trying to access vector for unknown neighbor rank!");
        return nullptr;
      }

      void print() {
        for (size_t i = 0; i < neighboringRanks_.size(); i++) {
          Dune::printvector(std::cout, *(neighbor_basis[i]), "remote vec from neighbor " + std::to_string(neighboringRanks_[i]), "");
        }
      }

      std::shared_ptr<Vector> localVector_;
      std::vector<std::shared_ptr<Vector> > neighbor_basis;
      std::vector<int> neighboringRanks_;
    };

    template<typename V>
    struct MultiGatherScatter
    {
      static typename V::value_type gather(const V& a, int i)
      {
        return (*a.localVector_)[i];
      }
      static void scatter (V& a, typename V::value_type v, int i, int proc)
      {
        (*a.getVectorForRank(proc))[i]=v;
      }
    };

    /*!
    * \brief This constructs a coarse space from a per-subdomain local basis.
    *
    * The per-subdomain coarse basis is communicated to each subdomain's neighbors,
    * a global coarse system based on those is constructed and distributed to
    * across all processes. In the process, the per-subdomain basis functions are
    * extended by zeros, resulting in a sparse system.
    */
    template<class GridView, class M, class X>
    class NonoverlappingSubdomainProjectedCoarseSpace : public CoarseSpace<X>
    {

      typedef int rank_type;
    public:
      typedef typename CoarseSpace<X>::COARSE_V COARSE_V;
      typedef typename CoarseSpace<X>::COARSE_M COARSE_M;
      typedef typename M::field_type field_type;

      /*! \brief Constructor.
      * \param gfs Grid function space.
      * \param AF_exterior_ Stiffness matrix of the problem to be solved.
      * \param subdomainbasis Per-subdomain coarse basis.
      * \param verbosity Verbosity.
      */
      NonoverlappingSubdomainProjectedCoarseSpace (NonoverlappingOverlapAdapter<GridView, X, M>& adapter, const GridView& gridView, const M& AF_exterior_, std::shared_ptr<SubdomainBasis<X> > subdomainbasis, int verbosity = 1)
       : adapter_(adapter),
         gridView_(gridView),
         AF_exterior_(AF_exterior_),
         verbosity_(verbosity),
         ranks_(gridView.comm().size()),
         my_rank_(gridView.comm().rank()),
         subdomainbasis_(subdomainbasis)
      {
        neighbor_ranks_ = adapter.findNeighboringRanks();//parallelhelper.getNeighborRanks();

        setup_coarse_system();
      }

    private:
      void setup_coarse_system() {

        // Barrier for proper time measurement
        gridView_.comm().barrier();
        if (my_rank_ == 0 && verbosity_ > 0) std::cout << "Matrix setup" << std::endl;
        Dune::Timer timer_setup;

        // Communicate local coarse space dimensions
        local_basis_sizes_.resize(ranks_);
        rank_type local_size = subdomainbasis_->basis_size();
        gridView_.comm().allgather(&local_size, 1, local_basis_sizes_.data());


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
        using Attribute = EPISAttribute;
        Dune::AllSet<Attribute> allAttribute;
        auto allinterface = std::shared_ptr<Dune::Interface>(new Dune::Interface());
        allinterface->build(*adapter_.getRemoteIndices(),allAttribute,allAttribute); // all to all communication
        auto communicator = std::shared_ptr<DuneWithRank::BufferedCommunicator>(new DuneWithRank::BufferedCommunicator());
        communicator->build<X>(*allinterface);


        MultiVectorBundle<GridView, X, M> bundle(adapter_);

        // Assemble local section of coarse matrix
        for (rank_type basis_index_remote = 0; basis_index_remote < max_local_basis_size; basis_index_remote++) {

          // Communicate one basis vectors of every subdomain to all of its neighbors in one go
          // If the current rank has already communicated all its basis vectors, just pass zeros
          if (basis_index_remote < local_basis_sizes_[my_rank_]) {
            bundle.localVector_ = subdomainbasis_->get_basis_vector(basis_index_remote);
          } else {
            bundle.localVector_ = std::make_shared<X>(adapter_.getExtendedSize());
          }
          communicator->forward<MultiGatherScatter<MultiVectorBundle<GridView, X, M>>>(bundle,bundle); // make function known in other subdomains
          if (verbosity_ > 2)
            bundle.print();

          // Compute local products of basis functions with discretization matrix
          if (basis_index_remote < local_basis_sizes_[my_rank_]) {
            auto basis_vector = *subdomainbasis_->get_basis_vector(basis_index_remote);
            X Atimesv(adapter_.getExtendedSize());
            AF_exterior_.mv(basis_vector, Atimesv);
            for (rank_type basis_index = 0; basis_index < local_basis_sizes_[my_rank_]; basis_index++) {
              field_type entry = *subdomainbasis_->get_basis_vector(basis_index)*Atimesv;
              local_rows[basis_index][neighbor_ranks_.size()].push_back(entry);
            }
          }

          // Compute products of discretization matrix with local and remote vectors
          for (std::size_t neighbor_id = 0; neighbor_id < neighbor_ranks_.size(); neighbor_id++) {
            if (basis_index_remote >= local_basis_sizes_[neighbor_ranks_[neighbor_id]])
              continue;

            std::shared_ptr<X> basis_vector = bundle.getVectorForRank(neighbor_ranks_[neighbor_id]); //*bundle.neighbor_basis[neighbor_id];

            if (verbosity_ > 2) {
              std::cout << "Product of neighbor " << neighbor_ranks_[neighbor_id] << "'s basis fct " << basis_index_remote <<
                           "with local basis fcts" << std::endl;
            }

            X Atimesv(adapter_.getExtendedSize());
            AF_exterior_.mv(*basis_vector, Atimesv);

            for (rank_type basis_index = 0; basis_index < local_basis_sizes_[my_rank_]; basis_index++) {
              field_type entry=0.0;
              std::vector<double> contribution((*subdomainbasis_->get_basis_vector(basis_index)).N());
              for(int vect_iterator=0; vect_iterator<(*subdomainbasis_->get_basis_vector(basis_index)).N(); vect_iterator++){
                for(int j=0; j<3; j++){
                  contribution[vect_iterator] = (*subdomainbasis_->get_basis_vector(basis_index))[vect_iterator][j] * Atimesv[vect_iterator][j];
                  entry+= (*subdomainbasis_->get_basis_vector(basis_index))[vect_iterator][j] * Atimesv[vect_iterator][j];
                }
              }



              // if(basis_index == 0 && basis_index_remote == 0 && ((my_rank_==0 && neighbor_ranks_[neighbor_id]==1)||(my_rank_==1 && neighbor_ranks_[neighbor_id]==0))){
              if(basis_index == 0 && basis_index_remote == 0){
                std::cout << entry << " | " << my_rank_ << " | " << neighbor_ranks_[neighbor_id] << std::endl;
                // for(int it = 0; it<(*subdomainbasis_->get_basis_vector(basis_index)).N(); it++){
                //   // if ((my_rank_==0 && it==363) || (my_rank_==1 && it==60)){
                //   // if ((my_rank_==0 && it==190) || (my_rank_==1 && it==386)){
                //   // if (my_rank_==1){
                //   //   int i1=it;
                //   //   // if (std::abs(Atimesv[it][0]) > 0.0001 && std::abs(Atimesv[it][1]) > 0.0001 && std::abs(Atimesv[it][2]) > 0.0001){
                //   //   // if (std::abs(Atimesv[it]) > 0.0001){
                //   //     std::cout << my_rank_ << " | i       | " << it << std::endl;
                //   //     std::cout << my_rank_ << " | contrib | " << contribution[it] << std::endl;
                //   //     std::cout << my_rank_ << " | Phi[i]  | " << (*subdomainbasis_->get_basis_vector(basis_index))[it] << std::endl;
                //   //     std::cout << my_rank_ << " | Psy*A[i]| " << Atimesv[it] << std::endl;
                //   //     std::cout << my_rank_ << " | Psy[i]  | " << (*basis_vector)[it] << std::endl;
                //   //     std::cout << my_rank_ << " | A[i][i]  | " << AF_exterior_[i1][i1] << std::endl;
                //   //     // std::cout << i1 << " || " << i1 << std::endl;
                //   //     // std::cout << my_rank_ << " | Aii[0]  | " << AF_exterior_[i1][i1][0][0] << ", " << AF_exterior_[i1][i1][0][1] << ", " << AF_exterior_[i1][i1][0][2] << std::endl;
                //   //     // std::cout << my_rank_ << " | Aii[1]  | " << AF_exterior_[i1][i1][1][0] << ", " << AF_exterior_[i1][i1][1][1] << ", " << AF_exterior_[i1][i1][1][2] << std::endl;
                //   //     // std::cout << my_rank_ << " | Aii[2]  | " << AF_exterior_[i1][i1][2][0] << ", " << AF_exterior_[i1][i1][2][1] << ", " << AF_exterior_[i1][i1][2][2] << std::endl;
                //   //     std::cout << std::endl;
                //   //   // }
                //   // }
                // }


              }
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
            gridView_.comm().broadcast(&couplings, 1, rank);

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
            gridView_.comm().broadcast(entries_pos, couplings, rank);

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
            gridView_.comm().broadcast(entries, couplings, rank);

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


        if (my_rank_ == 0 && verbosity_ > 0) std::cout << "Coarse matrix setup finished: M=" << timer_setup.elapsed() << std::endl;
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
          for (std::size_t i = 0; i < fine.N(); i++)
            buf_defect_local[basis_index] += (*subdomainbasis_->get_basis_vector(basis_index))[i] * fine[i];
        }

        MPI_Allgatherv(&buf_defect_local, local_basis_sizes_[my_rank_], MPITraits<field_type>::getType(), &buf_defect, recvcounts, displs, MPITraits<field_type>::getType(), gridView_.comm());
        for (rank_type basis_index = 0; basis_index < global_basis_size_; basis_index++) {
          restricted[basis_index] = buf_defect[basis_index];
        }
      }

      void prolongate (const COARSE_V& coarse, X& prolongated) const override {
        prolongated = 0.0;

        // Prolongate result
        for (rank_type basis_index = 0; basis_index < local_basis_sizes_[my_rank_]; basis_index++) {
          X local_result(*subdomainbasis_->get_basis_vector(basis_index));
          local_result *= coarse[my_basis_array_offset_ + basis_index];
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

      NonoverlappingOverlapAdapter<GridView, X, M> adapter_;

      const GridView& gridView_;

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

#endif
