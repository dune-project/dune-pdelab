
#ifndef DUNE_GENEO_SUBDOMAINPROJECTEDCOARSESPACE_HH
#define DUNE_GENEO_SUBDOMAINPROJECTEDCOARSESPACE_HH


#if HAVE_ARPACKPP

#include <dune/pdelab/boilerplate/pdelab.hh>

#include <dune/common/timer.hh>

#include <dune/pdelab/preconditioner/geneo/multicommdatahandle.hh>

#include <dune/pdelab/preconditioner/geneo/coarsespace.hh>
#include <dune/pdelab/preconditioner/geneo/subdomainbasis.hh>

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
    template<class GFS, class M, class X, int dim>
    class SubdomainProjectedCoarseSpace : public CoarseSpace<M,X>
    {

    public:
      typedef typename CoarseSpace<M,X>::COARSE_V COARSE_V;
      typedef typename CoarseSpace<M,X>::COARSE_M COARSE_M;

      /*! \brief Constructor.
      * \param gfs_ Grid function space.
      * \param AF_exterior Stiffness matrix of the problem to be solved.
      * \param subdomainbasis Per-subdomain coarse basis.
      * \param verbosity Verbosity.
      */
      SubdomainProjectedCoarseSpace (const GFS& gfs_, const M& AF_exterior, std::shared_ptr<SubdomainBasis<X> > subdomainbasis, int verbosity = 0)
       : gfs(gfs_), AF_exterior(AF_exterior),
          ranks(gfs.gridView().comm().size()),
          my_rank(gfs.gridView().comm().rank()),
          subdomainbasis_(subdomainbasis),
          _verbosity(verbosity)
      {

        // Find neighbors (based on parallelhelper.hh in PDELab)

        using RankVector = Dune::PDELab::Backend::Vector<GFS,int>;
        RankVector rank_partition(gfs, my_rank); // vector to identify unique decomposition
        //! Type for storing rank values.

        //! Type used to store owner rank values of all DOFs.

        Dune::InterfaceType _interiorBorder_all_interface;

        //! The actual communication interface used when algorithm requires All_All_Interface.
        Dune::InterfaceType _all_all_interface;

        // TODO: The following shortcut may never be fulfilled because we have no overlap?
        // Let's try to be clever and reduce the communication overhead by picking the smallest
        // possible communication interface depending on the overlap structure of the GFS.
        // FIXME: Switch to simple comparison as soon as dune-grid:1b3e83ec0 is reliably available!
        if (gfs.entitySet().partitions().value == Dune::Partitions::interiorBorder.value)
          {
            // The GFS only spans the interior and border partitions, so we can skip sending or
            // receiving anything else.
            _interiorBorder_all_interface = Dune::InteriorBorder_InteriorBorder_Interface;
            _all_all_interface = Dune::InteriorBorder_InteriorBorder_Interface;
          }
        else
          {
            // In general, we have to transmit more.
            _interiorBorder_all_interface = Dune::InteriorBorder_All_Interface;
            _all_all_interface = Dune::All_All_Interface;
          }
        Dune::PDELab::DisjointPartitioningDataHandle<GFS,RankVector> pdh(gfs,rank_partition);
        gfs.gridView().communicate(pdh,_interiorBorder_all_interface,Dune::ForwardCommunication);

        std::set<int> rank_set;
        for (int rank : rank_partition)
          if (rank != my_rank)
            rank_set.insert(rank);

        for (int rank : rank_set)
          neighbor_ranks.push_back(rank);

        setup_coarse_system();
      }

    private:
      void setup_coarse_system () {
        using Dune::PDELab::Backend::native;

        // Get local basis vectors
        auto local_basis = subdomainbasis_->local_basis;

        // Normalize basis vectors
        for (int i = 0; i < local_basis.size(); i++) {
          native(*(local_basis[i])) *= 1.0 / (native(*(local_basis[i])) * native(*(local_basis[i])));
        }

        gfs.gridView().comm().barrier();
        if (my_rank == 0) std::cout << "Matrix setup" << std::endl;
        Dune::Timer timer_setup;

        // Communicate local coarse space dimensions
        int buf_basis_sizes[ranks];
        int local_size = local_basis.size();
        MPI_Allgather(&local_size, 1, MPI_INT, &buf_basis_sizes, 1, MPI_INT, gfs.gridView().comm());
        local_basis_sizes = std::vector<int>(buf_basis_sizes, buf_basis_sizes + ranks);

        // Count coarse space dimensions
        global_basis_size = 0;
        for (int n : local_basis_sizes) {
          global_basis_size += n;
        }
        my_basis_array_offset = 0;
        for (int i = 0; i < my_rank; i++) {
          my_basis_array_offset += local_basis_sizes[i];
        }

        if (my_rank == 0) std::cout << "Global basis size B=" << global_basis_size << std::endl;

        int max_local_basis_size = 0;
        for (int rank = 0; rank < ranks; rank++) {
          if (local_basis_sizes[rank] > max_local_basis_size)
            max_local_basis_size = local_basis_sizes[rank];
        }


        coarse_system = std::make_shared<COARSE_M>(global_basis_size, global_basis_size, COARSE_M::row_wise);

        std::vector<std::vector<std::vector<double> > > local_rows(local_basis_sizes[my_rank]);
        for (int basis_index = 0; basis_index < local_basis_sizes[my_rank]; basis_index++) {
          local_rows[basis_index].resize(neighbor_ranks.size()+1);
        }

        for (int basis_index_remote = 0; basis_index_remote < max_local_basis_size; basis_index_remote++) {

          std::vector<std::shared_ptr<X> > neighbor_basis(neighbor_ranks.size()); // Local coarse space basis
          for (int i = 0; i < neighbor_basis.size(); i++) {
            neighbor_basis[i] = std::make_shared<X>(gfs, 0.0);
          }

          if (basis_index_remote < local_basis_sizes[my_rank]) {
            Dune::MultiCommDataHandle<GFS,X,int,dim> commdh(gfs, *local_basis[basis_index_remote], neighbor_basis, neighbor_ranks);
            gfs.gridView().communicate(commdh,Dune::All_All_Interface,Dune::ForwardCommunication);
          } else {
            X dummy(gfs, 0.0);
            Dune::MultiCommDataHandle<GFS,X,int,dim> commdh(gfs, dummy, neighbor_basis, neighbor_ranks);
            gfs.gridView().communicate(commdh,Dune::All_All_Interface,Dune::ForwardCommunication);
          }


          if (basis_index_remote < local_basis_sizes[my_rank]) {
            auto basis_vector = *local_basis[basis_index_remote];
            X Atimesv(gfs,0.0);
            native(AF_exterior).mv(native(basis_vector), native(Atimesv));
            for (int basis_index = 0; basis_index < local_basis_sizes[my_rank]; basis_index++) {
              double entry = *local_basis[basis_index]*Atimesv;
              local_rows[basis_index][neighbor_ranks.size()].push_back(entry);
            }
          }

          for (int neighbor_id = 0; neighbor_id < neighbor_ranks.size(); neighbor_id++) {
            if (basis_index_remote >= local_basis_sizes[neighbor_ranks[neighbor_id]])
              continue;

            auto basis_vector = *neighbor_basis[neighbor_id];
            X Atimesv(gfs,0.0);
            native(AF_exterior).mv(native(basis_vector), native(Atimesv));

            for (int basis_index = 0; basis_index < local_basis_sizes[my_rank]; basis_index++) {

              double entry = *local_basis[basis_index]*Atimesv;
              local_rows[basis_index][neighbor_id].push_back(entry);
            }
          }

        }

        auto setup_row = coarse_system->createbegin();
        int row_id = 0;
        for (int rank = 0; rank < ranks; rank++) {
          for (int basis_index = 0; basis_index < local_basis_sizes[rank]; basis_index++) {

            // Communicate number of entries in this row
            int couplings = 0;
            if (rank == my_rank) {
              couplings = local_basis_sizes[my_rank];
              for (int neighbor_id : neighbor_ranks) {
                couplings += local_basis_sizes[neighbor_id];
              }
            }
            MPI_Bcast(&couplings, 1, MPI_INT, rank, gfs.gridView().comm());

            // Communicate row's pattern
            int entries_pos[couplings];
            if (rank == my_rank) {
              int cnt = 0;
              for (int basis_index2 = 0; basis_index2 < local_basis_sizes[my_rank]; basis_index2++) {
                entries_pos[cnt] = my_basis_array_offset + basis_index2;
                cnt++;
              }
              for (int neighbor_id = 0; neighbor_id < neighbor_ranks.size(); neighbor_id++) {
                int neighbor_offset = basis_array_offset (neighbor_ranks[neighbor_id]);
                for (int basis_index2 = 0; basis_index2 < local_basis_sizes[neighbor_ranks[neighbor_id]]; basis_index2++) {
                  entries_pos[cnt] = neighbor_offset + basis_index2;
                  cnt++;
                }
              }
            }
            MPI_Bcast(&entries_pos, couplings, MPI_INT, rank, gfs.gridView().comm());

            // Communicate actual entries
            double entries[couplings];
            if (rank == my_rank) {
              int cnt = 0;
              for (int basis_index2 = 0; basis_index2 < local_basis_sizes[my_rank]; basis_index2++) {
                entries[cnt] = local_rows[basis_index][neighbor_ranks.size()][basis_index2];
                cnt++;
              }
              for (int neighbor_id = 0; neighbor_id < neighbor_ranks.size(); neighbor_id++) {
                int neighbor_offset = basis_array_offset (neighbor_ranks[neighbor_id]);
                for (int basis_index2 = 0; basis_index2 < local_basis_sizes[neighbor_ranks[neighbor_id]]; basis_index2++) {
                  entries[cnt] = local_rows[basis_index][neighbor_id][basis_index2];
                  cnt++;
                }
              }
            }
            MPI_Bcast(&entries, couplings, MPI_DOUBLE, rank, gfs.gridView().comm());

            // Build matrix row based on pattern
            for (int i = 0; i < couplings; i++)
              setup_row.insert(entries_pos[i]);
            ++setup_row;

            // Set matrix entries
            for (int i = 0; i < couplings; i++) {
              (*coarse_system)[row_id][entries_pos[i]] = entries[i];
            }

            row_id++;
          }
        }


        if (my_rank == 0) std::cout << "Matrix setup finished: M=" << timer_setup.elapsed() << std::endl;
      }

      double coarse_time = 0.0;

      int basis_array_offset (int rank) {
        int offset = 0;
        for (int i = 0; i < rank; i++) {
          offset += local_basis_sizes[i];
        }
        return offset;
      }

    public:

      std::shared_ptr<COARSE_V> restrict_defect (const X& d) const override {

        auto local_basis = subdomainbasis_->local_basis;

        using Dune::PDELab::Backend::native;
        std::shared_ptr<COARSE_V> coarse_defect = std::make_shared<COARSE_V>(global_basis_size,global_basis_size);

        int recvcounts[ranks];
        int displs[ranks];
        for (int rank = 0; rank < ranks; rank++) {
          displs[rank] = 0;
        }
        for (int rank = 0; rank < ranks; rank++) {
          recvcounts[rank] = local_basis_sizes[rank];
          for (int i = rank+1; i < ranks; i++)
            displs[i] += local_basis_sizes[rank];
        }

        double buf_defect[global_basis_size];
        double buf_defect_local[local_basis_sizes[my_rank]];

        for (int basis_index = 0; basis_index < local_basis_sizes[my_rank]; basis_index++) {
          buf_defect_local[basis_index] = 0.0;
          for (int i = 0; i < native(d).N(); i++)
            buf_defect_local[basis_index] += native(*local_basis[basis_index])[i] * native(d)[i];
        }

        MPI_Allgatherv(&buf_defect_local, local_basis_sizes[my_rank], MPI_DOUBLE, &buf_defect, recvcounts, displs, MPI_DOUBLE, gfs.gridView().comm());
        for (int basis_index = 0; basis_index < global_basis_size; basis_index++) {
          (*coarse_defect)[basis_index] = buf_defect[basis_index];
        }
        return coarse_defect;
      }

      std::shared_ptr<X> prolongate_defect (const COARSE_V& v0) const override {
        auto local_basis = subdomainbasis_->local_basis;

        using Dune::PDELab::Backend::native;
        auto v = std::make_shared<X>(gfs, 0.0);

        // Prolongate result
        for (int basis_index = 0; basis_index < local_basis_sizes[my_rank]; basis_index++) {
          X local_result(*local_basis[basis_index]);
          native(local_result) *= v0[my_basis_array_offset + basis_index];
          *v += local_result;
        }
        return v;
      }

      std::shared_ptr<COARSE_M> get_coarse_system () override {
        return coarse_system;
      }

      int basis_size() override {
        return global_basis_size;
      }

    private:

      const GFS& gfs;
      const M& AF_exterior;

      std::vector<int> neighbor_ranks;

      std::shared_ptr<SubdomainBasis<X> > subdomainbasis_;

      int ranks, my_rank;
      int _verbosity;

      std::vector<int> local_basis_sizes; // Dimensions of local coarse space per subdomain
      int my_basis_array_offset; // Start of local basis functions in a consecutive global ordering
      int global_basis_size; // Dimension of entire coarse space

      std::shared_ptr<COARSE_M> coarse_system; // Coarse space matrix
    };
  }
}
#endif

#endif //DUNE_GENEO_SUBDOMAINPROJECTEDCOARSESPACE_HH
