#ifndef DUNE_PDELAB_BACKEND_ISTL_GENEO_MULTICOMMDATAHANDLE_HH
#define DUNE_PDELAB_BACKEND_ISTL_GENEO_MULTICOMMDATAHANDLE_HH

#include <dune/pdelab/gridfunctionspace/entityindexcache.hh>

namespace Dune {
  namespace PDELab {

    /*!
     * \brief Gather/scatter communication that passes a single function from each
     * subdomain to all its neighbors.
     */
    template<typename GFS, typename RankIndex, typename V>
    class MultiCommGatherScatter
    {

      typedef Dune::PDELab::EntityIndexCache<GFS> IndexCache;
      typedef typename V::template LocalView<IndexCache> LocalView;

    public:

      typedef std::size_t size_type;

      template<typename MessageBuffer, typename Entity, typename LocalView>
      bool gather(MessageBuffer& buff, const Entity& e, LocalView& local_view) const
      {
        // Write values
        for (std::size_t i = 0; i < local_view.size(); ++i) {
          buff.write (local_view[i]);
        }
        return false;
      }

      template<typename MessageBuffer, typename Entity, typename LocalView>
      bool scatter(MessageBuffer& buff, std::size_t n, const Entity& e, LocalView& local_view) const
      {
        RankIndex remote_rank = buff.senderRank();

        // Get neighbor index from rank
        int remote_index = std::distance(_neighbor_ranks.begin(), std::find(_neighbor_ranks.begin(), _neighbor_ranks.end(), remote_rank));
        if (remote_index == static_cast<RankIndex>(_neighbor_ranks.size()))
          DUNE_THROW(Exception,"Received remote rank " << remote_rank << ", but it's not in the given neighbor set!");

        // Get values
        auto target_view = LocalView(*_target_vectors[remote_index]);
        _index_cache.update(e);
        target_view.bind(_index_cache);


        if (target_view.cache().gridFunctionSpace().entitySet().partitions().contains(e.partitionType())) {
          if (target_view.size() != n)
            DUNE_THROW(Exception,"size mismatch in GridFunctionSpace data handle, have " << target_view.size() << "DOFs, but received " << n);

          for (size_type i = 0; i < target_view.size(); ++i) {
            typename LocalView::ElementType x;
            buff.read(x);
            target_view[i] = x;
          }
          target_view.unbind();
          return true;

        } else {

          if (target_view.size() != 0)
            DUNE_THROW(Exception,"expected no DOFs in partition '" << e.partitionType() << "', but have " << target_view.size());

          for (size_type i = 0; i < target_view.size(); ++i) {
            typename LocalView::ElementType dummy;
            buff.read(dummy);
          }
          target_view.unbind();
          return false;
        }
      }

      // enhanced scatter with support for function spaces with different structure on sender and receiver side
      template<typename MessageBuffer, typename Offsets, typename Entity, typename LocalView>
      bool scatter(MessageBuffer& buff, const Offsets& remote_offsets, const Offsets& local_offsets, const Entity& e, LocalView& local_view) const
      {
        if (local_view.cache().gridFunctionSpace().entitySet().partitions().contains(e.partitionType()))
        {
          // the idea here is this:
          // the compile time structure of the overall function space (and its ordering) will be identical on both sides
          // of the communication, but one side may be missing information on some leaf orderings, e.g. because the DOF
          // belongs to a MultiDomain subdomain that only has an active grid cell on one side of the communication.
          // So we step through the leaves and simply ignore any block where one of the two sides is of size 0.
          // Otherwise, it's business as usual: we make sure that the sizes from both sides match and then process all
          // data with the DOF-local gather / scatter functor.

          RankIndex remote_rank = buff.senderRank();

          // Get neighbor index from rank
          int remote_index = std::distance(_neighbor_ranks.begin(), std::find(_neighbor_ranks.begin(), _neighbor_ranks.end(), remote_rank));
          if (remote_index == static_cast<RankIndex>(_neighbor_ranks.size()))
            DUNE_THROW(Exception,"Received remote rank " << remote_rank << ", but it's not in the given neighbor set!");

          // Get values
          auto target_view = LocalView(*_target_vectors[remote_index]);
          _index_cache.update(e);
          target_view.bind(_index_cache);


          size_type remote_i = 0;
          size_type local_i = 0;
          bool needs_commit = false;
          for (size_type block = 1; block < local_offsets.size(); ++block)
          {
            // we didn't get any data - just ignore
            if (remote_offsets[block] == remote_i)
            {
              local_i = local_offsets[block];
              continue;
            }

            // we got data for DOFs we don't have - drain buffer
            if (local_offsets[block] == local_i)
            {
              for (; remote_i < remote_offsets[block]; ++remote_i)
                {
                  typename LocalView::ElementType dummy;
                  buff.read(dummy);
                }
              continue;
            }

            if (remote_offsets[block] - remote_i != local_offsets[block] - local_i)
              DUNE_THROW(Exception,"size mismatch in GridFunctionSpace data handle block " << block << ", have " << local_offsets[block] - local_i << "DOFs, but received " << remote_offsets[block] - remote_i);

            for (; local_i < local_offsets[block]; ++local_i) { // TODO: Correct?
              typename LocalView::ElementType x;
              buff.read(x);
              target_view[local_i] = x;
            }

            remote_i = remote_offsets[block];
            needs_commit = true;
          }
          target_view.unbind();
          return needs_commit;
        }
        else
        {
          if (local_view.size() != 0)
            DUNE_THROW(Exception,"expected no DOFs in partition '" << e.partitionType() << "', but have " << local_view.size());

          for (std::size_t i = 0; i < remote_offsets.back(); ++i)
          {
            typename LocalView::ElementType dummy;
            buff.read(dummy);
          }
          return false;
        }
      }

      /*!
       * \param gfs Grid function space to be operated on.
       * \param rank  The MPI rank of the current process.
       * \param target_vectors Vectors which the received vectors will be written to.
       * \param neighbor_ranks List of ranks of neighboring subdomains defining the order in which results are written to target_vectors.
       */
      MultiCommGatherScatter(const GFS& gfs, RankIndex rank,
                             std::vector<std::shared_ptr<V> > target_vectors,
                             std::vector<RankIndex> neighbor_ranks)
        : _index_cache(gfs)
        , _rank(rank)
        , _target_vectors(target_vectors)
        , _neighbor_ranks(neighbor_ranks)
      {}

    private:

      mutable IndexCache _index_cache;

      const RankIndex _rank;
      std::vector<std::shared_ptr<V> > _target_vectors;
      std::vector<RankIndex> _neighbor_ranks;
    };

    /*!
     * \brief Gather/scatter communication that passes a single function from each
     * subdomain to all its neighbors. The functions are received individually
     * without applying any further operations.
     */
    template<class GFS, class V, typename RankIndex>
    class MultiCommDataHandle
    : public Dune::PDELab::GFSDataHandle<GFS,V,
                                         MultiCommGatherScatter<GFS, RankIndex, V>,
                                         Dune::PDELab::DOFDataCommunicationDescriptor<typename V::ElementType,true>>
    {
      typedef Dune::PDELab::GFSDataHandle<GFS,V,
                                          MultiCommGatherScatter<GFS, RankIndex, V>,
                                          Dune::PDELab::DOFDataCommunicationDescriptor<typename V::ElementType,true>> BaseT;

    public:

      /*!
       * \param gfs_ Grid function space to be operated on.
       * \param v_ The current rank's vector to be passed to its neighbors.
       * \param target_vectors Vectors which the received vectors will be written to.
       * \param neighbor_ranks List of ranks of neighboring subdomains defining the order in which results are written to target_vectors.
       */
      MultiCommDataHandle(const GFS& gfs_, V& v_, std::vector<std::shared_ptr<V> > target_vectors, std::vector<RankIndex> neighbor_ranks)
      : BaseT(gfs_,v_,MultiCommGatherScatter<GFS, RankIndex, V>(gfs_, gfs_.gridView().comm().rank(), target_vectors, neighbor_ranks),
              Dune::PDELab::DOFDataCommunicationDescriptor<typename V::ElementType,true>())
      {
      }
    };
  }
}

#endif //DUNE_PDELAB_BACKEND_ISTL_GENEO_MULTICOMMDATAHANDLE_HH
