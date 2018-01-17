#ifndef DUNE_MULTICOMMDATAHANDLE_HH
#define DUNE_MULTICOMMDATAHANDLE_HH

#include <dune/pdelab/gridfunctionspace/entityindexcache.hh>

namespace Dune {

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

    template<typename MessageBuffer, typename Entity, typename LocalView>
    bool gather(MessageBuffer& buff, const Entity& e, LocalView& local_view) const
    {
      // Write rank
      buff.write((double)_rank);

      // Write values
      for (std::size_t i = 0; i < local_view.size(); ++i) {
        buff.write (local_view[i]);
      }
      return true;
    }

    template<typename MessageBuffer, typename Entity, typename LocalView>
    bool scatter(MessageBuffer& buff, std::size_t n, const Entity& e, LocalView& local_view) const
    {
      double received_rank;
      buff.read(received_rank); // read in original type!
      RankIndex remote_rank = received_rank;

      // Get neighbor index from rank
      int remote_index = -1;
      for (int i = 0; i < _neighbor_ranks.size(); i++) {
        if (_neighbor_ranks[i] == remote_rank) {
          remote_index = i;
          break;
        }
      }
      if (remote_index == -1)
        DUNE_THROW(Exception,"Received remote rank " << remote_rank << ", but it's not in the given neighbor set!");

      // Get values
      auto target_view = LocalView(*_target_vectors[remote_index]);
      _index_cache.update(e);
      target_view.bind(_index_cache);


      if (target_view.cache().gridFunctionSpace().entitySet().partitions().contains(e.partitionType())) {
        if (target_view.size() != n - 1) // 1st entry for rank, n-1 for DOFs
          DUNE_THROW(Exception,"size mismatch in GridFunctionSpace data handle, have " << target_view.size() << "DOFs, but received " << n);

        for (std::size_t i = 0; i < target_view.size(); ++i) {
          typename LocalView::ElementType x;
          buff.read(x);
          target_view[i] = x;
        }
        target_view.unbind();
        return true;

      } else {

        if (target_view.size() != 0)
          DUNE_THROW(Exception,"expected no DOFs in partition '" << e.partitionType() << "', but have " << target_view.size());

        for (std::size_t i = 0; i < target_view.size(); ++i) {
          typename LocalView::ElementType dummy;
          buff.read(dummy);
        }
        target_view.unbind();
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
     : _rank(rank), _target_vectors(target_vectors), _neighbor_ranks(neighbor_ranks), _index_cache(gfs)
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
  template<class GFS, class V, typename RankIndex,int dim>
  class MultiCommDataHandle
  : public Dune::PDELab::GFSDataHandle<GFS,V,
                                       MultiCommGatherScatter<GFS, RankIndex, V>,
                                       Dune::PDELab::EntityDataCommunicationDescriptor<typename V::ElementType>>
  {
    typedef Dune::PDELab::GFSDataHandle<GFS,V,
                                        MultiCommGatherScatter<GFS, RankIndex, V>,
                                        Dune::PDELab::EntityDataCommunicationDescriptor<typename V::ElementType>> BaseT;

  public:

    /*!
     * \param gfs_ Grid function space to be operated on.
     * \param v_ The current rank's vector to be passed to its neighbors.
     * \param target_vectors Vectors which the received vectors will be written to.
     * \param neighbor_ranks List of ranks of neighboring subdomains defining the order in which results are written to target_vectors.
     */
    MultiCommDataHandle(const GFS& gfs_, V& v_, std::vector<std::shared_ptr<V> > target_vectors, std::vector<RankIndex> neighbor_ranks)
    : BaseT(gfs_,v_,MultiCommGatherScatter<GFS, RankIndex, V>(gfs_, gfs_.gridView().comm().rank(), target_vectors, neighbor_ranks),
            Dune::PDELab::EntityDataCommunicationDescriptor<typename V::ElementType>(dim+1)) // TODO: Size request derived from element!
    {
    }
  };

}

#endif //DUNE_MULTICOMMDATAHANDLE_HH
