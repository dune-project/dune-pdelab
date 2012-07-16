// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_ENTITYCONTAINERINDEXCACHE_HH
#define DUNE_PDELAB_ENTITYCONTAINERINDEXCACHE_HH

#include <dune/pdelab/gridfunctionspace/lfscontainerindexcache.hh>

namespace Dune {
  namespace PDELab {


    template<typename GFS>
    class EntityContainerIndexCache
    {

    public:

      typedef GFS GridFunctionSpace;
      typedef typename GFS::Ordering Ordering;
      typedef typename Ordering::Traits::ContainerIndex ContainerIndex;
      typedef ContainerIndex CI;
      typedef typename Ordering::Traits::DOFIndex DOFIndex;
      typedef DOFIndex DI;
      typedef std::size_t size_type;

      typedef std::vector<CI> CIVector;

      EntityContainerIndexCache(const GFS& gfs)
        : _gfs(gfs)
        , _container_indices(gfs.maxLocalSize())
        , _size(0)
      {}

      template<typename Entity>
      void update(const Entity& e)
      {
        if (!_gfs.dataHandleContains(GFS::Traits::GridView::dimension,Entity::codimension))
          {
            _size = 0;
            return;
          }

        // clear out existing state
        for (typename CIVector::iterator it = _container_indices.begin(); it != _container_indices.end(); ++it)
          it->clear();

        _size = _gfs.dataHandleContainerIndices(e,_container_indices);
      }

      /*
      const DI& dof_index(size_type i) const
      {
        return _lfs.dofIndex(i);
      }
      */

      const CI& container_index(size_type i) const
      {
        return _container_indices[i];
      }

      /*
      const CI& container_index(const DI& i) const
      {
        // look up DOFIndex i
        std::pair<typename CIMap::iterator,bool> r = _container_index_map.insert(std::make_pair(std::ref(i),CI()));

        // i did not exist in the cache, map it into the newly inserted container index
        if (r.second)
            _lfs.gridFunctionSpace().ordering().map_index(i.view(),r.first->second);

        // return cached container index
        return r.first->second;
      }

      bool constrained(size_type i) const
      {
        return _dof_flags[i] & DOF_CONSTRAINED;
      }

      bool dirichlet_constraint(size_type i) const
      {
        return _dof_flags[i] & DOF_DIRICHLET;
      }

      ConstraintsIterator constraints_begin(size_type i) const
      {
        assert(constrained(i));
        return _constraints_iterators[i].first;
      }

      ConstraintsIterator constraints_end(size_type i) const
      {
        assert(constrained(i));
        return _constraints_iterators[i].second;
      }
      */

      const GridFunctionSpace& gridFunctionSpace() const
      {
        return _gfs;
      }

      size_type size() const
      {
        return _size;
      }

    private:

      const GFS& _gfs;
      CIVector _container_indices;
      size_type _size;

    };

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ENTITYCONTAINERINDEXCACHE_HH
