// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_ENTITYINDEXCACHE_HH
#define DUNE_PDELAB_ENTITYINDEXCACHE_HH

#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>

namespace Dune {
  namespace PDELab {


    template<typename GFS, bool map_dof_indices = false>
    class EntityIndexCache
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
      typedef std::vector<DI> DIVector;

      EntityIndexCache(const GFS& gfs)
        : _gfs(gfs)
        , _container_indices(gfs.maxLocalSize())
        , _dof_indices(map_dof_indices ? gfs.maxLocalSize() : 0)
        , _size(0)
      {}

      template<typename Entity>
      void update(const Entity& e)
      {
        if (!_gfs.dataHandleContains(Entity::codimension))
          {
            _size = 0;
            return;
          }

        _size = _gfs.dataHandleIndices(e,_container_indices,_dof_indices,std::integral_constant<bool,map_dof_indices>());
      }

      const DI& dofIndex(size_type i) const
      {
        assert(map_dof_indices);
        return _dof_indices[i];
      }

      const CI& containerIndex(size_type i) const
      {
        return _container_indices[i];
      }

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
      DIVector _dof_indices;
      size_type _size;

    };

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ENTITYINDEXCACHE_HH
