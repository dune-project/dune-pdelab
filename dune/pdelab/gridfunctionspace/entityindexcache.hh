// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_ENTITYINDEXCACHE_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_ENTITYINDEXCACHE_HH

#include <dune/common/array.hh>
#include <dune/typetree/utility.hh>

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

    private:

      static const size_type leaf_count = TypeTree::TreeInfo<Ordering>::leafCount;

    public:

      typedef std::array<size_type,leaf_count + 1> Offsets;

      EntityIndexCache(const GFS& gfs)
        : _gfs(gfs)
        , _container_indices(gfs.maxLocalSize())
        , _dof_indices(map_dof_indices ? gfs.maxLocalSize() : 0)
      {
        std::fill(_offsets.begin(),_offsets.end(),0);
      }

      template<typename Entity>
      void update(const Entity& e)
      {
        std::fill(_offsets.begin(),_offsets.end(),0);
        if (!_gfs.dataHandleContains(Entity::codimension))
          return;

        _gfs.dataHandleIndices(e,_container_indices,_dof_indices,_offsets.begin(),std::integral_constant<bool,map_dof_indices>());
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
        return _offsets[leaf_count];
      }

      const Offsets& offsets() const
      {
        return _offsets;
      }

    private:

      const GFS& _gfs;
      CIVector _container_indices;
      DIVector _dof_indices;
      Offsets _offsets;

    };

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_ENTITYINDEXCACHE_HH
