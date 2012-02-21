// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_DIRECTLEAFLOCALORDERING_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_DIRECTLEAFLOCALORDERING_HH

#include <dune/pdelab/common/typetree/leafnode.hh>
#include <dune/pdelab/gridfunctionspace/orderingutility.hh>

#include <vector>

namespace Dune {
  namespace PDELab {

    template<typename FEM, typename GV, typename DI, typename CI>
    class DirectLeafLocalOrdering
      : public TypeTree::LeafNode
    {

      template<typename>
      friend class LeafGridViewOrdering;

    public:

      typedef LocalOrderingTraits<GV,DI,CI> Traits;

      void map_local_index(const typename Traits::SizeType geometry_type_index,
                           const typename Traits::SizeType entity_index,
                           typename Traits::TreeIndexView mi,
                           typename Traits::ContainerIndex& ci) const
      {
        DUNE_THROW(NotImplemented,"not implemented");
      }

      template<typename ItIn, typename ItOut>
      void map_indices(const ItIn begin, const ItIn end, ItOut out) const
      {
        // don't do anything - this is handled by the specialized GridViewOrdering
      }

      typename Traits::SizeType size(const typename Traits::SizeType geometry_type_index, const typename Traits::SizeType entity_index) const
      {
        typedef typename Traits::SizeType size_type;
        if (_fixed_size)
          return _gt_dof_offsets[geometry_type_index];
        else if (_gt_used[geometry_type_index])
          {
            const size_type index = _gt_entity_offsets[geometry_type_index] + entity_index;
            return _entity_dof_offsets[index+1] - _entity_dof_offsets[index];
          }
        else
          return 0;
      }

      typename Traits::SizeType size(const typename Traits::SizeType geometry_type_index, const typename Traits::SizeType entity_index, const typename Traits::SizeType child_index) const
      {
        DUNE_THROW(NotImplemented,"not implemented");
      }

      typename Traits::SizeType offset(const typename Traits::SizeType geometry_type_index, const typename Traits::SizeType entity_index, const typename Traits::SizeType child_index) const
      {
        assert(child_index == 0);
        return 0;
      }

      explicit DirectLeafLocalOrdering(const shared_ptr<const FEM>& fem, const GV& gv)
        : _fem(fem)
        , _gv(gv)
        , _fixed_size(false)
        , _container_blocked(false)
      {}

      const typename Traits::GridView& gridView() const
      {
        return _gv;
      }

      const FEM& finiteElementMap() const
      {
        return *_fem;
      }

    private:

      void update_a_priori_fixed_size()
      {
        _fixed_size = _fem->fixedSize();
      }

      typedef std::vector<GeometryType> GTVector;

      void update_fixed_size(const GTVector& geom_types)
      {
        assert(_fixed_size);

        typedef typename Traits::SizeType size_type;
        const size_type dim = Traits::GridView::dimension;
        _codim_used.assign(dim,false);
        _gt_used.assign(GlobalGeometryTypeIndex::size(dim),false);
        _gt_dof_offsets.assign(GlobalGeometryTypeIndex::size(dim),0);
        for (GTVector::const_iterator it = geom_types.begin(); it != geom_types.end(); ++it)
          {
            size_type size = _fem->size(*it);
            _gt_dof_offsets[GlobalGeometryTypeIndex::index(*it)] = size;
            _gt_used[GlobalGeometryTypeIndex::index(*it)] = size > 0;
            _codim_used[dim - it->dim()] = _codim_used[dim - it->dim()] || (size > 0);
          }
      }

      typename Traits::SizeType maxLocalSize() const
      {
        return _fem->maxLocalSize();
      }

    protected:

      shared_ptr<const FEM> _fem;
      GV _gv;
      bool _fixed_size;
      const bool _container_blocked;

      std::vector<bool> _codim_used;
      std::vector<bool> _gt_used;

      std::vector<typename Traits::SizeType> _gt_entity_offsets;
      std::vector<typename Traits::SizeType> _gt_dof_offsets;
      std::vector<typename Traits::SizeType> _entity_dof_offsets;

    };


  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_DIRECTLEAFLOCALORDERING_HH
