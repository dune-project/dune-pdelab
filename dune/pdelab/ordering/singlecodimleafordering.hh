// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_ORDERING_SINGLECODIMLEAFORDERING_HH
#define DUNE_PDELAB_ORDERING_SINGLECODIMLEAFORDERING_HH

#include <dune/typetree/leafnode.hh>

#include <dune/pdelab/ordering/utility.hh>
#include <dune/pdelab/gridfunctionspace/tags.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup Ordering
    //! \{

    template<typename GV, typename DI, typename CI>
    class SingleCodimLeafOrdering
      : public TypeTree::LeafNode
    {

    public:

      typedef FlatContainerAllocationTag ContainerAllocationTag;

      typedef SimpleLFSCacheTag CacheTag;


      typedef SimpleOrderingTraits<DI,CI> Traits;

      typename Traits::ContainerIndex mapIndex(const typename Traits::DOFIndex& di) const
      {
        return di[0];
      }

      void mapIndex(typename Traits::DOFIndex di, typename Traits::ContainerIndex& ci) const
      {
        ci = di[0];
      }

      typename Traits::SizeType size() const
      {
        return _gv.size(0);
      }

      typename Traits::SizeType blockCount() const
      {
        return size();
      }

      typename Traits::SizeType maxLocalSize() const
      {
        return 1;
      }

      void update()
      {
      }

      SingleCodimLeafOrdering(const GV& gv)
        : _gv(gv)
      {
      }

      bool container_blocked() const
      {
        return false;
      }

    private:

      GV _gv;

    };


    template<typename GFS, typename Transformation>
    struct leaf_gfs_to_ordering_descriptor<GFS,Transformation,SingleCodimMapper>
    {

      static const bool recursive = false;

      typedef SingleCodimLeafOrdering<
        typename GFS::Traits::GridView,
        SimpleDOFIndex<typename GFS::Traits::SizeType>,
        SimpleContainerIndex<typename GFS::Traits::SizeType>
        > transformed_type;

      typedef std::shared_ptr<transformed_type> transformed_storage_type;

      static transformed_type transform(const GFS& gfs, const Transformation& t)
      {
        return transformed_type(gfs.gridView());
      }

      static transformed_storage_type transform_storage(std::shared_ptr<const GFS> gfs, const Transformation& t)
      {
        return std::make_shared<transformed_type>(gfs->gridView());
      }

    };

    //! \} group ordering

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ORDERING_SINGLECODIMLEAFORDERING_HH
