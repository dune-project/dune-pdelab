// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_LEAFORDERING_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_LEAFORDERING_HH

#include <cstddef>

#include <dune/pdelab/common/typetree/leafnode.hh>
#include <dune/pdelab/gridfunctionspace/orderingdynamicbase.hh>
#include <dune/pdelab/gridfunctionspace/gridviewordering.hh>
#include <dune/pdelab/gridfunctionspace/leaflocalordering.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup GridFunctionSpace
    //! \ingroup PDELab
    //! \{

    template<typename LeafGFS, typename Transformation>
    struct leaf_gfs_to_ordering_descriptor
    {

      static const bool recursive = false;

      typedef LeafLocalOrdering<LeafGFS,
                                typename Transformation::MultiIndex,
                                typename Transformation::ContainerIndex
                                > LocalOrdering;

      typedef LeafGridViewOrdering<typename LeafGFS::Traits::GridView,LocalOrdering> GridViewOrdering;

      typedef LeafOrdering<GridViewOrdering> transformed_type;
      typedef shared_ptr<transformed_type> transformed_storage_type;

      static transformed_type transform(const LeafGFS& s, const Transformation& t)
      {
        return transformed_type(make_tuple(make_shared<GridViewOrdering>(gfs.gridview(),make_tuple(make_shared<LocalOrdering>(gfs)))));
      }

      static transformed_storage_type transform_storage(shared_ptr<const LeafGFS> gfs, const Transformation& t)
      {
        return make_shared<transformed_type>(make_tuple(make_shared<GridViewOrdering>(gfs.gridview(),make_tuple(make_shared<LocalOrdering>(gfs)))));
      }

    };

    template<typename LeafGFS, typename Params>
    leaf_gfs_to_ordering_descriptor<LeafGFS,gfs_to_ordering<Params> >
    lookupNodeTransformation(LeafGFS* gfs, gfs_to_ordering<Params>* t, LeafGridFunctionSpaceTag tag);

   //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_LEAFORDERING_HH
