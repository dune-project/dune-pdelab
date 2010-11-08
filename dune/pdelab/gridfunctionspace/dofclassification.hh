// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_DOFCLASSIFICATION_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_DOFCLASSIFICATION_HH

#include <dune/common/exceptions.hh>

#include <dune/grid/common/gridenums.hh>

#include <dune/pdelab/gridfunctionspace/localvector.hh>

namespace Dune {
  namespace PDELab {

    //! classify the dofs according to the partition type of their entities
    /**
     * \code
#include <dune/pdelab/gridfunctionspace/dofclassification.hh>
     * \endcode
     * This determines for every dof whether that dof is in or on the border
     * of an interior, overlap, or ghost grid element.  The arrays \c
     * isInterior, \c isOverlap, and \c isGhost must have beem instanciated
     * with the correct GridFunctionSpace.  They are not otherwise required to
     * be initialized in any way.
     */
    template<class GFS>
    void classifyDofs
    ( const GFS& gfs,
      typename GFS::template VectorContainer<bool>::Type &isInterior,
      typename GFS::template VectorContainer<bool>::Type &isOverlap,
      typename GFS::template VectorContainer<bool>::Type &isGhost)
    {
      typedef typename GFS::Traits::GridViewType GV;
      typedef typename GV::template Codim<0>::Iterator Iterator;
      typedef typename GFS::template VectorContainer<bool>::Type V;

      isInterior = false;
      isOverlap = false;
      isGhost = false;

      typename GFS::LocalFunctionSpace lfs(gfs);
      LocalVector<bool,AnySpaceTag> lv(gfs.maxLocalSize(), true);

      const GV &gv = gfs.gridview();
      const Iterator &end = gv.template end<0>();
      for(Iterator it = gv.template begin<0>(); it != end; ++it) {
        lfs.bind(*it);
        switch(it->partitionType()) {
        case InteriorEntity: lfs.vwrite(lv, isInterior); break;
        case OverlapEntity:  lfs.vwrite(lv, isOverlap);  break;
        case GhostEntity:    lfs.vwrite(lv, isGhost);    break;
        default:
          DUNE_THROW(Exception, "Encountered unexpected PartitionType "
                     << it->partitionType());
        }
      }
    }

  } // end namespace PDELab
} // end namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_DOFCLASSIFICATION_HH
