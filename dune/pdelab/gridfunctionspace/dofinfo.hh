// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_DOFINFO_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_DOFINFO_HH

#include <cstddef>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dune/grid/common/genericreferenceelements.hh>
#include <dune/grid/common/gridenums.hh>

#include <dune/localfunctions/common/interfaceswitch.hh>
#include <dune/localfunctions/common/localkey.hh>

#include <dune/pdelab/backend/backendselector.hh>
#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>
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
      typename Dune::PDELab::BackendVectorSelector<GFS,
                                                   bool>::Type &isInterior,
      typename Dune::PDELab::BackendVectorSelector<GFS, bool>::Type &isOverlap,
      typename Dune::PDELab::BackendVectorSelector<GFS, bool>::Type &isGhost)
    {
      typedef typename GFS::Traits::GridViewType GV;
      typedef typename GV::template Codim<0>::Iterator Iterator;
      typedef typename Dune::PDELab::BackendVectorSelector<GFS, bool>::Type V;

      isInterior = false;
      isOverlap = false;
      isGhost = false;

      LocalFunctionSpace<GFS> lfs(gfs);
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

    //! get the "positions" of dofs
    /**
     * \code
#include <dune/pdelab/gridfunctionspace/dofclassification.hh>
     * \endcode
     * More precisely, determine the position of the center of their
     * associated entity.
     *
     * \note Does not work for systems currently.
     */
    template<class GFS>
    void getDofEntityPositions
    ( const GFS& gfs,
      typename Dune::PDELab::BackendVectorSelector<
        GFS,
        FieldVector<typename GFS::Traits::GridViewType::ctype,
                    GFS::Traits::GridViewType::dimensionworld>
        >::Type &dofPositions)
    {
      typedef typename GFS::Traits::GridViewType GV;
      typedef typename GV::template Codim<0>::Iterator Iterator;
      typedef typename GV::template Codim<0>::Geometry Geometry;
      typedef FieldVector<typename GV::ctype, GV::dimensionworld> Domain;
      typedef GenericReferenceElements<typename GV::ctype, GV::dimension>
        Refelems;
      typedef GenericReferenceElement<typename GV::ctype, GV::dimension>
        Refelem;
      typedef FiniteElementInterfaceSwitch<typename LocalFunctionSpace<GFS>::
                                           Traits::FiniteElementType> FESwitch;

      LocalFunctionSpace<GFS> lfs(gfs);
      LocalVector<Domain, AnySpaceTag> lv(gfs.maxLocalSize());

      const GV &gv = gfs.gridview();
      const Iterator &end = gv.template end<0>();
      for(Iterator it = gv.template begin<0>(); it != end; ++it) {
        const Refelem &refelem = Refelems::general(it->type());

        lfs.bind(*it);
        const typename FESwitch::Coefficients &coefficients =
          FESwitch::coefficients(lfs.finiteElement());
        const Geometry &geo = it->geometry();

        for(std::size_t i = 0; i < lfs.size(); ++i) {
          const LocalKey &key = coefficients.localKey(i);
          lv[i] = geo.global(refelem.position(key.subEntity(), key.codim()));
        }

        lfs.vwrite(lv, dofPositions);
      }
    }

  } // end namespace PDELab
} // end namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_DOFINFO_HH
