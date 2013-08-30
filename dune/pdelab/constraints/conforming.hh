// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_CONFORMINGCONSTRAINTS_HH
#define DUNE_PDELAB_CONFORMINGCONSTRAINTS_HH

#include <cstddef>
#include <algorithm>

#include <dune/common/exceptions.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/grid/common/grid.hh>

#include <dune/localfunctions/common/interfaceswitch.hh>

#include <dune/typetree/typetree.hh>

#include <dune/pdelab/common/geometrywrapper.hh>
#include <dune/pdelab/constraints/common/constraintsparameters.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/localfunctionspacetags.hh>
#include <dune/pdelab/gridfunctionspace/localvector.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup Constraints
    //! \ingroup FiniteElementMap
    //! \{

    //! Dirichlet Constraints construction
    // works in any dimension and on all element types
    class ConformingDirichletConstraints
    {
    public:
      enum { doBoundary = true };
      enum { doProcessor = false };
      enum { doSkeleton = false };
      enum { doVolume = false };

      //! boundary constraints
      /**
       * \tparam P   Parameter class, wich fulfills the DirichletConstraintsParameters interface
       * \tparam IG  intersection geometry
       * \tparam LFS local function space
       * \tparam T   TransformationType
       */
      template<typename P, typename IG, typename LFS, typename T>
      void boundary (const P& param, const IG& ig, const LFS& lfs, T& trafo) const
      {
        typedef FiniteElementInterfaceSwitch<
        typename LFS::Traits::FiniteElementType
          > FESwitch;
        typedef FieldVector<typename IG::ctype, IG::dimension-1> FaceCoord;

        const int face = ig.indexInInside();

        // find all local indices of this face
        Dune::GeometryType gt = ig.inside()->type();
        typedef typename IG::ctype DT;
        const int dim = IG::Entity::Geometry::dimension;
        const Dune::ReferenceElement<DT,dim>& refelem = Dune::ReferenceElements<DT,dim>::general(gt);

        const Dune::ReferenceElement<DT,dim-1> &
          face_refelem = Dune::ReferenceElements<DT,dim-1>::general(ig.geometry().type());

        // empty map means Dirichlet constraint
        typename T::RowType empty;

        const FaceCoord testpoint = face_refelem.position(0,0);

        // Abort if this isn't a Dirichlet boundary
        if (!param.isDirichlet(ig,testpoint))
          return;

        for (std::size_t i=0;
             i<std::size_t(FESwitch::coefficients(lfs.finiteElement()).size());
             i++)
          {
            // The codim to which this dof is attached to
            unsigned int codim =
              FESwitch::coefficients(lfs.finiteElement()).localKey(i).codim();

            if (codim==0) continue;

            for (int j=0; j<refelem.size(face,1,codim); j++){

              if (static_cast<int>(FESwitch::coefficients(lfs.finiteElement()).
                                   localKey(i).subEntity())
                  == refelem.subEntity(face,1,j,codim))
                trafo[lfs.dofIndex(i)] = empty;
            }
          }
      }
    };

    //! extend conforming constraints class by processor boundary
    class OverlappingConformingDirichletConstraints : public ConformingDirichletConstraints
    {
    public:
      enum { doProcessor = true };

      //! processor constraints
      /**
       * \tparam IG  intersection geometry
       * \tparam LFS local function space
       * \tparam T   TransformationType
       */
      template<typename IG, typename LFS, typename T>
      void processor (const IG& ig, const LFS& lfs, T& trafo) const
      {
        typedef FiniteElementInterfaceSwitch<
        typename LFS::Traits::FiniteElementType
          > FESwitch;

        // determine face
        const int face = ig.indexInInside();

        // find all local indices of this face
        Dune::GeometryType gt = ig.inside()->type();
        typedef typename IG::ctype DT;
        const int dim = IG::Entity::Geometry::dimension;

        const Dune::ReferenceElement<DT,dim>& refelem = Dune::ReferenceElements<DT,dim>::general(gt);

        // empty map means Dirichlet constraint
        typename T::RowType empty;

        // loop over all degrees of freedom and check if it is on given face
        for (size_t i=0; i<FESwitch::coefficients(lfs.finiteElement()).size();
             i++)
          {
            // The codim to which this dof is attached to
            unsigned int codim =
              FESwitch::coefficients(lfs.finiteElement()).localKey(i).codim();

            if (codim==0) continue;

            for (int j=0; j<refelem.size(face,1,codim); j++)
              if (FESwitch::coefficients(lfs.finiteElement()).localKey(i).
                  subEntity() == std::size_t(refelem.subEntity(face,1,j,codim)))
                trafo[lfs.dofIndex(i)] = empty;
          }
      }
    };

    //! extend conforming constraints class by processor boundary
    template<typename GV>
    class NonoverlappingConformingDirichletConstraints : public ConformingDirichletConstraints
    {
    public:
      enum { doVolume = true };

      //! volume constraints
      /**
       * \tparam EG  element geometry
       * \tparam LFS local function space
       * \tparam T   TransformationType
       */
      template<typename EG, typename LFS, typename T>
      void volume (const EG& eg, const LFS& lfs, T& trafo) const
      {
        typedef FiniteElementInterfaceSwitch<
        typename LFS::Traits::FiniteElementType
          > FESwitch;

        // nothing to do for interior entities
        if (eg.entity().partitionType()==Dune::InteriorEntity)
          return;

        typedef typename FESwitch::Coefficients Coefficients;
        const Coefficients& coeffs = FESwitch::coefficients(lfs.finiteElement());

        // empty map means Dirichlet constraint
        typename T::RowType empty;

        typedef typename LFS::Traits::GridFunctionSpaceType::Traits::BackendType B;

        const ReferenceElement<typename GV::ctype,GV::dimension>& ref_el =
          ReferenceElements<typename GV::ctype,GV::dimension>::general(eg.entity().type());

        // loop over all degrees of freedom and check if it is not owned by this processor
        for (size_t i = 0; i < coeffs.size(); ++i)
          {
            size_t codim = coeffs.localKey(i).codim();
            size_t sub_entity = coeffs.localKey(i).subEntity();

            size_t entity_index = _gv.indexSet().subIndex(eg.entity(),sub_entity,codim);
            size_t gt_index = GlobalGeometryTypeIndex::index(ref_el.type(sub_entity,codim));

            size_t index = _gt_offsets[gt_index] + entity_index;

            if (_ghosts[index])
              {
                trafo[lfs.dofIndex(i)] = empty;
              }
          }
      }

      template<class GFS>
      void compute_ghosts (const GFS& gfs)
      {
        std::fill(_gt_offsets.begin(),_gt_offsets.end(),0);

        typedef std::vector<GeometryType> GTVector;

        for (size_t codim = 0; codim <= GV::dimension; ++codim)
          {
            if (gfs.ordering().contains(codim))
              {
                const GTVector& geom_types = _gv.indexSet().geomTypes(codim);
                for (GTVector::const_iterator it = geom_types.begin(),
                       end = geom_types.end();
                     it != end;
                     ++it)
                  _gt_offsets[GlobalGeometryTypeIndex::index(*it) + 1] = _gv.indexSet().size(*it);
              }
          }

        std::partial_sum(_gt_offsets.begin(),_gt_offsets.end(),_gt_offsets.begin());

        _ghosts.assign(_gt_offsets.back(),true);

        typedef typename GV::template Codim<0>::
          template Partition<Interior_Partition>::Iterator Iterator;

        for(Iterator it = _gv.template begin<0, Interior_Partition>(),
              end = _gv.template end<0, Interior_Partition>();
            it != end;
            ++it)
          {
            const ReferenceElement<typename GV::ctype,GV::dimension>& ref_el =
              ReferenceElements<typename GV::ctype,GV::dimension>::general(it->type());

            for (size_t codim = 0; codim <= GV::dimension; ++codim)
              if (gfs.ordering().contains(codim))
                {
                  for (int i = 0; i < ref_el.size(codim); ++i)
                    {
                      size_t entity_index = _gv.indexSet().subIndex(*it,i,codim);
                      size_t gt_index = GlobalGeometryTypeIndex::index(ref_el.type(i,codim));
                      size_t index = _gt_offsets[gt_index] + entity_index;

                      _ghosts[index] = false;
                    }
                }
          }

      }

      NonoverlappingConformingDirichletConstraints(const GV& gv)
        : _gv(gv)
        , _rank(gv.comm().rank())
        , _gt_offsets(GlobalGeometryTypeIndex::size(GV::dimension) + 1)
      {}

    private:

      GV _gv;
      int _rank;
      std::vector<bool> _ghosts;
      std::vector<size_t> _gt_offsets;
    };
    //! \}

  }
}

#endif
