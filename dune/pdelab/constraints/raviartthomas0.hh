// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
#ifndef DUNE_PDELAB_CONSTRAINTS_RAVIARTTHOMAS0_HH
#define DUNE_PDELAB_CONSTRAINTS_RAVIARTTHOMAS0_HH

#include<dune/typetree/typetree.hh>

#include<dune/common/exceptions.hh>
#include<dune/geometry/referenceelements.hh>
#include<dune/geometry/type.hh>
#include<dune/pdelab/common/geometrywrapper.hh>

#include <dune/typetree/typetree.hh>

namespace Dune {
  namespace PDELab {

    //! Neumann Constraints construction, as needed for RT0
    class RT0Constraints
      : public TypeTree::LeafNode
    {
    public:
      enum{doBoundary=true};
      enum{doProcessor=false};
      enum{doSkeleton=false};
      enum{doVolume=false};


      //! boundary constraints
      /**
       * \tparam P   Parameter class, wich fulfills the FluxConstraintsParameters interface
       * \tparam IG  intersection geometry
       * \tparam LFS local function space
       * \tparam T   TransformationType
       */
      template<typename P, typename IG, typename LFS, typename T>
      void boundary (const P& p, const IG& ig, const LFS& lfs, T& trafo) const
      {
        const int face = ig.indexInInside();
        auto face_refelem = referenceElement(ig.geometry());
        const auto ip = face_refelem.position(0,0);
        if (p.isNeumann(ig,ip)) {
          typename T::RowType empty;              // need not interpolate
          trafo[lfs.dofIndex(face)]=empty;
        }
      }
    };

  }
}

#endif // DUNE_PDELAB_CONSTRAINTS_RAVIARTTHOMAS0_HH
