// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_P0CONSTRAINTS_HH
#define DUNE_PDELAB_P0CONSTRAINTS_HH

#include "../common/geometrywrapper.hh"

namespace Dune {
  namespace PDELab {

    // extend constraints class by processor boundary 
    class P0ParallelConstraints
    {
    public:
      enum{doBoundary=false};
      enum{doProcessor=true};
      enum{doSkeleton=false};
      enum{doVolume=false};

      // boundary constraints
      // IG : intersection geometry
      // LFS : local function space
      // T : TransformationType
      template<typename I, typename LFS, typename T>
      void processor (const Dune::PDELab::IntersectionGeometry<I>& ig, 
                      const LFS& lfs, T& trafo) const
      {
        typename T::RowType empty;
        typedef typename LFSU::Traits::SizeType size_type;
        for (size_type i=0; i<lfs.size(); i++)
          trafo[lfs.localIndex(i)] = empty;
      }

    };

  }
}

#endif
