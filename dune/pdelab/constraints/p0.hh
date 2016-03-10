// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_CONSTRAINTS_P0_HH
#define DUNE_PDELAB_CONSTRAINTS_P0_HH

#include "../common/geometrywrapper.hh"

namespace Dune {
  namespace PDELab {

    //! \addtogroup Constraints
    //! \ingroup FiniteElementMap
    //! \{

    //! Parallel P0 constraints for overlapping grids
    class P0ParallelConstraints
    {
    public:
      enum{doBoundary=false};
      enum{doProcessor=true};
      enum{doSkeleton=false};
      enum{doVolume=false};

      //! processor constraints
      /**
       * \tparam IG  intersection geometry
       * \tparam LFS local function space
       * \tparam T   TransformationType
       */
      template<typename I, typename LFS, typename T>
      void processor (const Dune::PDELab::IntersectionGeometry<I>& ig,
                      const LFS& lfs, T& trafo) const
      {
        typename T::RowType empty;
        typedef typename LFS::Traits::SizeType size_type;
        for (size_type i=0; i<lfs.size(); i++){
          trafo[lfs.dofIndex(i)] = empty;
        }
      }
    };
    //! \}

  }
}

#endif // DUNE_PDELAB_CONSTRAINTS_P0_HH
