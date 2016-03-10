// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_CONSTRAINTS_P0GHOST_HH
#define DUNE_PDELAB_CONSTRAINTS_P0GHOST_HH

#include "../common/geometrywrapper.hh"
#include<dune/grid/common/gridenums.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup Constraints
    //! \ingroup FiniteElementMap
    //! \{

    //! Parallel P0 constraints for nonoverlapping grids with ghosts
    class P0ParallelGhostConstraints
    {
    public:
      enum{doBoundary=false};
      enum{doProcessor=false};
      enum{doSkeleton=false};
      enum{doVolume=true};

      //! volume constraints
      /**
       * \tparam EG  element geometry
       * \tparam LFS local function space
       * \tparam T   TransformationType
       */

      template<typename P, typename EG, typename LFS, typename T>
      void volume (const P& param, const EG& eg, const LFS& lfs, T& trafo) const
      {
        // nothing to do for interior entities
        if (eg.entity().partitionType()==Dune::InteriorEntity)
          return;

        // constrain ghost entities
        else if  (eg.entity().partitionType()==Dune::GhostEntity){
          typename T::RowType empty;
          typedef typename LFS::Traits::SizeType size_type;
          for (size_type i=0; i<lfs.size(); i++){
            trafo[lfs.dofIndex(i)] = empty;
          }
        }

      }
    };
    //! \}

  }
}

#endif // DUNE_PDELAB_CONSTRAINTS_P0GHOST_HH
