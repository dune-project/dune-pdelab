// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_PATTERN_HH
#define DUNE_PDELAB_PATTERN_HH

#include<dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/static_assert.hh>

#include"../common/geometrywrapper.hh"
#include"../gridoperatorspace/gridoperatorspace.hh"
#include"../gridoperatorspace/gridoperatorspaceutilities.hh"


namespace Dune {
  namespace PDELab {

    //! sparsity pattern generator
    class FullVolumePattern
    {
    public:
      enum { doPatternBoundary = false }; // is this really needed?
      enum { doPatternSkeleton = false };
      enum { doPatternVolume = true };

      // define sparsity pattern of operator representation
      template<typename EG, typename LFSU, typename LFSV>
      void pattern_volume (const EG& eg, const LFSU& lfsu, const LFSV& lfsv, 
                           LocalSparsityPattern& pattern) const
      {
        for (int i=0; i<lfsv.size(); ++i)
          for (int j=0; j<lfsu.size(); ++j)
            pattern.push_back(SparsityLink(i,j));
      }
   };

    //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif
