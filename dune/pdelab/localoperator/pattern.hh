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

      // define sparsity pattern of operator representation
      template<typename LFSU, typename LFSV>
      void pattern_volume (const LFSU& lfsu, const LFSV& lfsv, 
                           LocalSparsityPattern& pattern) const
      {
        for (size_t i=0; i<lfsv.size(); ++i)
          for (size_t j=0; j<lfsu.size(); ++j)
            pattern.push_back(SparsityLink(i,j));
      }
   };

    //! sparsity pattern generator
    class FullSkeletonPattern
    {
    public:

      // define sparsity pattern connecting self and neighbor dofs
      template<typename LFSU, typename LFSV>
      void pattern_skeleton (const LFSU& lfsu_s, const LFSV& lfsv_s, const LFSU& lfsu_n, const LFSV& lfsv_n, 
                             LocalSparsityPattern& pattern_sn, LocalSparsityPattern& pattern_ns) const
      {
        for (int i=0; i<lfsv_s.size(); ++i)
          for (int j=0; j<lfsu_n.size(); ++j)
            pattern_sn.push_back(SparsityLink(i,j));
        // other half of the pattern made by the neighboring element
      }
   };

    //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif
