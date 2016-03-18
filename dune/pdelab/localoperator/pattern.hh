// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_LOCALOPERATOR_PATTERN_HH
#define DUNE_PDELAB_LOCALOPERATOR_PATTERN_HH

#include<dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

namespace Dune {
  namespace PDELab {

    //! sparsity pattern generator
    class FullVolumePattern
    {
    public:

      // define sparsity pattern of operator representation
      template<typename LFSU, typename LFSV, typename LocalPattern>
      void pattern_volume (const LFSU& lfsu, const LFSV& lfsv,
                           LocalPattern& pattern) const
      {
        for (size_t i=0; i<lfsv.size(); ++i)
          for (size_t j=0; j<lfsu.size(); ++j)
            pattern.addLink(lfsv,i,lfsu,j);
      }
   };

    //! sparsity pattern generator
    class FullSkeletonPattern
    {
    public:

      // define sparsity pattern connecting self and neighbor dofs
      template<typename LFSU, typename LFSV, typename LocalPattern>
      void pattern_skeleton (const LFSU& lfsu_s, const LFSV& lfsv_s, const LFSU& lfsu_n, const LFSV& lfsv_n,
                            LocalPattern& pattern_sn,
                            LocalPattern& pattern_ns) const
      {
        for (unsigned int i=0; i<lfsv_s.size(); ++i)
          for (unsigned int j=0; j<lfsu_n.size(); ++j)
            pattern_sn.addLink(lfsv_s,i,lfsu_n,j);

        for (unsigned int i=0; i<lfsv_n.size(); ++i)
          for (unsigned int j=0; j<lfsu_s.size(); ++j)
            pattern_ns.addLink(lfsv_n,i,lfsu_s,j);
      }
   };

    //! sparsity pattern generator
    class FullBoundaryPattern
    {
    public:

      // define sparsity pattern connecting dofs on boundary elements
      template<typename LFSU, typename LFSV, typename LocalPattern>
      void pattern_boundary(const LFSU& lfsu_s, const LFSV& lfsv_s,
                            LocalPattern& pattern_ss) const
      {
        for (unsigned int i=0; i<lfsv_s.size(); ++i)
          for (unsigned int j=0; j<lfsu_s.size(); ++j)
            pattern_ss.addLink(lfsv_s,i,lfsu_s,j);
      }
   };

    //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_LOCALOPERATOR_PATTERN_HH
