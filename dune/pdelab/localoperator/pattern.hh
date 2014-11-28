// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_PATTERN_HH
#define DUNE_PDELAB_PATTERN_HH

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
        if (_block_optimization)
          {
            if (lfsu.size() > 0 && lfsv.size() > 0)
              pattern.addLink(lfsv,0,lfsu,0);
          }
        else
          {
            for (size_t i=0; i<lfsv.size(); ++i)
              for (size_t j=0; j<lfsu.size(); ++j)
                pattern.addLink(lfsv,i,lfsu,j);
          }
      }

      FullVolumePattern(bool block_optimization = false)
        : _block_optimization(block_optimization)
      {}

    private:
      const bool _block_optimization;

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
        if (_block_optimization)
          {
            if (lfsv_s.size() > 0 && lfsu_n.size() > 0)
              pattern_sn.addLink(lfsv_s,0,lfsu_n,0);

            if (lfsv_n.size() > 0 && lfsu_s.size() > 0)
              pattern_ns.addLink(lfsv_n,0,lfsu_s,0);
          }
        else
          {
            for (unsigned int i=0; i<lfsv_s.size(); ++i)
              for (unsigned int j=0; j<lfsu_n.size(); ++j)
                pattern_sn.addLink(lfsv_s,i,lfsu_n,j);

            for (unsigned int i=0; i<lfsv_n.size(); ++i)
              for (unsigned int j=0; j<lfsu_s.size(); ++j)
                pattern_ns.addLink(lfsv_n,i,lfsu_s,j);
          }
      }

      FullSkeletonPattern(bool block_optimization = false)
        : _block_optimization(block_optimization)
      {}

    private:
      const bool _block_optimization;

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
        if (_block_optimization)
          {
            if (lfsu_s.size() > 0 && lfsv_s.size() > 0)
              pattern_ss.addLink(lfsv_s,0,lfsu_s,0);
          }
        else
          {
            for (unsigned int i=0; i<lfsv_s.size(); ++i)
              for (unsigned int j=0; j<lfsu_s.size(); ++j)
                pattern_ss.addLink(lfsv_s,i,lfsu_s,j);
          }
      }

      FullBoundaryPattern(bool block_optimization = false)
        : _block_optimization(block_optimization)
      {}

    private:
      const bool _block_optimization;

   };

    //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif
