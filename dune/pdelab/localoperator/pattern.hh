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


      template<typename Context>
      void volumePattern(Context& ctx) const
      {
        if (ctx.fastDG())
        {
          ctx.pattern().addLink(ctx.test().space(),0,ctx.trial().space(),0);
          return;
        }
        for (auto i : ctx.test().space())
          for (auto j : ctx.trial().space())
            ctx.pattern().addLink(ctx.test().space(),i,ctx.trial().space(),j);
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

      template<typename Context>
      void skeletonPattern(Context& ctx) const
      {
        auto& inside  = ctx.inside();
        auto& outside = ctx.outside();

        if (ctx.fastDG())
        {
          ctx.pattern(inside,outside).addLink(inside.test().space(),0,outside.trial().space(),0);
          ctx.pattern(outside,inside).addLink(outside.test().space(),0,inside.trial().space(),0);
          return;
        }

        auto& pattern_io = ctx.pattern(inside,outside);
        for (auto i : inside.test().space())
          for (auto j : outside.trial().space())
            pattern_io.addLink(inside.test().space(),i,outside.trial().space(),j);

        auto& pattern_oi = ctx.pattern(outside,inside);
        for (auto i : outside.test().space())
          for (auto j : inside.trial().space())
            pattern_oi.addLink(outside.test().space(),i,inside.trial().space(),j);
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

      template<typename Context>
      void boundaryPattern(Context& ctx) const
      {
        auto& inside  = ctx.inside();

        if (ctx.fastDG())
        {
          ctx.pattern(inside,inside).addLink(inside.test().space(),0,inside.trial().space(),0);
          return;
        }

        auto& pattern = ctx.pattern(inside,inside);
        for (auto i : inside.test().space())
          for (auto j : inside.trial().space())
            pattern.addLink(inside.test().space(),i,inside.trial().space(),j);

      }
   };

    //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_LOCALOPERATOR_PATTERN_HH
