// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_DEBUG_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_DEBUG_HH

#include <cstddef>
#include <ostream>
#include <vector>

#include <dune/grid/common/genericreferenceelements.hh>

#include <dune/localfunctions/common/interfaceswitch.hh>
#include <dune/localfunctions/common/localkey.hh>

#include <dune/pdelab/common/typetree/traversal.hh>
#include <dune/pdelab/common/typetree/visitor.hh>
#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/localvector.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup GridFunctionSpace grid function space
    //! \ingroup PDELab
    //! \{

    //=======================================
    // functions for producing debugging output
    //=======================================

    template<class LV>
    class DebugDofVectorByPosVisitor :
      public TypeTree::DefaultVisitor,
      public TypeTree::DynamicTraversal,
      public TypeTree::VisitTree
    {
      typedef typename LV::value_type CF;
      const LV &lv;
      std::vector<std::vector<std::vector<CF> > > &perEntityDofs;

      template<class T>
      static T& access(std::vector<T> &v, std::size_t i) {
        if(i >= v.size())
          v.resize(i+1);
        return v[i];
      }
    public:
      DebugDofVectorByPosVisitor
      ( const LV &lv_,
        std::vector<std::vector<std::vector<CF> > > &perEntityDofs_
      ) :
        lv(lv_), perEntityDofs(perEntityDofs_)
      { perEntityDofs.clear(); }

      template<typename LFS, typename TreePath>
      void leaf(const LFS &lfs, const TreePath &) const {
        typedef typename LFS::Traits::FiniteElementType FE;
        const FE &fe = lfs.finiteElement();

        typedef FiniteElementInterfaceSwitch<FE> FESwitch;
        typedef typename FESwitch::Coefficients Coeffs;
        const Coeffs &coeffs = FESwitch::coefficients(fe);

        std::vector<std::vector<std::vector<CF> > > nodeDofs;
        for(std::size_t i = 0; i < coeffs.size(); ++i) {
          const LocalKey &lk = coeffs.localKey(i);
          access(access(access(nodeDofs, lk.codim()), lk.subEntity()),
                 lk.index())
            = lv(lfs, i);
        }
        if(nodeDofs.size() > perEntityDofs.size())
          perEntityDofs.resize(nodeDofs.size());
        for(std::size_t c = 0; c < nodeDofs.size(); ++c) {
          if(nodeDofs[c].size() > perEntityDofs[c].size())
            perEntityDofs[c].resize(nodeDofs[c].size());
          for(std::size_t e = 0; e < nodeDofs[c].size(); ++e)
            perEntityDofs[c][e].insert(perEntityDofs[c][e].end(),
                                       nodeDofs[c][e].begin(),
                                       nodeDofs[c][e].end());
        }
      }
    };

    template<class GFS, class V>
    void debugDofVectorByPos(std::ostream &stream, const GFS &gfs, const V &v)
    {
      typedef typename GFS::Traits::GridViewType GV;
      const GV &gv = gfs.gridview();

      typedef LocalFunctionSpace<GFS> LFS;
      LFS lfs(gfs);

      typedef typename GV::template Codim<0>::Iterator Iterator;
      const Iterator &eend = gv.template end<0>();

      typedef typename V::ElementType CF;
      typedef LocalVector<CF> LV;
      LV vl(lfs.maxSize());

      for(Iterator eit = gv.template begin<0>(); eit != eend; ++eit) {
        lfs.bind(*eit);
        lfs.vread(v, vl);

        std::vector<std::vector<std::vector<CF> > > perEntityDofs;
        TypeTree::applyToTree
          ( lfs, DebugDofVectorByPosVisitor<LV>(vl, perEntityDofs) );

        typedef typename GV::ctype DF;
        static const std::size_t dim = GV::dimension;
        const GenericReferenceElement<DF, dim> &refelem =
          GenericReferenceElements<DF, dim>::general(eit->type());

        typedef typename GV::template Codim<0>::Geometry Geometry;
        const Geometry &geo = eit->geometry();

        for(std::size_t c = 0; c < perEntityDofs.size(); ++c)
          for(std::size_t se = 0; se < perEntityDofs[c].size(); ++se)
            if(perEntityDofs[c][se].size() > 0) {
              stream << refelem.type(se,c) << "@("
                     << geo.global(refelem.position(se,c)) << "):";
              for(std::size_t i = 0; i < perEntityDofs[c][se].size(); ++i)
                stream << " " << perEntityDofs[c][se][i];
              stream << "\n";
            }
      }
    }

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_DEBUG_HH
