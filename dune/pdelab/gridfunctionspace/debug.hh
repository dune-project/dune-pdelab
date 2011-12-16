// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_DEBUG_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_DEBUG_HH

#include <cstddef>
#include <ostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/static_assert.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/grid/common/gridenums.hh>

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

    struct DebugDofVectorFilterInterface {
      template<class Element>
      void bind(const Element &);

      template<class CV>
      void print(std::ostream &stream, std::size_t codim,
                 std::size_t subEntity, const CV &dofValues);
    };


    template<class Imp>
    class DebugDofVectorFilterDefaultPrint {
      Imp &asImp() { return static_cast<Imp&>(*this); }
      const Imp &asImp() const { return static_cast<const Imp&>(*this); }

    public:
      template<class CV>
      void print(std::ostream &stream, std::size_t codim,
                 std::size_t subEntity, const CV &dofValues)
      {
        if(asImp().doPrint(codim, subEntity, dofValues)) {
          stream << asImp().label(codim, subEntity, dofValues);
          for(std::size_t i = 0; i < dofValues.size(); ++i)
            stream << " " << dofValues[i];
          stream << "\n";
        }
      }

      template<class CV>
      void print(std::ostream &stream, std::size_t codim,
                 std::size_t subEntity, const CV &dofValues) const
      {
        if(asImp().doPrint(codim, subEntity, dofValues)) {
          stream << asImp().label();
          for(std::size_t i = 0; i < dofValues.size(); ++i)
            stream << " " << dofValues;
          stream << "\n";
        }
      }
    };

    template<class GV>
    class DebugDofVectorFilterUniqueInteriorBorder :
      public DebugDofVectorFilterDefaultPrint
        <DebugDofVectorFilterUniqueInteriorBorder<GV> >
    {
      typedef typename GV::Grid Grid;
      typedef typename Grid::LocalIdSet IdSet;
      typedef typename IdSet::IdType Id;
      typedef typename GV::template Codim<0>::EntityPointer EP;

      typedef typename GV::ctype DF;
      static const std::size_t dim = GV::dimension;

      const GV &gv;
      const IdSet &idSet;
      EP ep;

      std::vector<std::vector<PartitionType> > partitions;
      std::vector<std::vector<bool> > boundary;

      std::set<Id> alreadyPrinted;

      static inline PartitionType
      mergePartitionTypes(PartitionType a, PartitionType b) {
        dune_static_assert(InteriorEntity == 0 && BorderEntity == 1 &&
                           OverlapEntity == 2 && FrontEntity == 3 &&
                           GhostEntity == 4, "Bummer, our assumptions about "
                           "the PartitionType enumerators are no longer "
                           "valid");
        static const PartitionType pttable[5][5] = {
          { InteriorEntity, BorderEntity, BorderEntity,  FrontEntity,  FrontEntity  },
          { BorderEntity,   BorderEntity, BorderEntity,  BorderEntity, BorderEntity },
          { BorderEntity,   BorderEntity, OverlapEntity, FrontEntity,  FrontEntity  },
          { BorderEntity,   BorderEntity, FrontEntity,   FrontEntity,  FrontEntity  },
          { BorderEntity,   BorderEntity, FrontEntity,   FrontEntity,  GhostEntity  }
        };
        if(a >= 5 || a < 0)
          DUNE_THROW(Exception, "Unknown partition type " << a);
        if(b >= 5 || b < 0)
          DUNE_THROW(Exception, "Unknown partition type " << b);
        return pttable[a][b];
      }

    public:
      DebugDofVectorFilterUniqueInteriorBorder(const GV &gv_) :
        gv(gv_), idSet(gv.grid().localIdSet()),
        ep(gv.template begin<0>()),
        partitions(dim+1), boundary(dim+1)
      {
        partitions[0].resize(1);
        boundary[0].resize(1, false);
      }

      void bind(const typename GV::template Codim<0>::Entity &e) {
        ep = EP(e);

        const GenericReferenceElement<DF, dim> &refelem =
          GenericReferenceElements<DF, dim>::general(e.type());

        // partition for codim 0
        partitions[0][0] = e.partitionType();
        // boundary for codim 0 is always false, no deed to recompute

        // resize vectors and assign default values, as appropriate
        for(std::size_t c = 1; c <= GV::dimension; ++c) {
          std::size_t size = refelem.size(c);
          partitions[c].assign(size, partitions[0][0]);
          boundary[c].assign(size, false);
        }

        // partition and boundary for codim 1
        typedef typename GV::IntersectionIterator IIterator;
        const IIterator &iend = gv.iend(e);
        for(IIterator iit = gv.ibegin(e); iit != iend; ++iit) {
          std::size_t index = iit->indexInInside();
          if(iit->boundary())
            boundary[1][index] = true;
          if(iit->neighbor())
            partitions[1][index] =
              mergePartitionTypes(partitions[1][iit->indexInInside()],
                                  iit->outside()->partitionType());
        }

        // partition and boundary for higher codims
        // iterate over codim1-entities
        for(int i = 0; i < refelem.size(1); ++i)
          //iterate over codimensions
          for(std::size_t c = 2; c <= GV::dimension; ++c) {
            std::size_t size = refelem.size(i, 1, c);
            // iterate over higher-codimension entities within the
            // codim1-entity
            for(std::size_t s = 1; s < size; ++s) {
              std::size_t index = refelem.subEntity(i, 1, s, c);
              partitions[c][index] =
                mergePartitionTypes(partitions[c][index], partitions[1][i]);
            }
          }
      }

      template<class CV>
      bool doPrint(std::size_t codim, std::size_t subEntity,
                   const CV &dofValues)
      {
        if(dofValues.empty())
          return false;
        return alreadyPrinted.insert(idSet.subId(*ep, subEntity, codim))
          .second;
      }

      template<class CV>
      std::string label(std::size_t codim, std::size_t subEntity,
                        const CV &dofValues) const
      {
        const GenericReferenceElement<DF, dim> &refelem =
          GenericReferenceElements<DF, dim>::general(ep->type());

        std::ostringstream stream;
        stream << refelem.type(subEntity, codim) << "@("
               << ep->geometry().global(refelem.position(subEntity, codim))
               << ")[" << partitions[codim][subEntity];
        if(boundary[codim][subEntity])
          stream << ",boundary";
        stream << "]:";
        return stream.str();
      }
    };

    template<class GFS, class V>
    void debugDofVectorByPos(std::ostream &stream, const GFS &gfs, const V &v)
    {
      DebugDofVectorFilterUniqueInteriorBorder
        <typename GFS::Traits::GridViewType> filter(gfs.gridview());
      debugDofVector(stream, gfs, v, filter);
    }

    template<class GFS, class V, class Filter>
    void debugDofVector(std::ostream &stream, const GFS &gfs, const V &v,
                        Filter &filter)
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
        filter.bind(*eit);

        lfs.bind(*eit);
        lfs.vread(v, vl);

        std::vector<std::vector<std::vector<CF> > > perEntityDofs;
        TypeTree::applyToTree
          ( lfs, DebugDofVectorByPosVisitor<LV>(vl, perEntityDofs) );

        for(std::size_t c = 0; c < perEntityDofs.size(); ++c)
          for(std::size_t se = 0; se < perEntityDofs[c].size(); ++se)
            filter.print(stream, c, se, perEntityDofs[c][se]);
      }
    }

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_DEBUG_HH
