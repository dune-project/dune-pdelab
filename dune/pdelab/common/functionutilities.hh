// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_COMMON_FUNCTIONUTILITIES_HH
#define DUNE_PDELAB_COMMON_FUNCTIONUTILITIES_HH

#include <memory>

#include <dune/common/geometrytype.hh>

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/quadraturerules.hh>
#include <dune/grid/utility/hierarchicsearch.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup PDELab_Function Function
    //! \ingroup PDELab
    //! \{

    //! Integrate a GridFunction
    /**
     * Integrate a GridFunction over the domain given by the GridFunction's
     * GridView.
     *
     * \tparam GF Type of the GridFunction.
     * \param gf     The GridFunction object.
     * \param sum    Resulting integral.  There is no need to clear this
     *               variable before calling this function.
     * \param qorder Quadrature order to use.  If the GridFunction is
     *               element-wise polynomial, then this is the order of the
     *               highest-order monom needed to represent the function.
     */
    template<typename GF>
    void integrateGridFunction(const GF& gf,
                               typename GF::Traits::RangeType& sum,
                               unsigned qorder = 1) {
      typedef typename GF::Traits::GridViewType GV;
      typedef typename GV::template Codim<0>::Iterator EIterator;
      typedef typename GV::template Codim<0>::Geometry Geometry;
      typedef typename GF::Traits::RangeType Range;
      typedef typename GF::Traits::DomainFieldType DF;
      static const int dimD = GF::Traits::dimDomain;
      typedef Dune::QuadratureRule<DF,dimD> QR;
      typedef Dune::QuadratureRules<DF,dimD> QRs;
      typedef typename QR::const_iterator QIterator;

      sum = 0;
      Range val;
      const EIterator eend = gf.getGridView().template end<0>();
      for(EIterator eit = gf.getGridView().template begin<0>();
          eit != eend; ++eit) {
        const Geometry& geo = eit->geometry();
        Dune::GeometryType gt = geo.type();
        const QR& rule = QRs::rule(gt,qorder);
        const QIterator qend = rule.end();

        for (QIterator qit=rule.begin(); qit != qend; ++qit)
        {
          // evaluate the given grid functions at integration point
          gf.evaluate(*eit,qit->position(),val);

          // accumulate error
          val *= qit->weight() * geo.integrationElement(qit->position());
          sum += val;
        }
      }
    }

    //! Evaluate a GridFunction at a certain global coordinate
    /**
     * This class should work correctly even in a parallel setup.  We look for
     * the entity on containing the global coordinate on all ranks, then find
     * the rank that actually has the entity.  The GridFunction is then
     * evaluated only for that rank, and the result is communicated to all the
     * other ranks.
     *
     * \tparam GF Type of the GridFunction to evaluate.
     */
    template<typename GF>
    class GridFunctionProbe {
      typedef typename GF::Traits::GridViewType GV;
      typedef typename GV::template Codim<0>::EntityPointer EPtr;
      typedef typename GF::Traits::DomainType Domain;
      typedef typename GF::Traits::RangeType Range;

    public:
      //! Constructor
      /**
       * Construct the object and find the rank and the entity containing the
       * global coordinate \c xg.  The is currently not facility to recompute
       * this information later, e.g. after a loadBalance() on the grid.
       *
       * \param gf_ The GridFunction to probe.
       * \param xg  The global coordinate the evaluate at.
       */
      GridFunctionProbe(const GF& gf_, const Domain& xg)
        : gf(gf_), e(), xl(0), evalRank(gf.getGridView().comm().size())
      {
        int myRank = gf.getGridView().comm().rank();
        try {
          e.reset(new EPtr
                  (HierarchicSearch<typename GV::Grid, typename GV::IndexSet>
                   (gf.getGridView().grid(), gf.getGridView().indexSet()).
                   findEntity(xg)));
          // make sure only interior entities are accepted
          if((*e)->partitionType() == InteriorEntity)
            evalRank = myRank;
        }
        catch(const Dune::GridError&) { /* do nothing */ }
        evalRank = gf.getGridView().comm().min(evalRank);
        if(myRank == evalRank)
          xl = (*e)->geometry().local(xg);
        else
          e.reset();
      }

      //! evaluate the GridFunction and broadcast result to all ranks
      /**
       * \param val Store the result here.
       */
      void eval_all(Range& val) const {
        if(gf.getGridView().comm().rank() == evalRank)
          gf.evaluate(**e, xl, val);
        gf.getGridView().comm().broadcast(&val,1,evalRank);
      }

      //! evaluate the GridFunction and communicate result to the given rank
      /**
       * \param val  Store the result here.  This variable should be considered
       *             undefined on any rank other than the one given in \c rank
       *             after the function returns.
       * \param rank Rank of the process to communicate the result to.
       *
       * \note CollectiveCommunication does not provide any direct
       *       process-to-process communication, so currently this function is
       *       identical with eval_all(), with the \c rank parameter ignored.
       */
      void eval(Range& val, int rank = 0) const {
        eval_all(val);
      }

    private:
      const GF& gf;
      std::auto_ptr<EPtr> e;
      Domain xl;
      int evalRank;
    };

    //! \} Function

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_COMMON_FUNCTIONUTILITIES_HH
