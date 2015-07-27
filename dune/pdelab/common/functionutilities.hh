// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_COMMON_FUNCTIONUTILITIES_HH
#define DUNE_PDELAB_COMMON_FUNCTIONUTILITIES_HH

#include <limits>
#include <ostream>
#include <memory>

#include <dune/common/debugstream.hh>

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/type.hh>

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/utility/hierarchicsearch.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup PDELab_Function Function
    //! \ingroup PDELab
    //! \{

    //! Integrate a GridFunction
    /**
     * \code
#include <dune/pdelab/common/functionutilities.hh>
     * \endcode
     *
     * Integrate a GridFunction over the domain given by the GridFunction's
     * GridView.  In the parallel case, this function integrates over the
     * Interior_Partition only.  If the accumulated result over all processors
     * result is required, use something like
     * \code
integrateGridFunction(gf, sum);
sum = gf.getGridView().comm().sum(sum);
     * \endcode
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
      typedef typename GV::template Codim<0>::Geometry Geometry;
      typedef typename GF::Traits::RangeType Range;
      typedef typename GF::Traits::DomainFieldType DF;
      static const int dimD = GF::Traits::dimDomain;
      typedef Dune::QuadratureRule<DF,dimD> QR;
      typedef Dune::QuadratureRules<DF,dimD> QRs;

      sum = 0;
      Range val;
      for(const auto& cell : elements(gf.getGridView(),Dune::Partitions::interior)) {
        const Geometry& geo = cell.geometry();
        Dune::GeometryType gt = geo.type();
        const QR& rule = QRs::rule(gt,qorder);
        for (const auto& qip : rule) {
          // evaluate the given grid functions at integration point
          gf.evaluate(cell,qip.position(),val);

          // accumulate error
          val *= qip.weight() * geo.integrationElement(qip.position());
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
      typedef typename GV::template Codim<0>::Entity Entity;
      typedef typename GF::Traits::DomainType Domain;
      typedef typename GF::Traits::RangeType Range;

    public:
      //! Constructor
      /**
       * Construct the object and find the rank and the entity containing the
       * global coordinate \c xg.  The is currently not facility to recompute
       * this information later, e.g. after a loadBalance() on the grid.
       *
       * \param gf The GridFunction to probe, either as a reference, a
       *           pointer, or a shared_ptr.
       * \param xg The global coordinate the evaluate at.
       */
      template<class GFHandle>
      GridFunctionProbe(const GFHandle& gf, const Domain& xg)
      {
        setGridFunction(gf);
        xl = 0;
        evalRank = gfp->getGridView().comm().size();
        int myRank = gfp->getGridView().comm().rank();
        try {
          e.reset(new Entity
                  (HierarchicSearch<typename GV::Grid, typename GV::IndexSet>
                   (gfp->getGridView().grid(), gfp->getGridView().indexSet()).
                   template findEntity<Interior_Partition>(xg)));
          // make sure only interior entities are accepted
          if(e->partitionType() == InteriorEntity)
            evalRank = myRank;
        }
        catch(const Dune::GridError&) { /* do nothing */ }
        evalRank = gfp->getGridView().comm().min(evalRank);
        if(myRank == evalRank)
          xl = e->geometry().local(xg);
        else
          e.reset();
        if(myRank == 0 && evalRank == gfp->getGridView().comm().size())
          dwarn << "Warning: GridFunctionProbe at (" << xg << ") is outside "
                << "the grid" << std::endl;
      }

      //! Set a new GridFunction
      /**
       * This takes the GridFunction as a refence.  The referenced object must
       * be valid for as long a the GridFunctionProbe is evaluated, or until
       * setGridFunction() is called again.
       */
      void setGridFunction(const GF &gf) {
        gfsp.reset();
        gfp = &gf;
      }

      //! Set a new GridFunction
      /**
       * This takes the GridFunction as a pointer.  Ownership of the
       * GridFunction object is transferred to the GridFunctionProbe.  The
       * GridFunction object will be deleted by the GridFunctionProbe on
       * destruction, or the next time setGridFunction is called.
       */
      void setGridFunction(const GF *gf) {
        gfsp.reset(gf);
        gfp = gf;
      }

      //! Set a new GridFunction
      /**
       * This takes the GridFunction as a shared_ptr.  Ownership of the
       * GridFunction object is shared with other places in the program, that
       * hold a shared_ptr to the same GridFunction object.  The GridFunction
       * object will be deleted by the GridFunctionProbe on destruction, or
       * the next time setGridFunction is called, if this GridFunctionProbe
       * holds the last reference at that time.
       */
      void setGridFunction(const std::shared_ptr<const GF> &gf) {
        gfsp = gf;
        gfp = &*gf;
      }

      //! evaluate the GridFunction and broadcast result to all ranks
      /**
       * \param val Store the result here.
       *
       * \note If the GridFunctionProbe is outside the grid NaN will be stored
       *       in val.
       */
      void eval_all(Range& val) const {
        typedef typename GF::Traits::RangeFieldType RF;
        if(evalRank == gfp->getGridView().comm().size())
          val = std::numeric_limits<RF>::quiet_NaN();
        else {
          if(gfp->getGridView().comm().rank() == evalRank)
            gfp->evaluate(*e, xl, val);
          gfp->getGridView().comm().broadcast(&val,1,evalRank);
        }
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
       * \note If the GridFunctionProbe is outside the grid NaN will be stored
       *       in val.
       */
      void eval(Range& val, int rank = 0) const {
        eval_all(val);
      }

    private:
      std::shared_ptr<const GF> gfsp;
      const GF *gfp;
      std::shared_ptr<Entity> e;
      Domain xl;
      int evalRank;
    };

    //! \} Function

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_COMMON_FUNCTIONUTILITIES_HH
