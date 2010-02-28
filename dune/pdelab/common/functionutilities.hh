// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_COMMON_FUNCTIONUTILITIES_HH
#define DUNE_PDELAB_COMMON_FUNCTIONUTILITIES_HH

#include <dune/common/geometrytype.hh>

#include <dune/grid/common/quadraturerules.hh>

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

    //! \} Function

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_COMMON_FUNCTIONUTILITIES_HH
