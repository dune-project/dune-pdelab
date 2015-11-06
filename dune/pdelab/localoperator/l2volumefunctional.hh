// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_LOCALOPERATOR_L2VOLUMEFUNCTIONAL_HH
#define DUNE_PDELAB_LOCALOPERATOR_L2VOLUMEFUNCTIONAL_HH

#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dune/localfunctions/common/interfaceswitch.hh>

#include <dune/pdelab/localoperator/idefault.hh>
#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/localoperator/flags.hh>

namespace Dune {
  namespace PDELab {
    //! \addtogroup LocalOperator
    //! \ingroup PDELab
    //! \{

    /** \brief A local operator that tests a function against a test function space defined on the entire grid
     *
     * It computes
     * \f{align*}{
     *           b_i = \int_\Omega f \varphi_i\, dx
     * \f}
     * with conforming finite elements on all types of grids in any dimension
     * \tparam F grid function type giving f
     */
    template<typename F>
    class L2VolumeFunctional
    : public LocalOperatorDefaultFlags
    {
    public:
      // residual assembly flags
      enum { doLambdaVolume = true };

      /** \brief Constructor
       *
       * \param f Function f to integrate over
       * \param quadOrder Order of the quadrature rule used for integrating over the element
       */
      L2VolumeFunctional (const F& f, unsigned int quadOrder)
        : f_(f), quadOrder_(quadOrder)
      {}

      // volume integral depending only on test functions
      template<typename EG, typename LFSV, typename R>
      void lambda_volume (const EG& eg, const LFSV& lfsv, R& r) const
      {
        // domain and range field type
        typedef FiniteElementInterfaceSwitch<typename LFSV::Traits::FiniteElementType> FESwitch;
        typedef BasisInterfaceSwitch<typename FESwitch::Basis> BasisSwitch;

        typedef typename BasisSwitch::DomainField DF;
        typedef typename BasisSwitch::RangeField RF;
        typedef typename BasisSwitch::Range Range;

        // dimensions
        static const int dimLocal = EG::Geometry::mydimension;

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const Dune::QuadratureRule<DF,dimLocal>& rule =
          Dune::QuadratureRules<DF,dimLocal>::rule(gt,quadOrder_);

        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dimLocal>::const_iterator it =
               rule.begin(); it!=rule.end(); ++it)
          {
            // evaluate shape functions
            std::vector<Range> phi(lfsv.size());
            FESwitch::basis(lfsv.finiteElement()).
              evaluateFunction(it->position(),phi);

            // evaluate right hand side parameter function
            typename F::Traits::RangeType y(0.0);
            f_.evaluate(eg.entity(),it->position(),y);

            // integrate f
            RF factor = r.weight() * it->weight() * eg.geometry().integrationElement(it->position());
            for (size_t i=0; i<lfsv.size(); i++)
              r.rawAccumulate(lfsv,i,y*phi[i]*factor);
          }
      }

    protected:
      const F& f_;

      // Quadrature rule order
      unsigned int quadOrder_;
    };

    //! \} group LocalOperator
  } // namespace PDELab
} // namespace Dune

#endif
