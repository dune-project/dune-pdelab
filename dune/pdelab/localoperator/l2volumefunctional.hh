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

#include<dune/pdelab/common/quadraturerules.hh>
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
        // Define types
        using FESwitch = FiniteElementInterfaceSwitch<typename LFSV::Traits::FiniteElementType>;
        using BasisSwitch = BasisInterfaceSwitch<typename FESwitch::Basis>;
        using Range = typename BasisSwitch::Range;

        // Reference to cell
        const auto& cell = eg.entity();

        // get geometry
        auto geo = eg.geometry();

        // Initialize vectors outside for loop
        std::vector<Range> phi(lfsv.size());
        typename F::Traits::RangeType y(0.0);

        // loop over quadrature points
        for (const auto& ip : quadratureRule(geo,quadOrder_))
          {
            // evaluate shape functions
            FESwitch::basis(lfsv.finiteElement()).
              evaluateFunction(ip.position(),phi);

            // evaluate right hand side parameter function
            f_.evaluate(cell,ip.position(),y);

            // integrate f
            auto factor = r.weight() * ip.weight() * geo.integrationElement(ip.position());
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

#endif // DUNE_PDELAB_LOCALOPERATOR_L2VOLUMEFUNCTIONAL_HH
