// -*- tab-width: 2; indent-tabs-mode: nil -*-
// vi: set et ts=2 sw=2 sts=2:
#ifndef DUNE_PDELAB_LOCALOPERATOR_DGINTERIORPENALTYPARAMETER_HH
#define DUNE_PDELAB_LOCALOPERATOR_DGINTERIORPENALTYPARAMETER_HH

#include <dune/common/parametertreeparser.hh>

namespace Dune {
  namespace PDELab {

    /** \brief Default implementation of the interior penalty factor.

        It computes the face diameter from the formula derived by
        Paul Houston et al.
        There is an overpenalized version, OverPenalizedInteriorPenalty,
        which is implemented in a different class due to efficiency reasons.

        \warning The scaling with the viscosity (in the Navier-Stokes case)
                 is not included in here any more!!!
    */
    template<typename RF>
    class DefaultInteriorPenalty
    {
    private :
      RF sigma;
    public :
      DefaultInteriorPenalty(const Dune::ParameterTree& config)
      {
        sigma = config.get<RF>("sigma");
      }

      template<typename GEO, typename IGEO, typename OGEO>
      RF getFaceIP(const GEO& geo, const IGEO& igeo, const OGEO& ogeo) const
      {
        using std::min;
        // volume of face divided by volume of element => 1/h_F
        return sigma * geo.volume()/min(igeo.volume(), ogeo.volume());
      }

      template<typename GEO, typename IGEO>
      RF getFaceIP(const GEO& geo, const IGEO& igeo) const
      {
        // volume of face divided by volume of element => 1/h_F
        return sigma * geo.volume()/igeo.volume();
      }
    }; // end class DefaultInteriorPenalty

    /** \brief Implementation of overpenalized interior penalty.

        It computes the face diameter from the formula derived by
        Paul Houston et al. and uses overpenalization.

        \warning The scaling with the viscosity (in the Navier-Stokes case)
                 is not included in here any more!!!
    */
    template<typename RF>
    class OverPenalizedInteriorPenalty
    {
    private :
      RF sigma;
      RF beta;
    public :
      OverPenalizedInteriorPenalty(const Dune::ParameterTree& config)
      {
        sigma = config.get<RF>("sigma");
        beta = config.get<RF>("beta");
      }

      template<typename GEO, typename IGEO, typename OGEO>
      RF getFaceIP(const GEO& geo, const IGEO& igeo, const OGEO& ogeo) const
      {
        using std::pow;
        using std::min;
        // volume of face divided by volume of element => 1/h_F
        return sigma * pow(geo.volume()/min(igeo.volume(),ogeo.volume()), beta);
      }

      template<typename GEO, typename IGEO>
      RF getFaceIP(const GEO& geo, const IGEO& igeo) const
      {
        using std::pow;
        // volume of face divided by volume of element => 1/h_F
        return sigma * pow(geo.volume()/igeo.volume(), beta);
      }
    }; // end class OverPenalizedInteriorPenalty

  } // end namespace PDELab
} // end namespace Dune
#endif
