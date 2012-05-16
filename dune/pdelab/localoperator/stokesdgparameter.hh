#ifndef DUNE_PDELAB_LOCALOPERATOR_STOKESDGPARAMETER_HH
#define DUNE_PDELAB_LOCALOPERATOR_STOKESDGPARAMETER_HH

#include <dune/common/parametertreeparser.hh>

#include <dune/pdelab/common/geometrywrapper.hh>

#include "dgparameter.hh"
#include "stokesparameter.hh"

namespace Dune {
  namespace PDELab {

    /** \brief Parameter class for local operator StokesDG

        \tparam F Momentum source term function
        \tparam B Boundary condition function
        \tparam V Dirichlet velocity boundary condition function
        \tparam P Dirichlet pressure boundary condition function
        \tparam IP A class providing the interior penalty factor for each face
    */
    template<typename GV, typename RF_, typename F, typename B, typename V, typename P,
             typename IP = DefaultInteriorPenalty<typename V::Traits::RangeFieldType> >
    class NavierStokesDGParameters
      : public NavierStokesParameters<GV,RF,F,B,V,P>
    {
    private:

      typedef NavierStokesParameters<GV,RF,F,B,V,P> Base;

      void initFromString(const std::string & method)
      {
        std::string s = method;
        std::transform(s.begin(), s.end(), s.begin(), tolower);

        // nipg (epsilon=1) 2d p1 -> Klaus sagt sollte auch sigma 1 klappen
        if (s.find("nipg") != std::string::npos)
          {
            epsilon = 1;
            return;
          }
        // sipg (epsilon=-1) 2d p1 -> Klaus sagt sigma=3.9irgendwas
        if (s.find("sipg") != std::string::npos)
          {
            epsilon = -1;
            return;
          }
        // obb sigma = 0, epsilon =
        if (s == "obb")
          {
            epsilon = 1;
            return;
          }
        // extract parameters
        {
          double sigma, beta;
          if (3 == sscanf(s.c_str(), "%d %lg %lg", &epsilon, &sigma, &beta))
            return;
        }
        DUNE_THROW(Dune::Exception, "Unknown DG type " << method);
      }


    public:

      //! Traits class
      typedef typename Base::Traits Traits;

      //! Common range field type
      typedef typename Traits::RangeFieldType RF;

      StokesDGParameters(const std::string & method, const RF mu_,
                         F & ff_, B & bf_, V & vf_, P & pf_, IP & ip_)
        : Base(mu,rho,f,b,v,j)
        , ip(ip_)
      {
        initFromString(method);
      }

      StokesDGParameters(const Dune::ParameterTree & configuration,
                         F & ff_, B & bf_, V & vf_, P & pf_, IP & ip_)
        : Base(configuration,f,b,v,j)
        , ip(ip_)
        , epsilon(configuration.get<int>("epsilon"))
      {}


      //! Rescaling factor for the incompressibility equation
      typename Traits::RangeFieldType
      incompressibilityScaling ( typename Traits::RangeFieldType  dt ) const
      {
        typename P::Traits::RangeType y(1.0 / dt);
        return y;
      }

      /** \brief Interior penalty parameter.

          \return The coefficient of the interior penalty term.
      */
      template<typename I>
      typename Traits::RangeFieldType
      getFaceIP(const I & ig) const
      {
        return ip.getFaceIP(ig);
      }

      /** \brief Interior penalty parameter.

          \return The coefficient of the interior penalty term.
      */
      template<typename I>
      typename Traits::RangeFieldType
      getFaceIP(const I & ig, const typename Traits::IntersectionDomainType& ) const
      {
        return ip.getFaceIP(ig);
      }

      //! Return the symmetry factor epsilon for this IP
      //! discretization
      int
      epsilonIPSymmetryFactor()
      {
        return epsilon;
      }

    private:

      IP & ip;              // Interior penalty
      int epsilon;          // IP symmetry factor
    };


    /** \brief Parameter class for local operator NavierStokesDG

        \tparam F Momentum source term function
        \tparam B Boundary condition function
        \tparam V Dirichlet velocity boundary condition function
        \tparam P Dirichlet pressure boundary condition function
        \tparam IP A class providing the interior penalty factor for each face
    */
    template<typename GV, typename RF, typename F, typename B, typename V, typename J,
             typename IP = DefaultInteriorPenalty<typename V::Traits::RangeFieldType> >
    class NavierStokesDGParameters
      : public StokesDGParameters<GV,RF,F,B,V,J,IP>
    {

      //! Type of base class
      typedef StokesDGParameters<GV,F,B,V,P,IP> Base;

    public:

      NavierStokesDGParameters(const std::string & method, const RF mu_,
                         F & ff_, B & bf_, V & vf_, P & pf_, IP & ip_)
        : Base(method,mu,rho,f,b,v,j,ip)
      {}

      StokesDGParameters(const Dune::ParameterTree & configuration,
                         F & ff_, B & bf_, V & vf_, P & pf_, IP & ip_)
        : Base(configuration,f,b,v,j,ip)
      {}

    };


    namespace NavierStokesDGImp{
      /** \brief Compile-time switch for the boundary slip factor

          If the parameter class declares
          static const bool enable_variable_slip = true;
          it may also define a function

          template<typename IntersectionGeometry>
          typename Traits::RangeFieldType
          boundarySlip (const IntersectionGeometry& ig,
          const typename Traits::IntersectionDomainType& x) const

          which returns a value in [0,1] and may be used to define a
          smooth transition from a slip to a no-slip boundary.

          @{
      */
      template< typename PRM, typename Dummy = void>
      struct VariableBoundarySlipSwitch{

        template<typename IntersectionGeometry>
        static typename PRM::Traits::RangeFieldType
        boundarySlip
        (const PRM & ,
         const IntersectionGeometry& ,
         const typename PRM::Traits::IntersectionDomainType& )
        {
          return 1.0;
        }

      };

      template< typename PRM>
      struct VariableBoundarySlipSwitch
      <PRM,typename Dune::enable_if<PRM::enable_variable_slip>::type>
      {
        template<typename IntersectionGeometry>
        static typename PRM::Traits::RangeFieldType
        boundarySlip
        (const PRM & prm,
         const IntersectionGeometry& ig,
         const typename PRM::Traits::IntersectionDomainType& x)
        {
          return prm.boundarySlip(ig,x);
        }

      };
      /** @} */
    };


  }
}

#endif
