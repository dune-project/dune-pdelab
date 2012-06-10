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
        \tparam J Neumann stress boundary function (vector- or scalar-valued).
                  Scalar values will be interpreted as the magnitude of a vector
                  oriented in outer normal direction.
                  For prescribed pressure you can use $J=p \cdot n$.
        \tparam IP A class providing the interior penalty factor for each face
    */
    template<typename GV, typename RF, typename F, typename B, typename V, typename J,
             typename IP = DefaultInteriorPenalty<typename V::Traits::RangeFieldType> >
    class StokesDGParameters
      : public NavierStokesDefaultParameters<GV,RF,F,B,V,J>
    {
    private:

      typedef NavierStokesDefaultParameters<GV,RF,F,B,V,J> Base;

      void initFromString(const std::string & method)
      {
        std::string s = method;
        std::transform(s.begin(), s.end(), s.begin(), tolower);

        // nipg (epsilon=1) 2d p1 -> Klaus sagt sollte auch sigma 1 klappen
        if (s.find("nipg") != std::string::npos)
          {
            _epsilon = 1;
            return;
          }
        // sipg (epsilon=-1) 2d p1 -> Klaus sagt sigma=3.9irgendwas
        if (s.find("sipg") != std::string::npos)
          {
            _epsilon = -1;
            return;
          }
        // obb sigma = 0, epsilon =
        if (s == "obb")
          {
            _epsilon = 1;
            return;
          }
        // extract parameters
        {
          double sigma, beta;
          if (3 == sscanf(s.c_str(), "%d %lg %lg", &_epsilon, &sigma, &beta))
            return;
        }
        DUNE_THROW(Dune::Exception, "Unknown DG type " << method);
      }

    protected:
      StokesDGParameters(const std::string & method, const RF mu, const RF rho,
                         F & f, B & b, V & v, J & j, IP & ip)
        : Base(mu,rho,f,b,v,j),
          _ip(ip)
      {
        initFromString(method);
      }

      
    public:

      //! Traits class
      typedef typename Base::Traits Traits;

      StokesDGParameters(const std::string & method, const RF mu,
                         F & f, B & b, V & v, J & j, IP & ip)
        : Base(mu,1.0,f,b,v,j),
          _ip(ip)
      {
        initFromString(method);
      }

      StokesDGParameters(const Dune::ParameterTree & configuration,
                         F & f, B & b, V & v, J & j, IP & ip)
        : Base(configuration,f,b,v,j),
          _ip(ip),
          _epsilon(configuration.get<int>("epsilon"))
      {}


      //! Rescaling factor for the incompressibility equation
      typename Traits::RangeField
      incompressibilityScaling ( typename Traits::RangeField  dt ) const
      {
        typename J::Traits::Range y(1.0 / dt);
        return y;
      }

      /** \brief Interior penalty parameter.

          \return The coefficient of the interior penalty term.
      */
      template<typename I>
      typename Traits::RangeField
      getFaceIP(const I & ig) const
      {
        return _ip.getFaceIP(ig);
      }

      /** \brief Interior penalty parameter.

          \return The coefficient of the interior penalty term.
      */
      template<typename I>
      typename Traits::RangeField
      getFaceIP(const I & ig, const typename Traits::IntersectionDomain& ) const
      {
        return _ip.getFaceIP(ig);
      }

      //! Return the symmetry factor epsilon for this IP
      //! discretization
      int
      epsilonIPSymmetryFactor()
      {
        return _epsilon;
      }

    private:

      IP & _ip;              // Interior penalty
      int _epsilon;          // IP symmetry factor
    };


    /** \brief Parameter class for local operator NavierStokesDG

        \tparam F Momentum source term function
        \tparam B Boundary condition function
        \tparam V Dirichlet velocity boundary condition function
        \tparam J Neumann stress boundary function (vector- or scalar-valued).
                  Scalar values will be interpreted as the magnitude of a vector
                  oriented in outer normal direction.
                  For prescribed pressure you can use $J=p \cdot n$.
        \tparam IP A class providing the interior penalty factor for each face
    */
    template<typename GV, typename RF, typename F, typename B, typename V, typename J,
             typename IP = DefaultInteriorPenalty<typename V::Traits::RangeField> >
    class NavierStokesDGParameters
      : public StokesDGParameters<GV,RF,F,B,V,J,IP>
    {

      //!  of base class
      typedef StokesDGParameters<GV,F,B,V,J,IP> Base;

    public:

      NavierStokesDGParameters(const std::string & method, const RF mu, const RF rho,
                               F & f, B & b, V & v, J & j, IP & ip)
        : Base(method,mu,rho,f,b,v,j,ip)
      {}

      NavierStokesDGParameters(const Dune::ParameterTree & configuration,
                               F & f, B & b, V & v, J & j, IP & ip)
        : Base(configuration,f,b,v,j,ip)
      {}

    };


    namespace NavierStokesDGImp{
      /** \brief Compile-time switch for the boundary slip factor

          If the parameter class declares
          static const bool enable_variable_slip = true;
          it may also define a function

          template<typename IntersectionGeometry>
          typename Traits::RangeField
          boundarySlip (const IntersectionGeometry& ig,
          const typename Traits::IntersectionDomain& x) const

          which returns a value in [0,1] and may be used to define a
          smooth transition from a slip to a no-slip boundary.

          @{
      */
      template< typename PRM, typename Dummy = void>
      struct VariableBoundarySlipSwitch{

        template<typename IntersectionGeometry>
        static typename PRM::Traits::RangeField
        boundarySlip
        (const PRM & ,
         const IntersectionGeometry& ,
         const typename PRM::Traits::IntersectionDomain& )
        {
          return 1.0;
        }

      };

      template< typename PRM>
      struct VariableBoundarySlipSwitch
      <PRM,typename Dune::enable_if<PRM::enable_variable_slip>::type>
      {
        template<typename IntersectionGeometry>
        static typename PRM::Traits::RangeField
        boundarySlip
        (const PRM & prm,
         const IntersectionGeometry& ig,
         const typename PRM::Traits::IntersectionDomain& x)
        {
          return prm.boundarySlip(ig,x);
        }

      };
      /** @} */
    };


  }
}

#endif
