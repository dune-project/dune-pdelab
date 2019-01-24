// -*- tab-width: 2; indent-tabs-mode: nil -*-
// vi: set et ts=2 sw=2 sts=2:

#ifndef DUNE_PDELAB_INSTATIONARY_ONESTEPMETHODS_HH
#define DUNE_PDELAB_INSTATIONARY_ONESTEPMETHODS_HH

#include <cmath>

#include <dune/common/exceptions.hh>
#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>
#include <dune/common/typetraits.hh>

namespace Dune {
  namespace PDELab { namespace OneStep {

    /**
     *  @addtogroup OneStepMethod
     *  @{
     */
    //! Base parameter class for time stepping scheme parameters
    /**
     * \tparam R C++ type of the floating point parameters
     */
    template<typename R>
    class Method
    {
    public:

      using Real = R;
      using Index = int;

      using value_type = Real;
      using index_type = Index;


      /*! \brief Return true if method is implicit
       */
      bool implicit() const
      {
        return _implicit;
      }

      Index stages() const
      {
        return _weights.rows();
      }

      Real weight(Index current_stage, Index old_stage) const
      {
        checkInput(current_stage,old_stage);
        return _weights[current_stage - 1][old_stage];
      }

      Real timeDerivativeWeight(Index current_stage, Index old_stage) const
      {
        checkInput(current_stage,old_stage);
        return _time_derivative_weights[current_stage - 1][old_stage];
      }

      Real timeStepFraction(Index stage) const
      {
        assert(stage > 0);
        assert(stage <= stages());
        return _time_step_fractions[stage];
      }

      bool active(Index current_stage, Index old_stage) const
      {
        checkInput(current_stage,old_stage);
        using std::abs;
        return abs(_weights[current_stage][old_stage]) > _activity_limit;
      }

      bool timeDerivativeActive(Index current_stage, Index old_stage) const
      {
        checkInput(current_stage,old_stage);
        using std::abs;
        return abs(_time_derivative_weights[current_stage][old_stage]) > _activity_limit;
      }

      std::string name() const
      {
        return _name;
      }

      void setName(const std::string& name)
      {
        _name = name;
      }

      R activityLimit() const
      {
        return _activity_limit;
      }

      void setActivityLimit(R activity_limit)
      {
        _activity_limit = activity_limit;
      }

      Method(
        const std::string& name,
        bool implicit,
        DynamicMatrix<Real> time_derivative_weights,
        DynamicMatrix<Real> weights,
        DynamicVector<Real> time_step_fractions
        )
        : _implicit(implicit)
        , _weights(weights)
        , _time_derivative_weights(time_derivative_weights)
        , _time_step_fractions(time_step_fractions)
        , _name(name)
      {}

    private:

      void checkInput(Index current_stage, Index old_stage) const
      {
        assert(current_stage > 0);
        assert(current_stage <= stages());
        assert(old_stage >= 0);
        assert(old_stage <= current_stage);
      }

      bool _implicit;
      R _activity_limit = 1e-6;
      DynamicMatrix<Real> _time_derivative_weights;
      DynamicMatrix<Real> _weights;
      DynamicVector<Real> _time_step_fractions;
      std::string _name;

    };



    /**
     * \brief Parameters to turn the OneStepMethod into an
     * one step theta method.
     *
     * For theta=0 this parameter class can be used with the
     * ExplicitOneStepMethod
     * \tparam R C++ type of the floating point parameters
     */
    template<typename R>
    Method<R> oneStepTheta(R theta)
    {
      // FIXME: Make this work for non-standard types
      static_assert(std::is_floating_point_v<R>,"R must be a floating point type");
      using namespace std::literals;
      return Method<R>{
        "one step theta"s,
        theta > 0,
        {
          { -1.0, 1.0 }
        },
        {
          { 1.0 - theta, theta }
        },
        { 0.0, 1.0 }
      };
    }

    /**
     * \brief Parameters to turn the ExplicitOneStepMethod into an
     * explicit Euler method.
     *
     * \tparam R C++ type of the floating point parameters
     */
    template<typename R>
    Method<R> explicitEuler(R = 0.0)
    {
      auto method = oneStepTheta<R>(0.0);
      method.setName("explicit Euler");
      return method;
    }


    /**
     * \brief Parameters to turn the OneStepMethod into an
     * implicit Euler method.
     *
     * \tparam R C++ type of the floating point parameters
     */
    template<typename R>
    Method<R> implicitEuler(R = 0.0)
    {
      auto method = oneStepTheta<R>(1.0);
      method.setName("implicit Euler");
      return method;
    }



    /**
     * \brief Parameters to turn the ExplicitOneStepMethod into a
     * Heun scheme.
     *
     * \tparam R C++ type of the floating point parameters
     */
    template<typename R>
    Method<R> heun(R = 0.0)
    {
      // FIXME: Make this work for non-standard types
      static_assert(std::is_floating_point_v<R>,"R must be a floating point type");
      using namespace std::literals;
      return Method<R>{
        "Heun"s,
        false,
        {
          { -1.0,  1.0, 0.0 },
          { -0.5, -0.5, 1.0 }
        },
        {
          { 1.0, 0.0, 0.0 },
          { 0.0, 0.5, 0.0 }
        },
        { 0.0, 1.0, 1.0 }
      };
    }


    /**
     * \brief Parameters to turn the ExplicitOneStepMethod into a
     * third order strong stability preserving (SSP) scheme.
     *
     * \tparam R C++ type of the floating point parameters
     */
    template<typename R>
    Method<R> shu3(R = 0.0)
    {
      // FIXME: Make this work for non-standard types
      static_assert(std::is_floating_point_v<R>,"R must be a floating point type");
      using namespace std::literals;
      return Method<R>{
        "Shu's third order method"s,
        false,
        {
          { -1.0    ,  1.0 ,  0.0    , 0.0 },
          { -0.75   , -0.25,  1.0    , 0.0 },
          { -1.0/3.0,  0.0 , -2.0/3.0, 1.0 }
        },
        {
          { 1.0, 0.0 , 0.0    , 0.0 },
          { 0.0, 0.25, 0.0    , 0.0 },
          { 0.0, 0.0 , 2.0/3.0, 0.0 }
        },
        { 0.0, 1.0, 0.5, 1.0 }
      };
    }


    /**
     * \brief Parameters to turn the ExplicitOneStepMethod into a
     * classical fourth order Runge-Kutta method
     *
     * \tparam R C++ type of the floating point parameters
     */
    template<typename R>
    Method<R> rk4(R = 0.0)
    {
      // FIXME: Make this work for non-standard types
      static_assert(std::is_floating_point_v<R>,"R must be a floating point type");
      using namespace std::literals;
      return Method<R>{
        "RK4"s,
        false,
        {
          { -1.0, 1.0, 0.0, 0.0, 0.0 },
          { -1.0, 0.0, 1.0, 0.0, 0.0 },
          { -1.0, 0.0, 0.0, 1.0, 0.0 },
          { -1.0, 0.0, 0.0, 0.0, 1.0 }
        },
        {
          { 0.5    , 0.0    , 0.0    , 0.0    , 0.0 },
          { 0.0    , 0.5    , 0.0    , 0.0    , 0.0 },
          { 0.0    , 0.0    , 1.0    , 0.0    , 0.0 },
          { 1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0, 0.0 }
        },
        { 0.0, 0.5, 0.5, 1.0, 1.0 }
      };
    }


    /**
     * \brief Parameters to turn the OneStepMethod into an
     * Alexander scheme.
     *
     * \tparam R C++ type of the floating point parameters
     */
    template<typename R>
    Method<R> alexander2(R = 0.0)
    {
      // FIXME: Make this work for non-standard types
      static_assert(std::is_floating_point_v<R>,"R must be a floating point type");
      using namespace std::literals;

      using std::sqrt;
      R alpha = 1.0 - 0.5 * sqrt(2.0);

      return Method<R>{
        "Alexander (order 2)"s,
        true,
        {
          { -1.0, 1.0, 0.0 },
          { -1.0, 0.0, 1.0 }
        },
        {
          { 0.0, alpha    , 0.0   },
          { 0.0, 1.0-alpha, alpha }
        },
        { 0.0, alpha, 1.0 }
      };
    }


    /**
     * \brief Parameters to turn the OneStepMethod into a
     * fractional step theta scheme.
     *
     * \tparam R C++ type of the floating point parameters
     */
    template<typename R>
    Method<R> fractionalStep(R = 0.0)
    {
      // FIXME: Make this work for non-standard types
      static_assert(std::is_floating_point_v<R>,"R must be a floating point type");
      using namespace std::literals;

      using std::sqrt;
      R theta  = 1.0 - 0.5*sqrt(2.0);
      R thetap = 1.0 - 2.0*theta;
      R alpha  = 2.0 - sqrt(2.0);
      R beta   = 1.0 - alpha;

      return Method<R>{
        "Fractional step theta"s,
        true,
        {
          { -1.0,  1.0,  0.0, 0.0 },
          {  0.0, -1.0,  1.0, 0.0 },
          {  0.0,  0.0, -1.0, 1.0 }
        },
        {
          { beta*theta, alpha*theta , 0.0        , 0.0         },
          { 0.0       , alpha*thetap, alpha*theta, 0.0         },
          { 0.0       , 0.0         , beta*theta , alpha*theta }
        },
        { 0.0, theta, 1.0-theta, 1.0 }
      };
    }


    /**

     * \brief Parameters to turn the OneStepMethod into an
     * Alexander3 scheme.
     *
     * \tparam R C++ type of the floating point parameters
     */
    template<typename R>
    Method<R> alexander3(R = 0.0)
    {
      // FIXME: Make this work for non-standard types
      static_assert(std::is_floating_point_v<R>,"R must be a floating point type");
      using namespace std::literals;

      R alpha = 0.4358665215;

      // Newton iteration for alpha
      for (int i=1; i<=10; i++)
      {
        alpha = alpha - (alpha*(alpha*alpha-3.0*(alpha-0.5))-1.0/6.0)/(3.0*alpha*(alpha-2.0)+1.5);
      }

      R tau2 = (1.0+alpha)*0.5;
      R b1 = -(6.0*alpha*alpha -16.0*alpha + 1.0)*0.25;
      R b2 = (6*alpha*alpha - 20.0*alpha + 5.0)*0.25;

      return Method<R>{
        "Alexander (claims order 3)"s,
        true,
        {
          { -1.0, 1.0, 0.0, 0.0 },
          { -1.0, 0.0, 1.0, 0.0 },
          { -1.0, 0.0, 0.0, 1.0 }
        },
        {
          { 0.0, alpha     , 0.0  , 0.0   },
          { 0.0, tau2-alpha, alpha, 0.0   },
          { 0.0, b1        , b2   , alpha }
        },
        { 0.0, alpha, tau2, 1.0 }
      };
    }

  } // end namespace OneStep
  } // end namespace PDELab
} // end namespace Dune
#endif // DUNE_PDELAB_INSTATIONARY_ONESTEPMETHODS_HH
