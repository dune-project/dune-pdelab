// -*- tab-width: 8; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_LOCALOPERATOR_LINEARACOUSTICSPARAMETER_HH
#define DUNE_PDELAB_LOCALOPERATOR_LINEARACOUSTICSPARAMETER_HH

#include<vector>

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/geometry/type.hh>
#include<dune/geometry/referenceelements.hh>
#include<dune/geometry/quadraturerules.hh>
#include<dune/pdelab/common/geometrywrapper.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/functionutilities.hh>

namespace Dune {
  namespace PDELab {

    /** \brief Traits class for linear acoustics parameters
     *
     * A class supplying parameters to a linear-acoustics local
     * operator has to define a public traits class exporting the needed
     * types and constants.
     */
    template<typename GV, typename RF>
    struct LinearAcousticsParameterTraits
    {
      //! \brief the grid view
      typedef GV GridViewType;

      //! \brief Enum for domain dimension
      enum {
        //! \brief dimension of the domain
        dimDomain = GV::dimension
      };

      //! \brief Export type for domain field
      typedef typename GV::Grid::ctype DomainFieldType;

      //! \brief domain type
      typedef Dune::FieldVector<DomainFieldType,dimDomain> DomainType;

      //! \brief domain type
      typedef Dune::FieldVector<DomainFieldType,dimDomain-1> IntersectionDomainType;

      //! \brief Export type for range field
      typedef RF RangeFieldType;

      //! \brief range type
      typedef Dune::FieldVector<RF,GV::dimension> RangeType;

      //! \brief range type
      typedef Dune::FieldVector<RF,GV::dimension+1> StateType;

      //! grid types
      typedef typename GV::Traits::template Codim<0>::Entity ElementType;
      typedef typename GV::Intersection IntersectionType;
    };

    template<typename T>
    class LinearAcousticsInitialValueAdapter
      : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<typename T::Traits::GridViewType,
                                                                               typename T::Traits::RangeFieldType,
                                                                               T::Traits::dimDomain+1,Dune::FieldVector<typename T::Traits::RangeFieldType,T::Traits::dimDomain+1> >
                                              ,LinearAcousticsInitialValueAdapter<T> >
    {
    public:
      typedef Dune::PDELab::GridFunctionTraits<typename T::Traits::GridViewType,
                                               typename T::Traits::RangeFieldType,
                                               T::Traits::dimDomain+1,Dune::FieldVector<typename T::Traits::RangeFieldType,T::Traits::dimDomain+1> > Traits;

      //! constructor
      LinearAcousticsInitialValueAdapter (const typename Traits::GridViewType& g_, const T& t_) : g(g_), t(t_) {}

      //! \copydoc GridFunctionBase::evaluate()
      inline void evaluate (const typename Traits::ElementType& e,
                            const typename Traits::DomainType& x,
                            typename Traits::RangeType& y) const
      {
        y = t.u0(e,x);
      }

      inline const typename Traits::GridViewType& getGridView () const
      {
        return g;
      }

    private:
      typename Traits::GridViewType g;
      const T& t;
    };

  template<typename GV, typename RF>
  class LinearAcousticsModelProblem
  {
  public:
    typedef LinearAcousticsParameterTraits<GV,RF> Traits;

    LinearAcousticsModelProblem ()
      : pi(3.141592653589793238462643), time(0.0)
    {
    }

    //! speed of sound
    typename Traits::RangeFieldType
    c (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
    {
      return 340.0;
    }

    //! Dirichlet boundary condition value
    typename Traits::StateType
    g (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, const typename Traits::StateType& s) const
    {
      typename Traits::DomainType xglobal = is.geometry().global(x);
      if (xglobal[0]<1e-6)
        {
          typename Traits::StateType u(0.0);
          u[0] = s[0];
          u[1] = 1.224*(1+0.5*sin(2*pi*1500.0*time));
          return u;
        }
      if (xglobal[0]>1.0-1e-6)
        {
          typename Traits::StateType u(0.0);
          return u;
        }
      typename Traits::StateType u(0.0);
      u[0] = s[0];
      u[2] = 0.0;
      return u;
    }

    //! right hand side
    typename Traits::StateType
    q (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
    {
      typename Traits::StateType rhs(0.0);
      return rhs;
    }

    //! initial value
    typename Traits::StateType
    u0 (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
    {
      typename Traits::StateType u(0.0);
      return u;
    }

    //! set time for subsequent evaluation
    void setTime (RF t)
    {
      time = t;
    }

  private:
    double pi;
    RF time;
  };



  }
}

#endif // DUNE_PDELAB_LOCALOPERATOR_LINEARACOUSTICSPARAMETER_HH
