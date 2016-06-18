// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_LOCALOPERATOR_CONVECTIONDIFFUSIONPARAMETER_HH
#define DUNE_PDELAB_LOCALOPERATOR_CONVECTIONDIFFUSIONPARAMETER_HH

#include<vector>

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/geometry/type.hh>
#include<dune/geometry/referenceelements.hh>
#include<dune/geometry/quadraturerules.hh>
#include<dune/pdelab/common/geometrywrapper.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/functionutilities.hh>
#include<dune/pdelab/constraints/common/constraintsparameters.hh>

namespace Dune {
  namespace PDELab {

    /** \brief Traits class for convection diffusion parameters
     *
     * A class supplying parameters to a convection-diffusion local
     * operator has to define a public traits class exporting the needed
     * types and constants.
     */
    template<typename GV, typename RF>
    struct ConvectionDiffusionParameterTraits
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
      typedef Dune::FieldVector<RF,GV::dimensionworld> RangeType;

      //! \brief permeability tensor type
      typedef Dune::FieldMatrix<RangeFieldType,dimDomain,dimDomain> PermTensorType;

      //! grid types
      typedef typename GV::Traits::template Codim<0>::Entity ElementType;
      typedef typename GV::Intersection IntersectionType;
    };

    /** \brief Class to define the boundary condition types
     */
    struct ConvectionDiffusionBoundaryConditions
    {
      enum Type { Dirichlet=1, Neumann=-1, Outflow=-2, None=-3 }; // BC requiring constraints must be >0 if
      // constraints assembler coming with PDELab is used
    };

    /** \brief Parameter class for solving the linear convection-diffusion equation
     *
     * A parameter class for the linear convection-diffusion equation
     * \f{align*}{
     *   -\nabla\cdot(A(x) \nabla u) + b(x)\cdot \nabla u + c(x)u &=& f \mbox{ in } \Omega,  \\
     *                                                          u &=& g \mbox{ on } \partial\Omega_D \\
     *                            (b(x,u) - A(x)\nabla u) \cdot n &=& j \mbox{ on } \partial\Omega_N \\
     *                                    -(A(x)\nabla u) \cdot n &=& o \mbox{ on } \partial\Omega_O
     * \f}
     * Note:
     *  - This formulation is valid for velocity fields which are non-divergence free.
     *  - Outflow boundary conditions should only be set on the outflow boundary
     *
     * \tparam T a traits class defining the necessary types
     */
    template<typename GV, typename RF>
    class ConvectionDiffusionModelProblem
    {
      typedef ConvectionDiffusionBoundaryConditions::Type BCType;

    public:
      typedef ConvectionDiffusionParameterTraits<GV,RF> Traits;

      //! tensor diffusion coefficient
      typename Traits::PermTensorType
      A (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
      {
        typename Traits::PermTensorType I;
        for (std::size_t i=0; i<Traits::dimDomain; i++)
          for (std::size_t j=0; j<Traits::dimDomain; j++)
            I[i][j] = (i==j) ? 1 : 0;
        return I;
      }

      //! velocity field
      typename Traits::RangeType
      b (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
      {
        typename Traits::RangeType v(0.0);
        return v;
      }

      //! sink term
      typename Traits::RangeFieldType
      c (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
      {
        return 0.0;
      }

      //! source term
      typename Traits::RangeFieldType
      f (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
      {
        return 0.0;
      }

      //! boundary condition type function
      BCType
      bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
      {
        return ConvectionDiffusionBoundaryConditions::Dirichlet;
      }

      //! Dirichlet boundary condition value
      typename Traits::RangeFieldType
      g (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
      {
        typename Traits::DomainType xglobal = e.geometry().global(x);
        return xglobal.two_norm();
      }

      //! Neumann boundary condition
      typename Traits::RangeFieldType
      j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
      {
        return 0.0;
      }

      //! outflow boundary condition
      typename Traits::RangeFieldType
      o (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
      {
        return 0.0;
      }
    };


    /*! Adapter that extracts boundary condition type function from parameter class

      \tparam T  model of ConvectionDiffusionParameterInterface
    */
    template<typename T>
    class ConvectionDiffusionBoundaryConditionAdapter
      :
      public Dune::PDELab::FluxConstraintsParameters,
      public Dune::PDELab::DirichletConstraintsParameters   /*@\label{bcp:base}@*/
    {
      const T& t;

    public:

      ConvectionDiffusionBoundaryConditionAdapter(const typename T::Traits::GridViewType& gv_,
                                                  const T& t_ )
        : t( t_ )
      {}

      ConvectionDiffusionBoundaryConditionAdapter(const T& t_ )
        : t( t_ )
      {}

      template<typename I>
      bool isDirichlet(const I & ig               /*@\label{bcp:name}@*/
                       , const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord
                       ) const
      {
        return( t.bctype( ig.intersection(), coord )
                == ConvectionDiffusionBoundaryConditions::Dirichlet );
      }

      template<typename I>
      bool isNeumann(const I & ig,   /*@\label{bcp:name}@*/
                     const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord
                     ) const
      {
        return !isDirichlet( ig, coord );
      }

    };





    /*! Adapter that extracts the flux boundary conditions from the parameter class

      \tparam T  model of ConvectionDiffusionParameterInterface
    */
    template<typename T>
    class ConvectionDiffusionVelocityExtensionAdapter
      : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits
                                                      <typename T::Traits::GridViewType,
                                                       typename T::Traits::RangeFieldType,
                                                       T::Traits::GridViewType::dimension>,
                                                      ConvectionDiffusionVelocityExtensionAdapter<T> >
    {
    public:
      typedef Dune::PDELab::AnalyticGridFunctionTraits<typename T::Traits::GridViewType,
                                                       typename T::Traits::RangeFieldType,
                                                       T::Traits::GridViewType::dimension> Traits;
      typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,ConvectionDiffusionVelocityExtensionAdapter<T> > BaseT;


      //! constructor
      ConvectionDiffusionVelocityExtensionAdapter (const typename Traits::GridViewType& gv_, T& t_)
        : BaseT(gv_), gv(gv_), t(t_)
      {}

      inline void evaluateGlobal (const typename Traits::DomainType& x,
                                  typename Traits::RangeType& y) const
      {
        y = t.b(x);
      }

      //! \copydoc GridFunctionBase::evaluate()
      inline void evaluate (const typename Traits::ElementType& e,
                            const typename Traits::DomainType& x,
                            typename Traits::RangeType& y) const
      {
        y = t.b(e,x);
      }

      inline const typename Traits::GridViewType& getGridView () const
      {
        return gv;
      }

      inline void setTime(double time_)
      {
        t.setTime(time_);
      }

    private:
      const typename Traits::GridViewType gv;
      T& t;
    };





  /*! Adapter that extracts Dirichlet boundary conditions from parameter class

    \tparam T  model of ConvectionDiffusionParameterInterface
  */
  template<typename T>
  class ConvectionDiffusionDirichletExtensionAdapter
    : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<typename T::Traits::GridViewType,
                                                                             typename T::Traits::RangeFieldType,
                                                                             1,Dune::FieldVector<typename T::Traits::RangeFieldType,1> >
                                            ,ConvectionDiffusionDirichletExtensionAdapter<T> >
  {
  public:
    typedef Dune::PDELab::GridFunctionTraits<typename T::Traits::GridViewType,
                                             typename T::Traits::RangeFieldType,
                                             1,Dune::FieldVector<typename T::Traits::RangeFieldType,1> > Traits;

    //! constructor
    ConvectionDiffusionDirichletExtensionAdapter (const typename Traits::GridViewType& g_, T& t_)
    : g(g_), t(t_)
    {}

    //! \copydoc GridFunctionBase::evaluate()
    inline void evaluate (const typename Traits::ElementType& e,
                          const typename Traits::DomainType& x,
                          typename Traits::RangeType& y) const
    {
      y = t.g(e,x);
    }

    inline const typename Traits::GridViewType& getGridView () const
    {
      return g;
    }

    inline void setTime(double time_)
    {
      t.setTime(time_);
    }

  private:
    const typename Traits::GridViewType g;
    T& t;
  };


/*! Adapter that extracts gradient of exact solution from parameter class

  \tparam T  model of ConvectionDiffusionParameterInterface
*/
template<typename T>
class ConvectionDiffusionExactGradientAdapter
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<typename T::Traits::GridViewType,
                                                                           typename T::Traits::RangeFieldType,
                                                                           T::Traits::GridViewType::dimension,Dune::FieldVector<typename T::Traits::RangeFieldType,T::Traits::GridViewType::dimension> >
                                          ,ConvectionDiffusionExactGradientAdapter<T> >
{
public:
  typedef Dune::PDELab::GridFunctionTraits<typename T::Traits::GridViewType,
                                           typename T::Traits::RangeFieldType,
                                           T::Traits::GridViewType::dimension,Dune::FieldVector<typename T::Traits::RangeFieldType,T::Traits::GridViewType::dimension> > Traits;

  //! constructor
  ConvectionDiffusionExactGradientAdapter (const typename Traits::GridViewType& g_, const T& t_) : g(g_), t(t_) {}

  //! \copydoc GridFunctionBase::evaluate()
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    y = t.gradient(e,x);
  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return g;
  }

private:
  const typename Traits::GridViewType g;
  const T& t;
};
  }
}


#endif // DUNE_PDELAB_LOCALOPERATOR_CONVECTIONDIFFUSIONPARAMETER_HH
