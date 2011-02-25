// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_CONVECTIONDIFFUSIONPARAMETER_HH
#define DUNE_PDELAB_CONVECTIONDIFFUSIONPARAMETER_HH

#include<vector>

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/static_assert.hh>
#include<dune/common/geometrytype.hh>
#include<dune/grid/common/genericreferenceelements.hh>
#include<dune/grid/common/quadraturerules.hh>
#include<dune/pdelab/common/geometrywrapper.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/functionutilities.hh>

#include"diffusionparam.hh"

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
      enum Type { Dirichlet=1, Neumann=-1, Outflow=-2 }; // BC requiring constraints must be >0 if
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
        typename Traits::DomainType xglobal = is.geometry().global(x);
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
      : public Dune::PDELab::BoundaryGridFunctionBase<Dune::PDELab::
                                                      BoundaryGridFunctionTraits<typename T::Traits::GridViewType,int,1,
                                                                                 Dune::FieldVector<int,1> >,
                                                      ConvectionDiffusionBoundaryConditionAdapter<T> >
    {
      const typename T::Traits::GridViewType& gv;
      const T& t;

    public:
      typedef Dune::PDELab::BoundaryGridFunctionTraits<typename T::Traits::GridViewType,int,1,
                                                       Dune::FieldVector<int,1> > Traits;
      typedef Dune::PDELab::BoundaryGridFunctionBase<Traits,ConvectionDiffusionBoundaryConditionAdapter<T> > BaseT;

      ConvectionDiffusionBoundaryConditionAdapter (const typename T::Traits::GridViewType& gv_, const T& t_) : gv(gv_), t(t_) {}

      template<typename I>
      inline void evaluate (const Dune::PDELab::IntersectionGeometry<I>& ig, 
                            const typename Traits::DomainType& x,
                            typename Traits::RangeType& y) const
      {  
        y = t.bctype(ig.intersection(),x);
      }

      //! get a reference to the GridView
      inline const typename T::Traits::GridViewType& getGridView ()
      {
        return gv;
      }
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
    ConvectionDiffusionDirichletExtensionAdapter (const typename Traits::GridViewType& g_, const T& t_) : g(g_), t(t_) {}

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
  
  private:
    const typename Traits::GridViewType& g;
    const T& t;
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
  const typename Traits::GridViewType& g;
  const T& t;
};

    /** \brief Adapter to get ConvectionDiffusion parameter object from the old style separate parameter grid functions
     * 
     */
    template<typename K, typename A0, typename F, typename B, typename G, typename J>
    class ConvectionDiffusion_Diffusion_Adapter
    {
      typedef typename F::Traits::RangeFieldType RF;
      typedef typename F::Traits::GridViewType GV;

      typedef ConvectionDiffusionBoundaryConditions::Type BCType;

    public:
      typedef ConvectionDiffusionParameterTraits<GV,RF> Traits;

      //! constructor
      ConvectionDiffusion_Diffusion_Adapter (const K& k_, const A0& a0_, const F& f_, const B& b_, const G& g_, const J& j_) :
        k__(k_), a0__(a0_), f__(f_), b__(b_), g__(g_), j__(j_)
      {}

      //! tensor diffusion coefficient
      typename Traits::PermTensorType
      A (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
      {
        typename K::Traits::RangeType tensor(0.0);
        k__.evaluate(e,x,tensor);
        return tensor;
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
        typename A0::Traits::RangeType y;
        a0__.evaluate(e,x,y);
        return y;
      }

      //! source term
      typename Traits::RangeFieldType 
      f (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
      {
        typename F::Traits::RangeType y;
        f__.evaluate(e,x,y);
        return y;
      }

      //! boundary condition type function
      BCType
      bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
      {
        typename B::Traits::RangeType bctype;
        b__.evaluate(Dune::PDELab::IntersectionGeometry<typename Traits::IntersectionType>(is,0),x,bctype);
        if (DiffusionBoundaryCondition::isDirichlet(bctype)) return ConvectionDiffusionBoundaryConditions::Dirichlet;
        return ConvectionDiffusionBoundaryConditions::Neumann;
      }

      //! Dirichlet boundary condition value
      typename Traits::RangeFieldType 
      g (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
      {
        typename G::Traits::RangeType y;
        g__.evaluate(e,x,y);
        return y;
      }

      //! Neumann boundary condition
      typename Traits::RangeFieldType 
      j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
      {
        typename J::Traits::RangeType y;
        j__.evaluate(*(is.inside()),is.geometryInInside().global(x),y);
        return y;
      }

      //! outflow boundary condition
      typename Traits::RangeFieldType 
      o (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
      {
        return 0.0;
      }

    private:
      const K& k__;
      const A0& a0__;
      const F& f__;
      const B& b__;
      const G& g__;
      const J& j__;
    };
  }
}


#endif
