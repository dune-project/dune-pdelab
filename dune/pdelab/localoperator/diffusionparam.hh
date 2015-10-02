#ifndef DUNE_PDELAB_DIFFUSIONPARAM_HH
#define DUNE_PDELAB_DIFFUSIONPARAM_HH
#warning This file is deprecated and will be removed after the Dune-PDELab 2.4 release, use convectiondiffusionparameter.hh and the implementation therein instead!

#include <iostream>
#include <dune/common/deprecated.hh>

namespace Dune {
  namespace PDELab {
    //! \addtogroup LocalOperator
    //! \ingroup PDELab
    //! \{

    /** \brief Group types and methods to specify the boundary condition of a diffusion problem
     */
    struct DiffusionBoundaryCondition {
      //! Enum for the Boundary condition type of a diffusion problem
      enum Type {
        //! Neumann boundary condition (prescribed flux)
        Neumann = 0,
        //! Dirichlet boundary condition (prescribed value)
        Dirichlet = 1
      };

      //! Test for Dirichlet boundary condition
      static bool isDirichlet (Type i)
      {
        return (i == Dirichlet);
      }
      //! Test for Neumann boundary condition
      static bool isNeumann (Type i)
      {
        return (i == Neumann);
      }
    };

    /** \brief Adapter to get ConvectionDiffusion parameter object from the old style separate parameter grid functions
     *
     */
    template<typename K, typename A0, typename F, typename B, typename G, typename J>
    class DUNE_DEPRECATED_MSG("Deprecated in Dune-PDELab 2.4, use the parameter class ConvectionDiffusionModelProblem instead!") ConvectionDiffusion_Diffusion_Adapter
    {
      typedef typename F::Traits::RangeFieldType RF;
      typedef typename F::Traits::GridViewType GV;

      typedef ConvectionDiffusionBoundaryConditions::Type BCType;

    public:
      typedef ConvectionDiffusionParameterTraits<GV,RF> Traits;

      //! constructor
      DUNE_DEPRECATED_MSG("Deprecated in Dune-PDELab 2.4, use the parameter class ConvectionDiffusionModelProblem instead!")
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
    //! \}
  }
}

#endif // DUNE_PDELAB_DIFFUSIONPARAM_HH
