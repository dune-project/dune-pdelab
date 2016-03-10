// -*- tab-width: 8; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_LOCALOPERATOR_LINEARELASTICITYPARAMETER_HH
#define DUNE_PDELAB_LOCALOPERATOR_LINEARELASTICITYPARAMETER_HH

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

    /** \brief Traits class for linear elasticity parameters
     *
     * A class supplying parameters to a linear-elasticity local
     * operator has to define a public traits class exporting the needed
     * types and constants.
     */
    template<typename GV, typename RF>
    struct LinearElasticityParameterTraits
    {
      //! \brief the grid view
      typedef GV GridViewType;

      //! \brief Enum for domain dimension
      enum {
        //! \brief dimension of the domain
        dimDomain = GV::dimension,
        dimRange = GV::dimension
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
      typedef Dune::FieldVector<RF,dimRange> RangeType;

      //! grid types
      typedef typename GV::Traits::template Codim<0>::Entity ElementType;
      typedef typename GV::Intersection IntersectionType;
    };

    //! base class for linear elasticity parameter class
    template<class T, class Imp>
    class LinearElasticityParameterInterface
      : public Dune::PDELab::DirichletConstraintsParameters
    {
    public:
      typedef T Traits;

      //! volume force term
      void
      f (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
         typename Traits::RangeType & y) const
      {
        asImp().f(e,x,y);
      }

      template<typename I>
      bool isDirichlet(const I & ig,
        const typename Traits::IntersectionDomainType & coord
        ) const
      {
        return asImp().isDirichlet( ig, coord );
      }

      //! Dirichlet boundary condition value (displacement)
      void
      u (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
        typename Traits::RangeType & y) const
      {
        asImp().g(e,x,y);
      }

      // //! Neumann boundary condition (surface force)
      // void
      // g (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
      //   typename Traits::RangeType & y) const
      // {
      //   asImp().g(e,x,y);
      // }

      //! First Lame parameter
      typename Traits::RangeFieldType
      lambda(const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
      {
        return asImp().lambda(e,x);
      }

      //! Second Lame parameter (shear modulus)
      typename Traits::RangeFieldType
      mu(const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
      {
        return asImp().mu(e,x);
      }

    private:
      Imp& asImp () {return static_cast<Imp &> (*this);}
      const Imp& asImp () const {return static_cast<const Imp &>(*this);}
    };

    /*! Adapter that extracts Dirichlet boundary conditions from parameter class

      \tparam T  model of LinearElasticityParameterInterface
    */
    template<typename T>
    class LinearElasticityDirichletExtensionAdapter
      : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<typename T::Traits::GridViewType,
                                                                               typename T::Traits::RangeFieldType,
                                                                               T::Traits::dimRange, typename T::Traits::RangeType >
                                              ,LinearElasticityDirichletExtensionAdapter<T> >
    {
    public:
      typedef Dune::PDELab::GridFunctionTraits<typename T::Traits::GridViewType,
                                               typename T::Traits::RangeFieldType,
                                               T::Traits::dimRange, typename T::Traits::RangeType > Traits;

      //! constructor
      LinearElasticityDirichletExtensionAdapter (const typename Traits::GridViewType& g_, T& t_)
        : g(g_), t(t_)
      {}

      //! \copydoc GridFunctionBase::evaluate()
      inline void evaluate (const typename Traits::ElementType& e,
        const typename Traits::DomainType& x,
        typename Traits::RangeType& y) const
      {
        t.u(e,x,y);
      }

      inline const typename Traits::GridViewType& getGridView () const
      {
        return g;
      }

    private:
      const typename Traits::GridViewType g;
      T& t;
    };

  }
}

#endif // DUNE_PDELAB_LOCALOPERATOR_LINEARELASTICITYPARAMETER_HH
