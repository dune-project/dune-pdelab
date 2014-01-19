// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_CONVECTIONDIFFUSION_HH
#define DUNE_PDELAB_CONVECTIONDIFFUSION_HH

#include<vector>

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/static_assert.hh>
#include<dune/geometry/type.hh>

#include<dune/geometry/referenceelements.hh>
#include<dune/geometry/quadraturerules.hh>

#include <dune/pdelab/common/function.hh>
#include<dune/pdelab/constraints/common/constraintsparameters.hh>

#include"defaultimp.hh"
#include"pattern.hh"
#include"flags.hh"
#include"idefault.hh"


namespace Dune {
  namespace PDELab {
    //! \addtogroup LocalOperator
    //! \ingroup PDELab
    //! \{

    /** a local operator for solving the convection-diffusion equation
     *
     * \f{align*}{
     *   \nabla\cdot\{q(x,u) - D(x) v(u) \nabla w(u)\} &=& f(u) \mbox{ in } \Omega,  \\
     *                                            u &=& g \mbox{ on } \partial\Omega_D \\
     *         (q(x,u) - K(x)\nabla w(u)) \cdot \nu &=& j(u) \mbox{ on } \partial\Omega_N \\
     * \f}
     */

    //! traits class for two phase parameter class
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

    //! base class for parameter class
    template<class T, class Imp>
    class ConvectionDiffusionParameterInterface
    {
    public:
      typedef T Traits;

      //! source/reaction term
      typename Traits::RangeFieldType
      f (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
         typename Traits::RangeFieldType u) const
      {
        return asImp().f(e,x,u);
      }

      //! nonlinearity under gradient
      typename Traits::RangeFieldType
      w (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
         typename Traits::RangeFieldType u) const
      {
        return asImp().w(e,x,u);
      }

      //! scalar nonlinearity in diffusion coefficient
      typename Traits::RangeFieldType
      v (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
         typename Traits::RangeFieldType u) const
      {
        return asImp().v(e,x,u);
      }

      //! tensor diffusion coefficient
      typename Traits::PermTensorType
      D (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
      {
        return asImp().D(e,x);
      }

      //! nonlinear flux vector
      typename Traits::RangeType
      q (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
         typename Traits::RangeFieldType u) const
      {
        return asImp().q(e,x,u);
      }

      template<typename I>
      bool isDirichlet(
                       const I & intersection,               /*@\label{bcp:name}@*/
                       const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord
                       ) const
      {
        return asImp().isDirichlet( intersection, coord );
      }

      //! Dirichlet boundary condition value
      typename Traits::RangeFieldType
      g (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
      {
        return asImp().g(e,x);
      }

      //! Neumann boundary condition
      // Good: The dependence on u allows us to implement Robin type boundary conditions.
      // Bad: This interface cannot be used for mixed finite elements where the flux is the essential b.c.
      typename Traits::RangeFieldType
      j (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
         typename Traits::RangeFieldType u) const
      {
        return asImp().j(e,x);
      }

    private:
      Imp& asImp () {return static_cast<Imp &> (*this);}
      const Imp& asImp () const {return static_cast<const Imp &>(*this);}
    };


    /*! Adapter that extracts boundary condition type function from parameter class

      \tparam T  model of ConvectionDiffusionParameterInterface
    */
    template<typename T>
    class BCTypeParam_CD
      : public Dune::PDELab::DirichletConstraintsParameters   /*@\label{bcp:base}@*/
    {
      const typename T::Traits::GridViewType gv;
      const T& t;

    public:
      BCTypeParam_CD( const typename T::Traits::GridViewType& gv_, const T& t_ )
        : gv( gv_ ), t( t_ )
      {
      }

      template<typename I>
      bool isDirichlet(
                       const I & intersection,               /*@\label{bcp:name}@*/
                       const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord
                       ) const
      {
        return t.isDirichlet( intersection, coord );
      }
    };


    /*! Adapter that extracts Dirichlet boundary conditions from parameter class

      \tparam T  model of ConvectionDiffusionParameterInterface
    */
    template<typename T>
    class DirichletBoundaryCondition_CD
      : public GridFunctionBase<Dune::PDELab::GridFunctionTraits<typename T::Traits::GridViewType,
                                                                 typename T::Traits::RangeFieldType,
                                                                 1,Dune::FieldVector<typename T::Traits::RangeFieldType,1> >
                                ,DirichletBoundaryCondition_CD<T> >
    {
    public:
      typedef Dune::PDELab::GridFunctionTraits<typename T::Traits::GridViewType,
                                               typename T::Traits::RangeFieldType,
                                               1,Dune::FieldVector<typename T::Traits::RangeFieldType,1> > Traits;

      //! constructor
      DirichletBoundaryCondition_CD (const typename Traits::GridViewType& g_, const T& t_) : g(g_), t(t_) {}

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
      typename Traits::GridViewType g;
      const T& t;
    };


    /** a local operator for solving the convection-diffusion equation defined above
     *
     * with conforming finite elements on all types of grids in any dimension
     * \tparam T model of ConvectionDiffusionParameterInterface
     */
    template<typename T>
    class ConvectionDiffusion :
      public NumericalJacobianApplyVolume<ConvectionDiffusion<T> >,
      public NumericalJacobianApplyBoundary<ConvectionDiffusion<T> >,
      public NumericalJacobianVolume<ConvectionDiffusion<T> >,
      public NumericalJacobianBoundary<ConvectionDiffusion<T> >,
      public FullVolumePattern,
      public LocalOperatorDefaultFlags,
      public InstationaryLocalOperatorDefaultMethods<double>
    {
    public:
      // pattern assembly flags
      enum { doPatternVolume = true };

      // residual assembly flags
      enum { doAlphaVolume = true };
      enum { doAlphaBoundary = true };

      ConvectionDiffusion (T& param_, int intorder_=2)
        : param(param_), intorder(intorder_)
      {}

      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
        // domain and range field type
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::JacobianType JacobianType;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RangeType;

        typedef typename LFSU::Traits::SizeType size_type;

        // dimensions
        const int dim = EG::Geometry::dimension;

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

        // evaluate diffusion tensor at cell center, assume it is constant over elements
        typename T::Traits::PermTensorType tensor;
        Dune::FieldVector<DF,dim> localcenter = Dune::ReferenceElements<DF,dim>::general(gt).position(0,0);
        tensor = param.D(eg.entity(),localcenter);

        // evaluate nonlinearity w(x_i); we assume here it is a Lagrange basis!
        std::vector<typename T::Traits::RangeFieldType> w(lfsu.size());
        for (size_type i=0; i<lfsu.size(); i++)
          w[i] = param.w(eg.entity(),localcenter,x(lfsu,i));

        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // evaluate basis functions
            std::vector<RangeType> phi(lfsu.size());
            lfsu.finiteElement().localBasis().evaluateFunction(it->position(),phi);

            // evaluate u
            RF u=0.0;
            for (size_type i=0; i<lfsu.size(); i++)
              u += w[i]*phi[i];

            // evaluate source term
            typename T::Traits::RangeFieldType f = param.f(eg.entity(),it->position(),u);

            // evaluate flux term
            typename T::Traits::RangeType q = param.q(eg.entity(),it->position(),u);

            // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
            std::vector<JacobianType> js(lfsu.size());
            lfsu.finiteElement().localBasis().evaluateJacobian(it->position(),js);

            // transform gradients of shape functions to real element
            const typename EG::Geometry::JacobianInverseTransposed jac
              = eg.geometry().jacobianInverseTransposed(it->position());
            std::vector<Dune::FieldVector<RF,dim> > gradphi(lfsu.size());
            for (size_type i=0; i<lfsu.size(); i++)
              {
                gradphi[i] = 0.0;
                jac.umv(js[i][0],gradphi[i]);
              }

            // v(u) compute gradient of u
            Dune::FieldVector<RF,dim> vgradu(0.0);
            for (size_type i=0; i<lfsu.size(); i++)
              vgradu.axpy(w[i],gradphi[i]);
            vgradu *= param.v(eg.entity(),it->position(),u);

            // compute D * v(u) * gradient of u
            Dune::FieldVector<RF,dim> Dvgradu(0.0);
            tensor.umv(vgradu,Dvgradu);

            // integrate (K grad u)*grad phi_i + a_0*u*phi_i
            RF factor = it->weight() * eg.geometry().integrationElement(it->position());
            for (size_type i=0; i<lfsu.size(); i++)
              r.accumulate(lfsu,i,( Dvgradu*gradphi[i] - q*gradphi[i] - f*phi[i] )*factor);
          }
      }

      // boundary integral
      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_boundary (const IG& ig,
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           R& r_s) const
      {
        // domain and range field type
        typedef typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RangeType;

        typedef typename LFSV::Traits::SizeType size_type;

        // dimensions
        const int dim = IG::dimension;

        // select quadrature rule
        Dune::GeometryType gtface = ig.geometryInInside().type();
        const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,intorder);

        // evaluate nonlinearity w(x_i); we assume here it is a Lagrange basis!
        Dune::FieldVector<DF,dim-1> facecenterlocal = Dune::ReferenceElements<DF,dim-1>::general(gtface).position(0,0);
        Dune::FieldVector<DF,dim> facecenterinelement = ig.geometryInInside().global( facecenterlocal );
        std::vector<typename T::Traits::RangeFieldType> w(lfsu_s.size());
        for (size_type i=0; i<lfsu_s.size(); i++)
          w[i] = param.w(*(ig.inside()),facecenterinelement,x_s(lfsu_s,i));

        // loop over quadrature points and integrate normal flux
        for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // evaluate boundary condition type
            // skip rest if we are on Dirichlet boundary
            if( param.isDirichlet( ig.intersection(), it->position() ) )
              continue;

            // position of quadrature point in local coordinates of element
            Dune::FieldVector<DF,dim> local = ig.geometryInInside().global(it->position());

            // evaluate test shape functions
            std::vector<RangeType> phi(lfsv_s.size());
            lfsv_s.finiteElement().localBasis().evaluateFunction(local,phi);

            // evaluate u
            RF u=0.0;
            for (size_type i=0; i<lfsu_s.size(); i++)
              u += w[i]*phi[i];

            // evaluate flux boundary condition
            typename T::Traits::RangeFieldType j;
            j = param.j(*(ig.inside()),local,u);

            // integrate j
            RF factor = it->weight()*ig.geometry().integrationElement(it->position());
            for (size_type i=0; i<lfsv_s.size(); i++)
              r_s.accumulate(lfsu_s,i,j*phi[i]*factor);
          }
      }

      //! set time in parameter class
      void setTime (double t)
      {
        param.setTime(t);
      }

    private:
      T& param;
      int intorder;
    };

    //! \} group LocalOperator
  } // namespace PDELab
} // namespace Dune

#endif
