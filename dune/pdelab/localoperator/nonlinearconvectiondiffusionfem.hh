// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_LOCALOPERATOR_NONLINEARCONVECTIONDIFFUSIONFEM_HH
#define DUNE_PDELAB_LOCALOPERATOR_NONLINEARCONVECTIONDIFFUSIONFEM_HH

#include<vector>

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/geometry/type.hh>

#include<dune/geometry/referenceelements.hh>

#include<dune/pdelab/common/quadraturerules.hh>
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

    /** a local operator for solving the non-linear convection-diffusion equation with standard FEM
     *
     * \f{align*}{
     *   \nabla\cdot\{q(x,u) - D(x) v(u) \nabla w(u)\} &=& f(u) \mbox{ in } \Omega,  \\
     *                                            u &=& g \mbox{ on } \partial\Omega_D \\
     *         (q(x,u) - K(x)\nabla w(u)) \cdot \nu &=& j(u) \mbox{ on } \partial\Omega_N \\
     * \f}
     */

    //! traits class for two phase parameter class
    template<typename GV, typename RF>
    struct NonLinearConvectionDiffusionParameterTraits
    {
      //! \brief the grid view
      using GridViewType = GV;

      //! \brief Enum for domain dimension
      enum {
        //! \brief dimension of the domain
        dimDomain = GV::dimension
      };

      //! \brief Export type for domain field
      using DomainFieldType = typename GV::Grid::ctype;

      //! \brief domain type
      using DomainType = Dune::FieldVector<DomainFieldType,dimDomain>;

      //! \brief domain type
      using IntersectionDomainType = Dune::FieldVector<DomainFieldType,dimDomain-1>;

      //! \brief Export type for range field
      using RangeFieldType = RF;

      //! \brief range type
      using RangeType = Dune::FieldVector<RF,GV::dimensionworld>;

      //! \brief permeability tensor type
      using PermTensorType = Dune::FieldMatrix<RangeFieldType,dimDomain,dimDomain>;

      //! grid types
      using ElementType = typename GV::Traits::template Codim<0>::Entity;
      using IntersectionType = typename GV::Intersection;
    };

    //! base class for parameter class
    template<class T, class Imp>
    class NonLinearConvectionDiffusionParameterInterface
    {
    public:
      using Traits = T;

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
                       const Dune::FieldVector<typename I::ctype, I::mydimension> & coord
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

      \tparam T  model of NonLinearConvectionDiffusionParameterInterface
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
                       const Dune::FieldVector<typename I::ctype, I::mydimension> & coord
                       ) const
      {
        return t.isDirichlet( intersection, coord );
      }
    };


    /*! Adapter that extracts Dirichlet boundary conditions from parameter class

      \tparam T  model of NonLinearConvectionDiffusionParameterInterface
    */
    template<typename T>
    class DirichletBoundaryCondition_CD
      : public GridFunctionBase<Dune::PDELab::GridFunctionTraits<typename T::Traits::GridViewType,
                                                                 typename T::Traits::RangeFieldType,
                                                                 1,Dune::FieldVector<typename T::Traits::RangeFieldType,1> >
                                ,DirichletBoundaryCondition_CD<T> >
    {
    public:
      using Traits = Dune::PDELab::GridFunctionTraits<typename T::Traits::GridViewType,
                                                      typename T::Traits::RangeFieldType,
                                                      1,Dune::FieldVector<typename T::Traits::RangeFieldType,1> >;

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
     * \tparam T model of NonLinearConvectionDiffusionParameterInterface
     */
    template<typename T>
    class NonLinearConvectionDiffusionFEM :
      public NumericalJacobianApplyVolume<NonLinearConvectionDiffusionFEM<T> >,
      public NumericalJacobianApplyBoundary<NonLinearConvectionDiffusionFEM<T> >,
      public NumericalJacobianVolume<NonLinearConvectionDiffusionFEM<T> >,
      public NumericalJacobianBoundary<NonLinearConvectionDiffusionFEM<T> >,
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

      NonLinearConvectionDiffusionFEM (T& param_, int intorder_=2)
        : param(param_), intorder(intorder_)
      {}

      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
        // define types
        using RF = typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType;
        using JacobianType = typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::JacobianType;
        using RangeType = typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType;
        using size_type = typename LFSU::Traits::SizeType;

        // dimensions
        const int dim = EG::Geometry::mydimension;

        // Reference to cell
        const auto& cell = eg.entity();

        // select quadrature rule
        auto geo = eg.geometry();

        // evaluate diffusion tensor at cell center, assume it is constant over elements
        auto ref_el = referenceElement(geo);
        auto localcenter = ref_el.position(0,0);
        auto tensor = param.D(cell,localcenter);

        // evaluate nonlinearity w(x_i); we assume here it is a Lagrange basis!
        std::vector<typename T::Traits::RangeFieldType> w(lfsu.size());
        for (size_type i=0; i<lfsu.size(); i++)
          w[i] = param.w(cell,localcenter,x(lfsu,i));

        // Transformation
        typename EG::Geometry::JacobianInverseTransposed jac;

        // Initialize vectors outside for loop
        std::vector<RangeType> phi(lfsu.size());
        std::vector<JacobianType> js(lfsu.size());
        std::vector<Dune::FieldVector<RF,dim> > gradphi(lfsu.size());
        Dune::FieldVector<RF,dim> vgradu(0.0);
        Dune::FieldVector<RF,dim> Dvgradu(0.0);

        // loop over quadrature points
        for (const auto& ip : quadratureRule(geo,intorder))
          {
            // evaluate basis functions
            lfsu.finiteElement().localBasis().evaluateFunction(ip.position(),phi);

            // evaluate u
            RF u=0.0;
            for (size_type i=0; i<lfsu.size(); i++)
              u += w[i]*phi[i];

            // evaluate source term
            auto f = param.f(cell,ip.position(),u);

            // evaluate flux term
            auto q = param.q(cell,ip.position(),u);

            // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
            lfsu.finiteElement().localBasis().evaluateJacobian(ip.position(),js);

            // transform gradients of shape functions to real element
            jac = geo.jacobianInverseTransposed(ip.position());
            for (size_type i=0; i<lfsu.size(); i++)
              {
                gradphi[i] = 0.0;
                jac.umv(js[i][0],gradphi[i]);
              }

            // v(u) compute gradient of u
            vgradu = 0.0;
            for (size_type i=0; i<lfsu.size(); i++)
              vgradu.axpy(w[i],gradphi[i]);
            vgradu *= param.v(cell,ip.position(),u);

            // compute D * v(u) * gradient of u
            Dvgradu = 0.0;
            tensor.umv(vgradu,Dvgradu);

            // integrate (K grad u)*grad phi_i + a_0*u*phi_i
            auto factor = ip.weight() * geo.integrationElement(ip.position());
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
        // define types
        using RF = typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType;
        using RangeType = typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType;
        using size_type = typename LFSV::Traits::SizeType;

        // get inside cell entity
        const auto& cell_inside = ig.inside();

        // get geometry
        auto geo = ig.geometry();

        // Get geometry of intersection in local coordinates of cell_inside
        auto geo_in_inside = ig.geometryInInside();

        // evaluate nonlinearity w(x_i); we assume here it is a Lagrange basis!
        auto ref_el_in_inside = referenceElement(geo_in_inside);
        auto local_face_center = ref_el_in_inside.position(0,0);
        auto face_center_in_element = geo_in_inside.global(local_face_center);
        std::vector<typename T::Traits::RangeFieldType> w(lfsu_s.size());
        for (size_type i=0; i<lfsu_s.size(); i++)
          w[i] = param.w(cell_inside,face_center_in_element,x_s(lfsu_s,i));

        // Initialize vectors outside for loop
        std::vector<RangeType> phi(lfsv_s.size());

        // loop over quadrature points and integrate normal flux
        for (const auto& ip : quadratureRule(geo,intorder))
          {
            // evaluate boundary condition type
            // skip rest if we are on Dirichlet boundary
            if( param.isDirichlet( ig.intersection(), ip.position() ) )
              continue;

            // position of quadrature point in local coordinates of element
            auto local = geo_in_inside.global(ip.position());

            // evaluate test shape functions
            lfsv_s.finiteElement().localBasis().evaluateFunction(local,phi);

            // evaluate u
            RF u=0.0;
            for (size_type i=0; i<lfsu_s.size(); i++)
              u += w[i]*phi[i];

            // evaluate flux boundary condition
            auto j = param.j(cell_inside,local,u);

            // integrate j
            auto factor = ip.weight()*geo.integrationElement(ip.position());
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

#endif // DUNE_PDELAB_LOCALOPERATOR_NONLINEARCONVECTIONDIFFUSIONFEM_HH
