// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_STOKESDG_VECFEM_HH
#define DUNE_PDELAB_STOKESDG_VECFEM_HH

#warning This file is deprecated and will be removed after the DUNE-PDELab 2.4 release! Include the header dune/pdelab/localoperator/dgnavierstokesvelvecfem.hh instead!

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parametertree.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/localfunctions/common/interfaceswitch.hh>
#include <dune/pdelab/localoperator/idefault.hh>

#include "defaultimp.hh"
#include "pattern.hh"
#include "flags.hh"
#include "stokesdg.hh"

#ifndef VBLOCK
#define VBLOCK 0
#endif
#define PBLOCK (- VBLOCK + 1)

namespace Dune {
  namespace PDELab {

  template<class Basis, class Dummy = void>
  struct VectorBasisInterfaceSwitch {
    //! export vector type of the local coordinates
    typedef typename Basis::Traits::DomainLocal DomainLocal;
    //! export field type of the values
    typedef typename Basis::Traits::RangeField RangeField;
    //! export dimension of the values
    static const std::size_t dimRange = Basis::Traits::dimRange;

    //! Compute global jacobian matrix for vector valued bases
    template<typename Geometry>
    static void jacobian(const Basis& basis, const Geometry& geometry,
                         const DomainLocal& xl,
                         std::vector<FieldMatrix<RangeField, dimRange,
                                          Geometry::coorddimension> >& jac)
    {
      jac.resize(basis.size());
      basis.evaluateJacobian(xl, jac);
    }
  };

  //! Switch for uniform treatment of local and global basis classes
  template<class Basis>
  struct VectorBasisInterfaceSwitch<
    Basis, typename std::enable_if<Basis::Traits::dimDomain>::type
    >
  {
    //! export vector type of the local coordinates
    typedef typename Basis::Traits::DomainType DomainLocal;
    //! export field type of the values
    typedef typename Basis::Traits::RangeFieldType RangeField;
    //! export dimension of the values
    static const std::size_t dimRange = Basis::Traits::dimRange;

    //! Compute global jacobian matrix for vector valued bases
    template<typename Geometry>
    static void jacobian(const Basis& basis, const Geometry& geometry,
                         const DomainLocal& xl,
                         std::vector<FieldMatrix<RangeField, dimRange,
                                          Geometry::coorddimension> >& jac)
    {
      jac.resize(basis.size());

      std::vector<FieldMatrix<
      RangeField, dimRange, Geometry::coorddimension> > ljac(basis.size());
      basis.evaluateJacobian(xl, ljac);

      const typename Geometry::Jacobian& geojac =
        geometry.jacobianInverseTransposed(xl);

      for(std::size_t i = 0; i < basis.size(); ++i)
        for(std::size_t row=0; row < dimRange; ++row)
          geojac.mv(ljac[i][row], jac[i][row]);
    }
  };

    /** \brief A local operator for solving the stokes equation using
            a DG discretization and a vector valued finite element map
            for the velocity grid function space.

        \tparam F velocity source term function
        \tparam B boundary condition function
        \tparam V dirichlet velocity boundary condition function
        \tparam P dirichlet pressure boundary condition function
        \tparam IP a class providing the interior penalty factor for each face
    */
    template<typename F, typename B, typename V, typename P,
             typename IP = DefaultInteriorPenalty<typename V::Traits::RangeFieldType> >
    class DUNE_DEPRECATED_MSG("Deprecated in DUNE-PDELab 2.4, use DGNavierStokesVelVecFEM instead!") StokesDGVectorFEM :
      public LocalOperatorDefaultFlags,
      public FullSkeletonPattern, public FullVolumePattern,
      public Dune::PDELab::NumericalJacobianApplyVolume<StokesDGVectorFEM<F,B,V,P,IP> >,
      public Dune::PDELab::NumericalJacobianVolume<StokesDGVectorFEM<F,B,V,P,IP> >,
      public Dune::PDELab::NumericalJacobianApplySkeleton<StokesDGVectorFEM<F,B,V,P,IP> >,
      public Dune::PDELab::NumericalJacobianSkeleton<StokesDGVectorFEM<F,B,V,P,IP> >,
      public Dune::PDELab::NumericalJacobianApplyBoundary<StokesDGVectorFEM<F,B,V,P,IP> >,
      public Dune::PDELab::NumericalJacobianBoundary<StokesDGVectorFEM<F,B,V,P,IP> >,
      public InstationaryLocalOperatorDefaultMethods<double>
    {
      typedef StokesBoundaryCondition BC;
      typedef typename V::Traits::RangeFieldType RF;

    public:
      typedef IP InteriorPenaltyFactor;

      // pattern assembly flags
      enum { doPatternVolume = true };
      enum { doPatternSkeleton = true };

      // residual assembly flags
      enum { doAlphaVolume    = true };
      enum { doAlphaSkeleton  = true };
      enum { doAlphaBoundary  = true };

      DUNE_DEPRECATED_MSG("Deprecated in DUNE-PDELab 2.4, use DGNavierStokesVelVecFEM instead!")
      StokesDGVectorFEM (const Dune::ParameterTree & configuration,const IP & ip_factor_, const RF mu_,
                         const F & _f, const B & _b, const V & _v, const P & _p, int _qorder=4) :
        f(_f), b(_b), v(_v), p(_p), qorder(_qorder), mu(mu_), ip_factor(ip_factor_)
      {
        epsilon = configuration.get<int>("epsilon");
      }

      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
        // dimensions
        static const unsigned int dim = EG::Geometry::dimension;
        static const unsigned int dimw = EG::Geometry::dimensionworld;

        // subspaces
        static_assert
          ((LFSV::CHILDREN == 2), "You seem to use the wrong function space for StokesDGVectorFEM");

        typedef typename LFSV::template Child<VBLOCK>::Type LFSV_V;
        const LFSV_V& lfsv_v = lfsv.template child<VBLOCK>();
        const LFSV_V& lfsu_v = lfsu.template child<VBLOCK>();

        typedef typename LFSV::template Child<PBLOCK>::Type LFSV_P;
        const LFSV_P& lfsv_p = lfsv.template child<PBLOCK>();
        const LFSV_P& lfsu_p = lfsu.template child<PBLOCK>();

        // domain and range field type
        typedef FiniteElementInterfaceSwitch<typename LFSV_V::Traits::FiniteElementType > FESwitch_V;
        typedef BasisInterfaceSwitch<typename FESwitch_V::Basis > BasisSwitch_V;
        typedef VectorBasisInterfaceSwitch<typename FESwitch_V::Basis > VectorBasisSwitch_V;
        typedef FiniteElementInterfaceSwitch<typename LFSV_P::Traits::FiniteElementType > FESwitch_P;
        typedef BasisInterfaceSwitch<typename FESwitch_P::Basis > BasisSwitch_P;
        typedef typename BasisSwitch_V::DomainField DF;
        typedef typename BasisSwitch_V::RangeField RF;
        typedef typename BasisSwitch_V::Range Range_V;
        typedef typename BasisSwitch_P::Range Range_P;
        typedef typename LFSV::Traits::SizeType size_type;

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,qorder);

        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // Compute weight and coordinates
            const Dune::FieldVector<DF,dim> local = it->position();
            const Dune::FieldVector<DF,dimw> global = eg.geometry().global(local);
            const RF weight = it->weight() * eg.geometry().integrationElement(it->position());

            // values of velocity shape functions
            std::vector<Range_V> phi_v(lfsv_v.size());
            FESwitch_V::basis(lfsv_v.finiteElement()).evaluateFunction(local,phi_v);

            // values of velocity shape functions
            std::vector<Range_P> phi_p(lfsv_p.size());
            FESwitch_P::basis(lfsv_p.finiteElement()).evaluateFunction(local,phi_p);

            // evaluate jacobian of basis functions on reference element
            std::vector<Dune::FieldMatrix<RF,dim,dim> > jac_phi_v(lfsu_v.size());
            VectorBasisSwitch_V::jacobian
              (FESwitch_V::basis(lfsv_v.finiteElement()), eg.geometry(), it->position(), jac_phi_v);

            // compute divergence of test functions
            std::vector<RF> div_phi_v(lfsv_v.size(),0.0);
            for (size_type i=0; i<lfsv_v.size(); i++)
              for (size_type d=0; d<dim; d++)
                div_phi_v[i] += jac_phi_v[i][d][d];

            // compute velocity value, jacobian, and divergence
            Range_V val_u(0.0);
            Dune::FieldMatrix<RF,dim,dim> jac_u(0.0);
            RF div_u(0.0);
            for (size_type i=0; i<lfsu_v.size(); i++){
              val_u.axpy(x(lfsu_v,i),phi_v[i]);
              jac_u.axpy(x(lfsu_v,i),jac_phi_v[i]);
              div_u += x(lfsu_v,i) * div_phi_v[i];
            }

            // compute pressure value
            Range_P val_p(0.0);
            for (size_type i=0; i<lfsu_p.size(); i++)
              val_p.axpy(x(lfsu_p,i),phi_p[i]);

            // evaluate source term
            typename F::Traits::RangeType fval;
            f.evaluateGlobal(global,fval);

            {// Integrate (f*v)
              const RF factor = weight;
              for (size_type i=0; i<lfsv_v.size(); i++)
                r.accumulate(lfsv_v, i, (fval * phi_v[i]) * factor );
            }

            {// Integrate (mu * d_i u_j d_i v_j)
              const RF factor = mu * weight;
              for (size_type i=0; i<lfsv_v.size(); i++){
                RF dvdu(0); contraction(jac_u,jac_phi_v[i],dvdu);
                r.accumulate(lfsv_v, i, dvdu * factor);
              }
            }

            {// Integrate - p div v
              for (size_type i=0; i<lfsv_v.size(); i++)
                r.accumulate(lfsv_v, i, - div_phi_v[i] * val_p * weight);
            }

            {// Integrate - q div u
              for (size_type i=0; i<lfsv_p.size(); i++)
                r.accumulate(lfsv_p, i, - div_u * phi_p[i] * weight);
            }

          }
      }

      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_skeleton (const IG& ig,
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                           R& r_s, R& r_n) const
      {
        // dimensions
        static const unsigned int dim = IG::Geometry::dimension;
        static const unsigned int dimw = IG::Geometry::dimensionworld;

        // subspaces
        static_assert
          ((LFSV::CHILDREN == 2), "You seem to use the wrong function space for StokesDGVectorFEM");

        typedef typename LFSV::template Child<VBLOCK>::Type LFSV_V;
        const LFSV_V& lfsv_v_s = lfsv_s.template child<VBLOCK>();
        const LFSV_V& lfsu_v_s = lfsu_s.template child<VBLOCK>();
        const LFSV_V& lfsv_v_n = lfsv_n.template child<VBLOCK>();
        const LFSV_V& lfsu_v_n = lfsu_n.template child<VBLOCK>();

        typedef typename LFSV::template Child<PBLOCK>::Type LFSV_P;
        const LFSV_P& lfsv_p_s = lfsv_s.template child<PBLOCK>();
        const LFSV_P& lfsu_p_s = lfsu_s.template child<PBLOCK>();
        const LFSV_P& lfsv_p_n = lfsv_n.template child<PBLOCK>();
        const LFSV_P& lfsu_p_n = lfsu_n.template child<PBLOCK>();

        // domain and range field type
        typedef FiniteElementInterfaceSwitch<typename LFSV_V::Traits::FiniteElementType > FESwitch_V;
        typedef BasisInterfaceSwitch<typename FESwitch_V::Basis > BasisSwitch_V;
        typedef VectorBasisInterfaceSwitch<typename FESwitch_V::Basis > VectorBasisSwitch_V;
        typedef FiniteElementInterfaceSwitch<typename LFSV_P::Traits::FiniteElementType > FESwitch_P;
        typedef BasisInterfaceSwitch<typename FESwitch_P::Basis > BasisSwitch_P;
        typedef typename BasisSwitch_V::DomainField DF;
        typedef typename BasisSwitch_V::RangeField RF;
        typedef typename BasisSwitch_V::Range Range_V;
        typedef typename BasisSwitch_P::Range Range_P;
        typedef typename LFSV::Traits::SizeType size_type;

        // select quadrature rule
        Dune::GeometryType gt = ig.geometry().type();
        const Dune::QuadratureRule<DF,dim-1>& rule
          = Dune::QuadratureRules<DF,dim-1>::rule(gt,qorder);

        const typename IG::EntityPointer self = ig.inside();
        const typename IG::EntityPointer neighbor = ig.outside();

        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); it!=rule.end(); ++it){

          // position of quadrature point in adjacent elements
          const Dune::FieldVector<DF,dim> local_s = ig.geometryInInside().global(it->position());
          const Dune::FieldVector<DF,dim> local_n = ig.geometryInOutside().global(it->position());
          const Dune::FieldVector<DF,dim> global = ig.geometry().global(it->position());
          const RF weight = it->weight() * ig.geometry().integrationElement(it->position());

          // values of velocity shape functions
          std::vector<Range_V> phi_v_s(lfsv_v_s.size());
          std::vector<Range_V> phi_v_n(lfsv_v_n.size());
          FESwitch_V::basis(lfsv_v_s.finiteElement()).evaluateFunction(local_s,phi_v_s);
          FESwitch_V::basis(lfsv_v_n.finiteElement()).evaluateFunction(local_n,phi_v_n);

          // values of velocity shape functions
          std::vector<Range_P> phi_p_s(lfsv_p_s.size());
          std::vector<Range_P> phi_p_n(lfsv_p_n.size());
          FESwitch_P::basis(lfsv_p_s.finiteElement()).evaluateFunction(local_s,phi_p_s);
          FESwitch_P::basis(lfsv_p_n.finiteElement()).evaluateFunction(local_n,phi_p_n);

          // evaluate jacobian of basis functions on reference element
          std::vector<Dune::FieldMatrix<RF,dim,dim> > jac_phi_v_s(lfsu_v_s.size());
          std::vector<Dune::FieldMatrix<RF,dim,dim> > jac_phi_v_n(lfsu_v_n.size());
          VectorBasisSwitch_V::jacobian
            (FESwitch_V::basis(lfsv_v_s.finiteElement()), ig.inside()->geometry(), local_s, jac_phi_v_s);
          VectorBasisSwitch_V::jacobian
            (FESwitch_V::basis(lfsv_v_n.finiteElement()), ig.outside()->geometry(), local_n, jac_phi_v_n);

          // compute divergence of test functions
          std::vector<RF> div_phi_v_s(lfsv_v_s.size(),0.0);
          std::vector<RF> div_phi_v_n(lfsv_v_s.size(),0.0);
          for (size_type d=0; d<dim; d++){
            for (size_type i=0; i<lfsv_v_s.size(); i++)
              div_phi_v_s[i] += jac_phi_v_s[i][d][d];
            for (size_type i=0; i<lfsv_v_n.size(); i++)
              div_phi_v_n[i] += jac_phi_v_n[i][d][d];
          }

          // compute velocity value, jacobian, and divergence
          Range_V val_u_s(0.0);
          Range_V val_u_n(0.0);
          Dune::FieldMatrix<RF,dim,dim> jac_u_s(0.0);
          Dune::FieldMatrix<RF,dim,dim> jac_u_n(0.0);
          RF div_u_s(0.0);
          RF div_u_n(0.0);
          for (size_type i=0; i<lfsu_v_s.size(); i++){
            val_u_s.axpy(x_s(lfsu_v_s,i),phi_v_s[i]);
            jac_u_s.axpy(x_s(lfsu_v_s,i),jac_phi_v_s[i]);
            div_u_s += x_s(lfsu_v_s,i) * div_phi_v_s[i];
          }
          for (size_type i=0; i<lfsu_v_n.size(); i++){
            val_u_n.axpy(x_n(lfsu_v_n,i),phi_v_n[i]);
            jac_u_n.axpy(x_n(lfsu_v_n,i),jac_phi_v_n[i]);
            div_u_n += x_n(lfsu_v_n,i) * div_phi_v_n[i];
          }

          // compute pressure value
          Range_P val_p_s(0.0);
          Range_P val_p_n(0.0);
          for (size_type i=0; i<lfsu_p_s.size(); i++)
            val_p_s.axpy(x_s(lfsu_p_s,i),phi_p_s[i]);
          for (size_type i=0; i<lfsu_p_n.size(); i++)
            val_p_n.axpy(x_n(lfsu_p_n,i),phi_p_n[i]);

          // face normal vector
          const Dune::FieldVector<DF,dimw> normal = ig.unitOuterNormal(it->position());

          // compute jump in solution
          const Dune::FieldVector<DF,dimw> jump = val_u_s - val_u_n;

          {// update residual (flux term)

            // compute flux
            Dune::FieldVector<DF,dimw> flux(0.0);
            add_compute_flux(jac_u_s,normal,flux);
            add_compute_flux(jac_u_n,normal,flux);
            flux *= 0.5;

            const RF factor = weight * mu;
            for (size_t i=0; i<lfsv_v_s.size(); i++)
              r_s.accumulate(lfsv_v_s, i, -(flux * phi_v_s[i]) * factor);
            for (size_t i=0; i<lfsv_v_n.size(); i++)
              r_n.accumulate(lfsv_v_n, i,  (flux * phi_v_n[i]) * factor);
          }

          {// update residual (symmetry term)
            const RF factor = weight * mu;

            for (size_t i=0; i<lfsv_v_s.size(); i++){
              Dune::FieldVector<DF,dimw> flux(0.0);
              add_compute_flux(jac_phi_v_s[i],normal,flux);
              r_s.accumulate(lfsv_v_s, i, epsilon * 0.5 * (flux * jump) * factor);
            }
            for (size_t i=0; i<lfsv_v_n.size(); i++){
              Dune::FieldVector<DF,dimw> flux(0.0);
              add_compute_flux(jac_phi_v_n[i],normal,flux);
              r_n.accumulate(lfsv_v_n, i, epsilon * 0.5 * (flux * jump) * factor);
            }
          }

          {// interior penalty
            const RF factor = weight;
            const RF gamma = ip_factor.getFaceIP(ig);
            for (size_t i=0; i<lfsv_v_s.size(); i++)
              r_s.accumulate(lfsv_v_s,i, gamma * (jump * phi_v_s[i]) * factor);
            for (size_t i=0; i<lfsv_v_n.size(); i++)
              r_n.accumulate(lfsv_v_n,i, - gamma * (jump * phi_v_n[i]) * factor);
          }

          {// pressure and incompressibility
            const RF factor = weight;
            const RF val_p = 0.5 * (val_p_s + val_p_n);
            for (size_t i=0; i<lfsv_v_s.size(); i++)
              r_s.accumulate(lfsv_v_s, i, val_p * (phi_v_s[i] * normal) * factor);
            for (size_t i=0; i<lfsv_v_n.size(); i++)
              r_n.accumulate(lfsv_v_n, i, - val_p * (phi_v_n[i] * normal) * factor);

            for (size_t i=0; i<lfsv_p_s.size(); i++)
              r_s.accumulate(lfsv_p_s, i, 0.5 * phi_p_s[i] * (jump*normal) * factor);
            for (size_t i=0; i<lfsv_p_n.size(); i++)
              r_n.accumulate(lfsv_p_n, i, 0.5 * phi_p_n[i] * (jump*normal) * factor);
          }

        } // it - quadrature

      }

      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_boundary (const IG& ig,
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           R& r_s) const
      {
        // dimensions
        static const unsigned int dim = IG::Geometry::dimension;
        static const unsigned int dimw = IG::Geometry::dimensionworld;

        // subspaces
        static_assert
          ((LFSV::CHILDREN == 2), "You seem to use the wrong function space for StokesDGVectorFEM");

        typedef typename LFSV::template Child<VBLOCK>::Type LFSV_V;
        const LFSV_V& lfsv_v_s = lfsv_s.template child<VBLOCK>();
        const LFSV_V& lfsu_v_s = lfsu_s.template child<VBLOCK>();

        typedef typename LFSV::template Child<PBLOCK>::Type LFSV_P;
        const LFSV_P& lfsv_p_s = lfsv_s.template child<PBLOCK>();
        const LFSV_P& lfsu_p_s = lfsu_s.template child<PBLOCK>();

        // domain and range field type
        typedef FiniteElementInterfaceSwitch<typename LFSV_V::Traits::FiniteElementType > FESwitch_V;
        typedef BasisInterfaceSwitch<typename FESwitch_V::Basis > BasisSwitch_V;
        typedef VectorBasisInterfaceSwitch<typename FESwitch_V::Basis > VectorBasisSwitch_V;
        typedef FiniteElementInterfaceSwitch<typename LFSV_P::Traits::FiniteElementType > FESwitch_P;
        typedef BasisInterfaceSwitch<typename FESwitch_P::Basis > BasisSwitch_P;
        typedef typename BasisSwitch_V::DomainField DF;
        typedef typename BasisSwitch_V::RangeField RF;
        typedef typename BasisSwitch_V::Range Range_V;
        typedef typename BasisSwitch_P::Range Range_P;
        typedef typename LFSV::Traits::SizeType size_type;

        // select quadrature rule
        Dune::GeometryType gt = ig.geometry().type();
        const Dune::QuadratureRule<DF,dim-1>& rule
          = Dune::QuadratureRules<DF,dim-1>::rule(gt,qorder);

        const typename IG::EntityPointer self = ig.inside();

        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); it!=rule.end(); ++it){

          // position of quadrature point in adjacent elements
          const Dune::FieldVector<DF,dim> local_s = ig.geometryInInside().global(it->position());
          const Dune::FieldVector<DF,dim> global = ig.geometry().global(it->position());
          const RF weight = it->weight() * ig.geometry().integrationElement(it->position());

          // evaluate boundary condition type
          typename B::Traits::RangeType bctype;
          b.evaluate(ig,it->position(),bctype);

          // values of velocity shape functions
          std::vector<Range_V> phi_v_s(lfsv_v_s.size());
          FESwitch_V::basis(lfsv_v_s.finiteElement()).evaluateFunction(local_s,phi_v_s);

          // values of velocity shape functions
          std::vector<Range_P> phi_p_s(lfsv_p_s.size());
          FESwitch_P::basis(lfsv_p_s.finiteElement()).evaluateFunction(local_s,phi_p_s);

          // evaluate jacobian of basis functions on reference element
          std::vector<Dune::FieldMatrix<RF,dim,dim> > jac_phi_v_s(lfsu_v_s.size());
          VectorBasisSwitch_V::jacobian
            (FESwitch_V::basis(lfsv_v_s.finiteElement()), ig.inside()->geometry(), local_s, jac_phi_v_s);

          // compute divergence of test functions
          std::vector<RF> div_phi_v_s(lfsv_v_s.size(),0.0);

          for (size_type d=0; d<dim; d++){
            for (size_type i=0; i<lfsv_v_s.size(); i++)
              div_phi_v_s[i] += jac_phi_v_s[i][d][d];
          }

          // compute velocity value, jacobian, and divergence
          Range_V val_u_s(0.0);
          Dune::FieldMatrix<RF,dim,dim> jac_u_s(0.0);
          RF div_u_s(0.0);
          for (size_type i=0; i<lfsu_v_s.size(); i++){
            val_u_s.axpy(x_s(lfsu_v_s,i),phi_v_s[i]);
            jac_u_s.axpy(x_s(lfsu_v_s,i),jac_phi_v_s[i]);
            div_u_s += x_s(lfsu_v_s,i) * div_phi_v_s[i];
          }

          // compute pressure value
          Range_P val_p_s(0.0);
          for (size_type i=0; i<lfsu_p_s.size(); i++)
            val_p_s.axpy(x_s(lfsu_p_s,i),phi_p_s[i]);

          // face normal vector
          const Dune::FieldVector<DF,dimw> normal = ig.unitOuterNormal(it->position());

          if (bctype == BC::VelocityDirichlet) {

            // Compute jump relative to Dirichlet value
            typename V::Traits::RangeType u0;
            v.evaluateGlobal(global,u0);
            const Dune::FieldVector<DF,dimw> jump = val_u_s - u0;

            {// update residual (flux term)

              // compute flux
              Dune::FieldVector<DF,dimw> flux(0.0);
              add_compute_flux(jac_u_s,normal,flux);

              const RF factor = weight * mu;
              for (size_t i=0; i<lfsv_v_s.size(); i++)
                r_s.accumulate(lfsv_v_s, i, -(flux * phi_v_s[i]) * factor);
            }

            {// update residual (symmetry term)
              const RF factor = weight * mu;

              for (size_t i=0; i<lfsv_v_s.size(); i++){
                Dune::FieldVector<DF,dimw> flux(0.0);
                add_compute_flux(jac_phi_v_s[i],normal,flux);
                r_s.accumulate(lfsv_v_s, i, epsilon * 0.5 * (flux * jump) * factor);
              }
            }

            {// interior penalty
              const RF factor = weight;
              const RF gamma = ip_factor.getFaceIP(ig);
              for (size_t i=0; i<lfsv_v_s.size(); i++)
                r_s.accumulate(lfsv_v_s,i, gamma * (jump * phi_v_s[i]) * factor);
            }

            {// pressure and incompressibility
              const RF factor = weight;
              const RF val_p = val_p_s;
              for (size_t i=0; i<lfsv_v_s.size(); i++)
                r_s.accumulate(lfsv_v_s, i, val_p * (phi_v_s[i] * normal) * factor);
              for (size_t i=0; i<lfsv_p_s.size(); i++)
                r_s.accumulate(lfsv_p_s, i, 0.5 * phi_p_s[i] * (jump*normal) * factor);
            }

          } // DirichletVelocity

          if (bctype == BC::StressNeumann){
            typename P::Traits::RangeType p0;
            p.evaluateGlobal(global,p0);

            for (size_t i=0; i<lfsv_v_s.size(); i++){
              const RF val = p0 * (normal*phi_v_s[i]) * weight;
              r_s.accumulate(lfsv_v_s, i, val);
            }
          } // PressureDirichlet

        } // it - quadrature

      }

    private:

      template<class M, class RF>
      static void contraction(const M & a, const M & b, RF & v)
      {
        v = 0;
        const int an = a.N();
        const int am = a.M();
        for(int r=0; r<an; ++r)
          for(int c=0; c<am; ++c)
            v += a[r][c] * b[r][c];
      }

      template<class DU, class R>
      static void add_compute_flux(const DU & du, const R & n, R & result)
      {
        const int N = du.N();
        const int M = du.M();
        for(int r=0; r<N; ++r)
          for(int c=0; c<M; ++c)
            result[r] += du[r][c] * n[c];
      }

      const F& f;
      const B& b;
      const V& v;
      const P& p;
      // values for NIPG / NIPG
      int    epsilon;
      int    qorder;
      // physical parameters
      double mu;
      const IP & ip_factor;
    };

    //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif
