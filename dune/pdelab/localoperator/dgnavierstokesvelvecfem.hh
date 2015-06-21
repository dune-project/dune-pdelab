// -*- tab-width: 2; indent-tabs-mode: nil -*-
// vi: set et ts=2 sw=2 sts=2:

#ifndef DUNE_PDELAB_LOCALOPERATOR_DGNAVIERSTOKESVELVECFEM_HH
#define DUNE_PDELAB_LOCALOPERATOR_DGNAVIERSTOKESVELVECFEM_HH

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parametertree.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/localfunctions/common/interfaceswitch.hh>
#include <dune/pdelab/localoperator/idefault.hh>

#include <dune/pdelab/localoperator/defaultimp.hh>
#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/localoperator/dgnavierstokesparameter.hh>
#include <dune/pdelab/localoperator/navierstokesmass.hh>

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
    Basis, typename enable_if<Basis::Traits::dimDomain>::type
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

        \tparam PRM Parameter class for this local operator.

    */
    template<typename PRM>
    class DGNavierStokesVelVecFEM :
      public LocalOperatorDefaultFlags,
      public FullSkeletonPattern, public FullVolumePattern,
      public InstationaryLocalOperatorDefaultMethods<double>
    {
      typedef StokesBoundaryCondition BC;
      typedef typename PRM::Traits::RangeField RF;

      typedef InstationaryLocalOperatorDefaultMethods<double> InstatBase;
      typedef typename InstatBase::RealType Real;

      static const bool navier = PRM::assemble_navier;
      static const bool full_tensor = PRM::assemble_full_tensor;

    public :

      // pattern assembly flags
      enum { doPatternVolume = true };
      enum { doPatternSkeleton = true };

      // call the assembler for each face only once
      enum { doSkeletonTwoSided = false };

      // residual assembly flags
      enum { doAlphaVolume    = true };
      enum { doAlphaSkeleton  = true };
      enum { doAlphaBoundary  = true };
      enum { doLambdaVolume   = true };

      /** \brief Constructor

          \param [in] _prm                        Parameter class for this local operator
          \param [in] _superintegration_order     This number will be added to the order of
                                                  quadrature in every integration. It is
                                                  only needed, when one of the parameters (e.g
                                                  rho, mu) is not constant or the mappings from
                                                  the reference elements to the cells are
                                                  nonlinear. Boundary conditions are assumed to
                                                  have the same order as the corresponding
                                                  finite element.
      */
      DGNavierStokesVelVecFEM (PRM& _prm, int _superintegration_order=0) :
        prm(_prm), superintegration_order(_superintegration_order),
        current_dt(1.0)
      {}

      // Store current dt
      void preStep (RealType , RealType dt, int )
      {
        current_dt = dt;
      }

      // set time in parameter class
      void setTime(Real t)
      {
        InstatBase::setTime(t);
        prm.setTime(t);
      }

      //============================================
      // TODO
      // Finish implementing the alpha_*, jacobian_* methods.
      // What about the lambda_* methods?
      // Answer: Only implement the lambda_volume method for the source term.
      //         No lambda_boundary method.
      //============================================

      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
      } // end alpha_volume

      // jacobian of volume term
      template<typename EG, typename LFSU, typename X, typename LFSV,
               typename LocalMatrix>
      void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                            LocalMatrix& mat) const
      {
      } // end jacobian_volume

      // skeleton term, each face is only visited ONCE
      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_skeleton (const IG& ig,
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                           R& r_s, R& r_n) const
      {
      } // end alpha_skeleton

      // jacobian of skeleton term, each face is only visited ONCE
      template<typename IG, typename LFSU, typename X, typename LFSV,
               typename LocalMatrix>
      void jacobian_skeleton (const IG& ig,
                              const LFSU& lfsu_s, const X&, const LFSV& lfsv_s,
                              const LFSU& lfsu_n, const X&, const LFSV& lfsv_n,
                              LocalMatrix& mat_ss, LocalMatrix& mat_sn,
                              LocalMatrix& mat_ns, LocalMatrix& mat_nn) const
      {
      } // end jacobian_skeleton

      // boundary term
      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_boundary (const IG& ig,
                           const LFSU& lfsu, const X& x, const LFSV& lfsv,
                           R& r) const
      {
      } // end alpha_boundary

      // jacobian of boundary term
      template<typename IG, typename LFSU, typename X, typename LFSV,
               typename LocalMatrix>
      void jacobian_boundary (const IG& ig,
                              const LFSU& lfsu, const X& x, const LFSV& lfsv,
                              LocalMatrix& mat) const
      {
      } // end jacobian_boundary

    private :

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

      PRM& prm;
      const int superintegration_order;
      Real current_dt;
    }; // end class DGNavierStokesVelVecFEM

  } // end namespace PDELab
} // end namespace Dune
#endif
