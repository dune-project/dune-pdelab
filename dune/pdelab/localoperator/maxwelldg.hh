// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_LOCALOPERATOR_MAXWELLDG_HH
#define DUNE_PDELAB_LOCALOPERATOR_MAXWELLDG_HH

#include<vector>

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/indices.hh>

#include <dune/geometry/referenceelements.hh>

#include<dune/pdelab/common/quadraturerules.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/geometrywrapper.hh>
#include<dune/pdelab/finiteelement/localbasiscache.hh>
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/idefault.hh>
#include<dune/pdelab/localoperator/pattern.hh>

#include"maxwellparameter.hh"

namespace Dune {
  namespace PDELab {


    template<int dim>
    class MaxwellEigenvectors
    {};

    /** \brief provide matrix which contains column-wise the eigenvectors of
               maxwell problem

        \tparam dim the space dimension
    */
    template<>
    class MaxwellEigenvectors<3>
    {
      enum { dim = 3 };
    public:

      /**
         Returns eigenvalues in the order

         s, s, -s, -s, 0, 0  (s = 1/sqrt(\mu \epsilon)

         \param eps permittivity
         \param mu permeability
         \param RT matrix to be filled
      */
      template<typename T1, typename T2, typename T3>
      static void eigenvalues (T1 eps, T1 mu, const Dune::FieldVector<T2,2*dim>& e)
      {
        T1 s = 1.0/sqrt(mu*eps); //speed of light s = 1/sqrt(\mu \epsilon)
        e[0] = s;
        e[1] = s;
        e[2] = -s;
        e[3] = -s;
        e[4] = 0;
        e[5] = 0;
      }

      /**
         Returns a matrix with columnwise the eigenvectors. They correspond to
         the eigenvalues in the order

         s, s, -s, -s, 0, 0  (s = 1/sqrt(\mu \epsilon)

         \param eps permittivity
         \param mu permeability
         \param n unit outer normal vector
         \param R matrix to be filled
      */
      template<typename T1, typename T2, typename T3>
      static void eigenvectors (T1 eps, T1 mu, const Dune::FieldVector<T2,dim>& n, Dune::FieldMatrix<T3,2*dim,2*dim>& R)
      {
        T1 a=n[0], b=n[1], c=n[2];

        Dune::FieldVector<T2,dim> re, im;
        if (std::abs(c)<0.5)
          {
            re[0]=a*c; re[1]=b*c; re[2]=c*c-1;
            im[0]=-b;  im[1]=a;   im[2]=0;
          }
        else
          {
            re[0]=a*b; re[1]=b*b-1; re[2]=b*c;
            im[0]=c;  im[1]=0.0;   im[2]=-a;
          }

        // \lambda_0,1 = s
        R[0][0] =  re[0];   R[0][1] =  -im[0];
        R[1][0] =  re[1];   R[1][1] =  -im[1];
        R[2][0] =  re[2];   R[2][1] =  -im[2];
        R[3][0] =  im[0];   R[3][1] =  re[0];
        R[4][0] =  im[1];   R[4][1] =  re[1];
        R[5][0] =  im[2];   R[5][1] =  re[2];

        // \lambda_2,3 = -s
        R[0][2] =  im[0];  R[0][3] =  re[0];
        R[1][2] =  im[1];  R[1][3] =  re[1];
        R[2][2] =  im[2];  R[2][3] =  re[2];
        R[3][2] =  re[0];  R[3][3] =  -im[0];
        R[4][2] =  re[1];  R[4][3] =  -im[1];
        R[5][2] =  re[2];  R[5][3] =  -im[2];

        // \lambda_4,5 = 0
        R[0][4] =   a;  R[0][5] =   0;
        R[1][4] =   b;  R[1][5] =   0;
        R[2][4] =   c;  R[2][5] =   0;
        R[3][4] =   0;  R[3][5] =   a;
        R[4][4] =   0;  R[4][5] =   b;
        R[5][4] =   0;  R[5][5] =   c;

        // apply scaling
        T1 weps=sqrt(eps);
        T1 wmu=sqrt(mu);
        for (std::size_t i=0; i<3; i++)
          for (std::size_t j=0; j<6; j++)
            R[i][j] *= weps;
        for (std::size_t i=3; i<6; i++)
          for (std::size_t j=0; j<6; j++)
            R[i][j] *= wmu;

        return;

        // if (std::abs(std::abs(c)-1)<1e-10)
        //   {
        //     if (c>0)
        //       {
        //         // \lambda_0,1 = s
        //         R[0][0] =  0;   R[0][1] =  1;
        //         R[1][0] =  -1;  R[1][1] =  0;
        //         R[2][0] =  0;   R[2][1] =  0;
        //         R[3][0] =  1;   R[3][1] =  0;
        //         R[4][0] =  0;   R[4][1] =  1;
        //         R[5][0] =  0;   R[5][1] =  0;

        //         // \lambda_2,3 = -s
        //         R[0][2] =  -1; R[0][3] =  0;
        //         R[1][2] =  0;  R[1][3] =  1;
        //         R[2][2] =  0;  R[2][3] =  0;
        //         R[3][2] =  0;  R[3][3] =  1;
        //         R[4][2] =  1;  R[4][3] =  0;
        //         R[5][2] =  0;  R[5][3] =  0;
        //       }
        //     else
        //       {
        //         // \lambda_0,1 = s
        //         R[0][0] =  -1; R[0][1] =  0;
        //         R[1][0] =  0;  R[1][1] =  1;
        //         R[2][0] =  0;  R[2][1] =  0;
        //         R[3][0] =  0;  R[3][1] =  1;
        //         R[4][0] =  1;  R[4][1] =  0;
        //         R[5][0] =  0;  R[5][1] =  0;

        //         // \lambda_2,3 = -s
        //         R[0][2] =  0;   R[0][3] =  1;
        //         R[1][2] =  -1;  R[1][3] =  0;
        //         R[2][2] =  0;   R[2][3] =  0;
        //         R[3][2] =  1;   R[3][3] =  0;
        //         R[4][2] =  0;   R[4][3] =  1;
        //         R[5][2] =  0;   R[5][3] =  0;
        //       }
        //   }
        // else if (std::abs(std::abs(b)-1)<1e-10)
        //   {
        //     if (b>0)
        //       {
        //         // \lambda_0,1 = s
        //         R[0][0] =  -1;  R[0][1] =  0;
        //         R[1][0] =  0;   R[1][1] =  0;
        //         R[2][0] =  0;   R[2][1] =  1;
        //         R[3][0] =  0;   R[3][1] =  1;
        //         R[4][0] =  0;   R[4][1] =  0;
        //         R[5][0] =  1;   R[5][1] =  0;

        //         // \lambda_2,3 = -s
        //         R[0][2] =  0;  R[0][3] =  1;
        //         R[1][2] =  0;  R[1][3] =  0;
        //         R[2][2] =  -1; R[2][3] =  0;
        //         R[3][2] =  1;  R[3][3] =  0;
        //         R[4][2] =  0;  R[4][3] =  0;
        //         R[5][2] =  0;  R[5][3] =  1;
        //       }
        //     else
        //       {
        //         // \lambda_0,1 = s
        //         R[0][0] =  0;  R[0][1] =  1;
        //         R[1][0] =  0;  R[1][1] =  0;
        //         R[2][0] =  -1; R[2][1] =  0;
        //         R[3][0] =  1;  R[3][1] =  0;
        //         R[4][0] =  0;  R[4][1] =  0;
        //         R[5][0] =  0;  R[5][1] =  1;

        //         // \lambda_2,3 = -s
        //         R[0][2] =  -1;  R[0][3] =  0;
        //         R[1][2] =  0;   R[1][3] =  0;
        //         R[2][2] =  0;   R[2][3] =  1;
        //         R[3][2] =  0;   R[3][3] =  1;
        //         R[4][2] =  0;   R[4][3] =  0;
        //         R[5][2] =  1;   R[5][3] =  0;
        //       }
        //   }
        // else if (std::abs(std::abs(a)-1)<1e-10)
        //   {
        //     if (a>0)
        //       {
        //         // \lambda_0,1 = s
        //         R[0][0] =  0;   R[0][1] =  0;
        //         R[1][0] =  0;   R[1][1] =  1;
        //         R[2][0] =  -1;  R[2][1] =  0;
        //         R[3][0] =  0;   R[3][1] =  0;
        //         R[4][0] =  1;   R[4][1] =  0;
        //         R[5][0] =  0;   R[5][1] =  1;

        //         // \lambda_2,3 = -s
        //         R[0][2] =  0;  R[0][3] =  0;
        //         R[1][2] =  -1; R[1][3] =  0;
        //         R[2][2] =  0;  R[2][3] =  1;
        //         R[3][2] =  0;  R[3][3] =  0;
        //         R[4][2] =  0;  R[4][3] =  1;
        //         R[5][2] =  1;  R[5][3] =  0;
        //       }
        //     else
        //       {
        //         // \lambda_0,1 = s
        //         R[0][0] =  0;  R[0][1] =  0;
        //         R[1][0] =  -1; R[1][1] =  0;
        //         R[2][0] =  0;  R[2][1] =  1;
        //         R[3][0] =  0;  R[3][1] =  0;
        //         R[4][0] =  0;  R[4][1] =  1;
        //         R[5][0] =  1;  R[5][1] =  0;

        //         // \lambda_2,3 = -s
        //         R[0][2] =  0;   R[0][3] =  0;
        //         R[1][2] =  0;   R[1][3] =  1;
        //         R[2][2] =  -1;  R[2][3] =  0;
        //         R[3][2] =  0;   R[3][3] =  0;
        //         R[4][2] =  1;   R[4][3] =  0;
        //         R[5][2] =  0;   R[5][3] =  1;
        //       }
        //   }
        // else
        //   {
        //     DUNE_THROW(Dune::Exception,"need axiparallel grids for now!");

        //     // \lambda_0,1 = s
        //     R[0][0] =  b;          R[0][1] =  -(b*b+c*c);
        //     R[1][0] =  -a;         R[1][1] =  a*b;
        //     R[2][0] =  0;          R[2][1] =  a*c;
        //     R[3][0] =  a*c;        R[3][1] =  0;
        //     R[4][0] =  b*c;        R[4][1] =  -c;
        //     R[5][0] =  -(a*a+b*b); R[5][1] =  b;

        //     // \lambda_2,3 = -s
        //     R[0][2] =  -b;         R[0][3] =  -(b*b+c*c);
        //     R[1][2] =  a;          R[1][3] =  a*b;
        //     R[2][2] =  0;          R[2][3] =  a*c;
        //     R[3][2] =  a*c;        R[3][3] =  0;
        //     R[4][2] =  b*c;        R[4][3] =  c;
        //     R[5][2] =  -(a*a+b*b); R[5][3] =  -b;
        //   }

        // // \lambda_4,5 = 0
        // R[0][4] =   0;  R[0][5] =   a;
        // R[1][4] =   0;  R[1][5] =   b;
        // R[2][4] =   0;  R[2][5] =   c;
        // R[3][4] =   a;  R[3][5] =   0;
        // R[4][4] =   b;  R[4][5] =   0;
        // R[5][4] =   c;  R[5][5] =   0;

        // // apply scaling
        // T1 weps=sqrt(eps);
        // T1 wmu=sqrt(mu);
        // for (std::size_t i=0; i<3; i++)
        //   for (std::size_t j=0; j<6; j++)
        //     R[i][j] *= weps;
        // for (std::size_t i=3; i<6; i++)
        //   for (std::size_t j=0; j<6; j++)
        //     R[i][j] *= wmu;

        //std::cout << R << std::endl;

      }
    };

    /** \brief Spatial local operator for discontinuous Galerkin method for
               Maxwells Equations

        \f{align*}{
          - \nabla \times (\mu^{-1} B) + (\frac\sigma\epsilon) D &= j \\
          + \nabla \times (\epsilon^{-1} D)                      &= 0
        \f}

        with the state vector \f$u=(D,B)^T\f$ having 2*dim components

        - Assumes that the local function space is a power space
          with 2*dim identical components.
        - Actually assumes dim==3 in the flux computation
        - Assumes Galerkin method, i.e. U=V

        \tparam T parameter class
        \tparam FEM Finite Element Map needed to select the cache

        \warning This operator cannot deal with spatially varying
                 \f$\epsilon\f$ or \f$\mu\f$.  You will get a result, but it
                 will contain spurious reflections.  This is a known bug, see
                 [issue
                 #41](https://gitlab.dune-project.org/pdelab/dune-pdelab/issues/41).
    */
    template<typename T, typename FEM>
    class DGMaxwellSpatialOperator :
      public NumericalJacobianApplyVolume<DGMaxwellSpatialOperator<T,FEM> >,
      public NumericalJacobianVolume<DGMaxwellSpatialOperator<T,FEM> >,
      public NumericalJacobianApplySkeleton<DGMaxwellSpatialOperator<T,FEM> >,
      public NumericalJacobianSkeleton<DGMaxwellSpatialOperator<T,FEM> >,
      public NumericalJacobianApplyBoundary<DGMaxwellSpatialOperator<T,FEM> >,
      public NumericalJacobianBoundary<DGMaxwellSpatialOperator<T,FEM> >,
      public FullSkeletonPattern,
      public FullVolumePattern,
      public LocalOperatorDefaultFlags,
      public InstationaryLocalOperatorDefaultMethods<typename T::Traits::RangeFieldType>
    {
      enum { dim = T::Traits::GridViewType::dimension };

    public:
      // pattern assembly flags
      enum { doPatternVolume = true };
      enum { doPatternSkeleton = true };

      // residual assembly flags
      enum { doAlphaVolume  = true };
      enum { doAlphaSkeleton  = true };
      enum { doAlphaBoundary  = true };
      enum { doLambdaVolume  = true };

      // ! constructor
      DGMaxwellSpatialOperator (T& param_, int overintegration_=0)
        : param(param_), overintegration(overintegration_), cache(20)
      {
      }

      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
        // Define types
        using namespace Indices;
        using DGSpace = TypeTree::Child<LFSV,0>;
        using RF = typename DGSpace::Traits::FiniteElementType::Traits
          ::LocalBasisType::Traits::RangeFieldType;
        using size_type = typename DGSpace::Traits::SizeType;

        // paranoia check number of number of components
        static_assert(TypeTree::StaticDegree<LFSV>::value == dim*2,
                      "need exactly dim*2 components!");

        // get local function space that is identical for all components
        const auto& dgspace = child(lfsv,_0);

        // Reference to cell
        const auto& cell = eg.entity();

        // Get geometry
        auto geo = eg.geometry();

        // evaluate parameters (assumed constant per element)
        auto ref_el = referenceElement(geo);
        auto localcenter = ref_el.position(0,0);
        auto mu = param.mu(cell,localcenter);
        auto eps = param.eps(cell,localcenter);
        auto sigma = param.sigma(cell,localcenter);
        auto muinv = 1.0/mu;
        auto epsinv = 1.0/eps;

        //std::cout << "alpha_volume center=" << eg.geometry().center() << std::endl;

        // Transformation
        typename EG::Geometry::JacobianInverseTransposed jac;

        // Initialize vectors outside for loop
        Dune::FieldVector<RF,dim*2> u(0.0);
        std::vector<Dune::FieldVector<RF,dim> > gradphi(dgspace.size());

        // loop over quadrature points
        const int order = dgspace.finiteElement().localBasis().order();
        const int intorder = overintegration+2*order;
        for (const auto &qp : quadratureRule(geo,intorder))
          {
            // evaluate basis functions
            const auto& phi = cache[order].evaluateFunction
              (qp.position(), dgspace.finiteElement().localBasis());

            // evaluate state vector u
            for (size_type k=0; k<dim*2; k++){ // for all components
              u[k] = 0.0;
              for (size_type j=0; j<dgspace.size(); j++) // for all basis functions
                u[k] += x(lfsv.child(k),j)*phi[j];
            }
            //std::cout << "  u at " << qp.position() << " : " << u << std::endl;

            // evaluate gradient of basis functions (we assume Galerkin method lfsu=lfsv)
            const auto& js = cache[order].evaluateJacobian
              (qp.position(), dgspace.finiteElement().localBasis());

            // compute global gradients
            jac = geo.jacobianInverseTransposed(qp.position());
            for (size_type i=0; i<dgspace.size(); i++)
              jac.mv(js[i][0],gradphi[i]);

            // integrate
            auto factor = qp.weight() * geo.integrationElement(qp.position());

            Dune::FieldMatrix<RF,dim*2,dim> F;
            static_assert(dim == 3, "Sorry, the following flux implementation "
                          "can only work for dim == 3");
            F[0][0] = 0;            F[0][1] = -muinv*u[5];  F[0][2] = muinv*u[4];
            F[1][0] = muinv*u[5];   F[1][1] = 0;            F[1][2] = -muinv*u[3];
            F[2][0] =-muinv*u[4];   F[2][1] = muinv*u[3];   F[2][2] = 0;
            F[3][0] = 0;            F[3][1] = epsinv*u[2];  F[3][2] = -epsinv*u[1];
            F[4][0] = -epsinv*u[2]; F[4][1] = 0;            F[4][2] = epsinv*u[0];
            F[5][0] = epsinv*u[1];  F[5][1] = -epsinv*u[0]; F[5][2] = 0;

            // for all components of the system
            for (size_type i=0; i<dim*2; i++)
              // for all test functions of this component
              for (size_type k=0; k<dgspace.size(); k++)
                // for all dimensions
                for (size_type j=0; j<dim; j++)
                  r.accumulate(lfsv.child(i),k,-F[i][j]*gradphi[k][j]*factor);

            // for the first half of the system
            for (size_type i=0; i<dim; i++)
              // for all test functions of this component
              for (size_type k=0; k<dgspace.size(); k++)
                r.accumulate(lfsv.child(i),k,(sigma/eps)*u[i]*phi[k]*factor);

            // std::cout << "  residual: ";
            // for (size_type i=0; i<r.size(); i++) std::cout << r[i] << " ";
            // std::cout << std::endl;
          }
      }

      // skeleton integral depending on test and ansatz functions
      // each face is only visited ONCE!
      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_skeleton (const IG& ig,
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                           R& r_s, R& r_n) const
      {
        using std::max;
        using std::sqrt;

        // Define types
        using namespace Indices;
        using DGSpace = TypeTree::Child<LFSV,0>;
        using DF = typename DGSpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType;
        using RF = typename DGSpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType;
        using size_type = typename DGSpace::Traits::SizeType;

        // get local function space that is identical for all components
        const auto& dgspace_s = child(lfsv_s,_0);
        const auto& dgspace_n = child(lfsv_n,_0);

        // References to inside and outside cells
        const auto& cell_inside = ig.inside();
        const auto& cell_outside = ig.outside();

        // Get geometries
        auto geo = ig.geometry();
        auto geo_inside = cell_inside.geometry();
        auto geo_outside = cell_outside.geometry();

       // Geometry of intersection in local coordinates of inside_cell and outside_cell
        auto geo_in_inside = ig.geometryInInside();
        auto geo_in_outside = ig.geometryInOutside();

        // evaluate speed of sound (assumed constant per element)
        auto ref_el_inside = referenceElement(geo_inside);
        auto ref_el_outside = referenceElement(geo_outside);
        auto local_inside = ref_el_inside.position(0,0);
        auto local_outside = ref_el_outside.position(0,0);
        // This is wrong -- this approach (with A- and A+) does not allow
        // position-dependent eps and mu, so we should not allow the parameter
        // class to specify them in a position-dependent manner.  See my
        // (Jorrit Fahlke) dissertation on how to do it right.
        auto mu_s = param.mu(cell_inside,local_inside);
        auto mu_n = param.mu(cell_outside,local_outside);
        auto eps_s = param.eps(cell_inside,local_inside);
        auto eps_n = param.eps(cell_outside,local_outside);
        //auto sigma_s = param.sigma(ig.inside(),local_inside);
        //auto sigma_n = param.sigma(ig.outside(),local_outside);

        // normal: assume faces are planar
        const auto& n_F = ig.centerUnitOuterNormal();

        // compute A+ (outgoing waves)
        Dune::FieldMatrix<DF,dim*2,dim*2> R_s;
        MaxwellEigenvectors<dim>::eigenvectors(eps_s,mu_s,n_F,R_s);
        Dune::FieldMatrix<DF,dim*2,dim*2> Dplus_s(0.0);
        Dplus_s[0][0] = 1.0/sqrt(eps_s*mu_s);
        Dplus_s[1][1] = 1.0/sqrt(eps_s*mu_s);
        Dune::FieldMatrix<DF,dim*2,dim*2> Aplus_s(R_s);
        Aplus_s.rightmultiply(Dplus_s);
        R_s.invert();
        Aplus_s.rightmultiply(R_s);

        // compute A- (incoming waves)
        Dune::FieldMatrix<DF,dim*2,dim*2> R_n;
        MaxwellEigenvectors<dim>::eigenvectors(eps_n,mu_n,n_F,R_n);
        Dune::FieldMatrix<DF,dim*2,dim*2> Dminus_n(0.0);
        Dminus_n[2][2] = -1.0/sqrt(eps_n*mu_n);
        Dminus_n[3][3] = -1.0/sqrt(eps_n*mu_n);
        Dune::FieldMatrix<DF,dim*2,dim*2> Aminus_n(R_n);
        Aminus_n.rightmultiply(Dminus_n);
        R_n.invert();
        Aminus_n.rightmultiply(R_n);

        // Initialize vectors outside for loop
        Dune::FieldVector<RF,dim*2> u_s(0.0);
        Dune::FieldVector<RF,dim*2> u_n(0.0);
        Dune::FieldVector<RF,dim*2> f(0.0);

        // std::cout << "alpha_skeleton center=" << ig.geometry().center() << std::endl;

        // loop over quadrature points
        const int order_s = dgspace_s.finiteElement().localBasis().order();
        const int order_n = dgspace_n.finiteElement().localBasis().order();
        const int intorder = overintegration+1+2*max(order_s,order_n);
        for (const auto& qp : quadratureRule(geo,intorder))
          {
            // position of quadrature point in local coordinates of elements
            const auto& iplocal_s = geo_in_inside.global(qp.position());
            const auto& iplocal_n = geo_in_outside.global(qp.position());

            // evaluate basis functions
            const auto& phi_s = cache[order_s].evaluateFunction(iplocal_s,dgspace_s.finiteElement().localBasis());
            const auto& phi_n = cache[order_n].evaluateFunction(iplocal_n,dgspace_n.finiteElement().localBasis());

            // evaluate u from inside and outside
            for (size_type i=0; i<dim*2; i++){ // for all components
              u_s[i] = 0.0;
              for (size_type k=0; k<dgspace_s.size(); k++) // for all basis functions
                u_s[i] += x_s(lfsv_s.child(i),k)*phi_s[k];
            }
            for (size_type i=0; i<dim*2; i++){ // for all components
              u_n[i] = 0.0;
              for (size_type k=0; k<dgspace_n.size(); k++) // for all basis functions
                u_n[i] += x_n(lfsv_n.child(i),k)*phi_n[k];
            }

            // compute numerical flux at integration point
            f = 0.0;
            Aplus_s.umv(u_s,f);
            // std::cout << "  after A_plus*u_s  " << f << std::endl;
            Aminus_n.umv(u_n,f);
            // std::cout << "  after A_minus*u_n " << f << std::endl;

            //std::cout << "f=" << f << " u_s=" << u_s << " u_n=" << u_n << std::endl;

            // integrate
            auto factor = qp.weight() * geo.integrationElement(qp.position());
            for (size_type k=0; k<dgspace_s.size(); k++) // loop over all vector-valued (!) basis functions (with identical components)
              for (size_type i=0; i<dim*2; i++) // loop over all components
                r_s.accumulate(lfsv_s.child(i),k,f[i]*phi_s[k]*factor);
            for (size_type k=0; k<dgspace_n.size(); k++) // loop over all vector-valued (!) basis functions (with identical components)
              for (size_type i=0; i<dim*2; i++) // loop over all components
                r_n.accumulate(lfsv_n.child(i),k,-f[i]*phi_n[k]*factor);
          }

        // std::cout << "  residual_s: ";
        // for (auto v : r_s) std::cout << v << " ";
        // std::cout << std::endl;
        // std::cout << "  residual_n: ";
        // for (auto v : r_n) std::cout << v << " ";
        // std::cout << std::endl;
      }

      // skeleton integral depending on test and ansatz functions
      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_boundary (const IG& ig,
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           R& r_s) const
      {
        // Define types
        using namespace Indices;
        using DGSpace = TypeTree::Child<LFSV,0>;
        using DF = typename DGSpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType;
        using RF = typename DGSpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType;
        using size_type = typename DGSpace::Traits::SizeType;

        // get local function space that is identical for all components
        const auto& dgspace_s = child(lfsv_s,_0);

        // References to inside cell
        const auto& cell_inside = ig.inside();

        // Get geometries
        auto geo = ig.geometry();
        auto geo_inside = cell_inside.geometry();

        // Geometry of intersection in local coordinates of inside_cell
        auto geo_in_inside = ig.geometryInInside();

        // evaluate speed of sound (assumed constant per element)
        auto ref_el_inside = referenceElement(geo_inside);
        auto local_inside = ref_el_inside.position(0,0);
        auto mu_s = param.mu(cell_inside,local_inside);
        auto eps_s = param.eps(cell_inside,local_inside);
        //auto sigma_s = param.sigma(ig.inside(),local_inside );

        // normal: assume faces are planar
        const auto& n_F = ig.centerUnitOuterNormal();

        // compute A+ (outgoing waves)
        Dune::FieldMatrix<DF,dim*2,dim*2> R_s;
        MaxwellEigenvectors<dim>::eigenvectors(eps_s,mu_s,n_F,R_s);
        Dune::FieldMatrix<DF,dim*2,dim*2> Dplus_s(0.0);
        Dplus_s[0][0] = 1.0/sqrt(eps_s*mu_s);
        Dplus_s[1][1] = 1.0/sqrt(eps_s*mu_s);
        Dune::FieldMatrix<DF,dim*2,dim*2> Aplus_s(R_s);
        Aplus_s.rightmultiply(Dplus_s);
        R_s.invert();
        Aplus_s.rightmultiply(R_s);

        // compute A- (incoming waves)
        Dune::FieldMatrix<DF,dim*2,dim*2> R_n;
        MaxwellEigenvectors<dim>::eigenvectors(eps_s,mu_s,n_F,R_n);
        Dune::FieldMatrix<DF,dim*2,dim*2> Dminus_n(0.0);
        Dminus_n[2][2] = -1.0/sqrt(eps_s*mu_s);
        Dminus_n[3][3] = -1.0/sqrt(eps_s*mu_s);
        Dune::FieldMatrix<DF,dim*2,dim*2> Aminus_n(R_n);
        Aminus_n.rightmultiply(Dminus_n);
        R_n.invert();
        Aminus_n.rightmultiply(R_n);

        // Initialize vectors outside for loop
        Dune::FieldVector<RF,dim*2> u_s(0.0);
        Dune::FieldVector<RF,dim*2> u_n;
        Dune::FieldVector<RF,dim*2> f(0.0);

        // std::cout << "alpha_boundary center=" << ig.geometry().center() << std::endl;

        // loop over quadrature points
        const int order_s = dgspace_s.finiteElement().localBasis().order();
        const int intorder = overintegration+1+2*order_s;
        for(const auto &qp : quadratureRule(geo,intorder))
          {
            // position of quadrature point in local coordinates of elements
            const auto& iplocal_s = geo_in_inside.global(qp.position());

            // evaluate basis functions
            const auto& phi_s = cache[order_s].evaluateFunction
              (iplocal_s,dgspace_s.finiteElement().localBasis());

            // evaluate u from inside and outside
            for (size_type i=0; i<dim*2; i++){ // for all components
              u_s[i] = 0.0;
              for (size_type k=0; k<dgspace_s.size(); k++) // for all basis functions
                u_s[i] += x_s(lfsv_s.child(i),k)*phi_s[k];
            }
            // std::cout << "  u_s " << u_s << std::endl;

            // evaluate boundary condition
            u_n = (param.g(ig.intersection(),qp.position(),u_s));
            // std::cout << "  u_n " << u_n << " bc: " << param.g(ig.intersection(),qp.position(),u_s) << std::endl;

            // compute numerical flux at integration point
            f = 0.0;
            Aplus_s.umv(u_s,f);
            // std::cout << "  after A_plus*u_s  " << f << std::endl;
            Aminus_n.umv(u_n,f);
            // std::cout << "  after A_minus*u_n " << f << std::endl;

            // integrate
            auto factor = qp.weight() * geo.integrationElement(qp.position());
            for (size_type k=0; k<dgspace_s.size(); k++) // loop over all vector-valued (!) basis functions (with identical components)
              for (size_type i=0; i<dim*2; i++) // loop over all components
                r_s.accumulate(lfsv_s.child(i),k,f[i]*phi_s[k]*factor);
          }

        // std::cout << "  residual_s: ";
        // for (size_type i=0; i<r_s.size(); i++) std::cout << r_s[i] << " ";
        // std::cout << std::endl;
      }

      // volume integral depending only on test functions
      template<typename EG, typename LFSV, typename R>
      void lambda_volume (const EG& eg, const LFSV& lfsv, R& r) const
      {
        // Define types
        using namespace Indices;
        using DGSpace = TypeTree::Child<LFSV,0>;
        using size_type = typename DGSpace::Traits::SizeType;

        // Get local function space that is identical for all components
        const auto& dgspace = child(lfsv,_0);

        // Reference to cell
        const auto& cell = eg.entity();

        // Get geometry
        auto geo = eg.geometry();

        // loop over quadrature points
        const int order_s = dgspace.finiteElement().localBasis().order();
        const int intorder = overintegration+2*order_s;
        for (const auto &qp : quadratureRule(geo,intorder))
          {
            // evaluate right hand side
            auto j = param.j(cell,qp.position());

            // evaluate basis functions
            const auto& phi = cache[order_s].evaluateFunction(qp.position(),dgspace.finiteElement().localBasis());

            // integrate
            auto factor = qp.weight() * geo.integrationElement(qp.position());
            for (size_type k=0; k<dim*2; k++) // for all components
              for (size_type i=0; i<dgspace.size(); i++) // for all test functions of this component
                r.accumulate(lfsv.child(k),i,-j[k]*phi[i]*factor);
          }
      }

      //! set time in parameter class
      void setTime (typename T::Traits::RangeFieldType t)
      {
      }

      //! to be called once before each time step
      void preStep (typename T::Traits::RangeFieldType time, typename T::Traits::RangeFieldType dt,
                    int stages)
      {
      }

      //! to be called once before each stage
      void preStage (typename T::Traits::RangeFieldType time, int r)
      {
      }

      //! to be called once at the end of each stage
      void postStage ()
      {
      }

      //! to be called once before each stage
      typename T::Traits::RangeFieldType suggestTimestep (typename T::Traits::RangeFieldType dt) const
      {
        return dt;
      }

    private:
      T& param;
      int overintegration;
      typedef typename FEM::Traits::FiniteElementType::Traits::LocalBasisType LocalBasisType;
      typedef Dune::PDELab::LocalBasisCache<LocalBasisType> Cache;
      std::vector<Cache> cache;
    };



    /** a local operator for the mass operator of a vector valued lfs (L_2 integral)
     *
     * \f{align*}{
     \int_\Omega uv dx
     * \f}
     *
     * This is pretty much the same as PowerL2.  Differences are:
     * \item this operator evaluates the child basis only once for child 0 an
     *       reuses that for the other components
     * \item this operator caches basis evaluations
     */
    template<typename T, typename FEM>
    class DGMaxwellTemporalOperator :
      public NumericalJacobianApplyVolume<DGMaxwellTemporalOperator<T,FEM> >,
      public LocalOperatorDefaultFlags,
      public InstationaryLocalOperatorDefaultMethods<typename T::Traits::RangeFieldType>
    {
      enum { dim = T::Traits::GridViewType::dimension };
    public:
      // pattern assembly flags
      enum { doPatternVolume = true };

      // residual assembly flags
      enum { doAlphaVolume = true };

      DGMaxwellTemporalOperator (T& param_, int overintegration_=0)
        : param(param_), overintegration(overintegration_), cache(20)
      {}

      // define sparsity pattern of operator representation
      template<typename LFSU, typename LFSV, typename LocalPattern>
      void pattern_volume (const LFSU& lfsu, const LFSV& lfsv,
                           LocalPattern& pattern) const
      {
        // paranoia check number of number of components
        static_assert(TypeTree::StaticDegree<LFSU>::value==TypeTree::StaticDegree<LFSV>::value, "need U=V!");
        static_assert(TypeTree::StaticDegree<LFSV>::value==dim*2, "need exactly dim*2 components!");

        for (size_t k=0; k<TypeTree::degree(lfsv); k++)
          for (size_t i=0; i<lfsv.child(k).size(); ++i)
            for (size_t j=0; j<lfsu.child(k).size(); ++j)
              pattern.addLink(lfsv.child(k),i,lfsu.child(k),j);
      }

      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
        // Define types
        using namespace Indices;
        using DGSpace = TypeTree::Child<LFSV,0>;
        using RF = typename DGSpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType;
        using size_type = typename DGSpace::Traits::SizeType;

        // get local function space that is identical for all components
        const auto& dgspace = child(lfsv,_0);

        // Get geometry
        auto geo = eg.geometry();

        // Initialize vectors outside for loop
        Dune::FieldVector<RF,dim*2> u(0.0);

        // loop over quadrature points
        const int order = dgspace.finiteElement().localBasis().order();
        const int intorder = overintegration+2*order;
        for (const auto& qp : quadratureRule(geo,intorder))
          {
            // evaluate basis functions
            const auto& phi = cache[order].evaluateFunction(qp.position(),dgspace.finiteElement().localBasis());

            // evaluate u
            for (size_type k=0; k<dim*2; k++){ // for all components
              u[k] = 0.0;
              for (size_type j=0; j<dgspace.size(); j++) // for all basis functions
                u[k] += x(lfsv.child(k),j)*phi[j];
            }

            // integrate
            auto factor = qp.weight() * geo.integrationElement(qp.position());
            for (size_type k=0; k<dim*2; k++) // for all components
              for (size_type i=0; i<dgspace.size(); i++) // for all test functions of this component
                r.accumulate(lfsv.child(k),i,u[k]*phi[i]*factor);
          }
      }

      // jacobian of volume term
      template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                            M& mat) const
      {
        // Define types
        using namespace Indices;
        using DGSpace = TypeTree::Child<LFSV,0>;
        using size_type = typename DGSpace::Traits::SizeType;

        // Get local function space that is identical for all components
        const auto& dgspace = child(lfsv,_0);

        // Get geometry
        auto geo = eg.geometry();

        // Loop over quadrature points
        const int order = dgspace.finiteElement().localBasis().order();
        const int intorder = overintegration+2*order;
        for (const auto& qp : quadratureRule(geo,intorder))
          {
            // Evaluate basis functions
            const auto& phi = cache[order].evaluateFunction(qp.position(),dgspace.finiteElement().localBasis());

            // Integrate
            auto factor = qp.weight() * geo.integrationElement(qp.position());
            for (size_type k=0; k<dim*2; k++) // for all components
              for (size_type i=0; i<dgspace.size(); i++) // for all test functions of this component
                for (size_type j=0; j<dgspace.size(); j++) // for all ansatz functions of this component
                  mat.accumulate(lfsv.child(k),i,lfsu.child(k),j,phi[j]*phi[i]*factor);
          }
      }

    private:
      T& param;
      int overintegration;
      typedef typename FEM::Traits::FiniteElementType::Traits::LocalBasisType LocalBasisType;
      typedef Dune::PDELab::LocalBasisCache<LocalBasisType> Cache;
      std::vector<Cache> cache;
    };

  }
}

#endif // DUNE_PDELAB_LOCALOPERATOR_MAXWELLDG_HH
