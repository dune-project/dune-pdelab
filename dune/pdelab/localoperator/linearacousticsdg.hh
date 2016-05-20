// -*- tab-width: 8; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_LOCALOPERATOR_LINEARACOUSTICSDG_HH
#define DUNE_PDELAB_LOCALOPERATOR_LINEARACOUSTICSDG_HH

#include<vector>

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/geometry/referenceelements.hh>

#include<dune/pdelab/common/quadraturerules.hh>
#include<dune/pdelab/common/referenceelements.hh>
#include<dune/pdelab/common/geometrywrapper.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/idefault.hh>
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/finiteelement/localbasiscache.hh>

#include"linearacousticsparameter.hh"

namespace Dune {
  namespace PDELab {


    template<int dim>
    class LinearAcousticsEigenvectors
    {};

    /** \brief provide matrix which contains rowwise the eigenvectors of linear acoustics problem
        \tparam dim the space dimension
    */
    template<>
    class LinearAcousticsEigenvectors<1>
    {
      enum { dim = 1 };
    public:
      /**
         \param c speed of sound
         \param n unit outer normal vector
         \param RT matrix to be filled
      */
      template<typename T1, typename T2, typename T3>
      static void eigenvectors_transposed (T1 c, const Dune::FieldVector<T2,dim>& n, Dune::FieldMatrix<T3,dim+1,dim+1>& RT)
      {
        RT[0][0] =  1; RT[0][1] = c*n[0];
        RT[1][0] = -1; RT[1][1] = c*n[0];
      }

      template<typename T1, typename T2, typename T3>
      static void eigenvectors (T1 c, const Dune::FieldVector<T2,dim>& n, Dune::FieldMatrix<T3,dim+1,dim+1>& RT)
      {
        RT[0][0] =  1; RT[1][0] = c*n[0];
        RT[0][1] = -1; RT[1][1] = c*n[0];
      }

    };

    /** \brief provide matrix which contains rowwise the eigenvectors of linear acoustics problem
        \tparam dim the space dimension
    */
    template<>
    class LinearAcousticsEigenvectors<2>
    {
      enum { dim = 2 };
    public:
      /**
         \param c speed of sound
         \param n unit outer normal vector
         \param RT matrix to be filled
      */
      template<typename T1, typename T2, typename T3>
      static void eigenvectors_transposed (T1 c, const Dune::FieldVector<T2,dim>& n, Dune::FieldMatrix<T3,dim+1,dim+1>& RT)
      {
        RT[0][0] =  0; RT[0][1] =  -n[1];  RT[0][2] = n[0];
        RT[1][0] =  1; RT[1][1] = c*n[0];  RT[1][2] = c*n[1];
        RT[2][0] = -1; RT[2][1] = c*n[0];  RT[2][2] = c*n[1];
      }

      template<typename T1, typename T2, typename T3>
      static void eigenvectors (T1 c, const Dune::FieldVector<T2,dim>& n, Dune::FieldMatrix<T3,dim+1,dim+1>& RT)
      {
        RT[0][0] =  0; RT[1][0] =  -n[1];  RT[2][0] = n[0];
        RT[0][1] =  1; RT[1][1] = c*n[0];  RT[2][1] = c*n[1];
        RT[0][2] = -1; RT[1][2] = c*n[0];  RT[2][2] = c*n[1];
      }
    };

    /** \brief provide matrix which contains rowwise the eigenvectors of linear acoustics problem
        \tparam dim the space dimension
    */
    template<>
    class LinearAcousticsEigenvectors<3>
    {
      enum { dim = 3 };
    public:
      /**
         \param c speed of sound
         \param n unit outer normal vector
         \param RT matrix to be filled
      */
      template<typename T1, typename T2, typename T3>
      static void eigenvectors_transposed (T1 c, const Dune::FieldVector<T2,dim>& n, Dune::FieldMatrix<T3,dim+1,dim+1>& RT)
      {
        Dune::FieldVector<T2,dim> s;
        s[0] = n[1]-n[2]; s[1] = n[2]-n[0]; s[2] = n[0]-n[1];
        if (s.two_norm()<1e-5)
          {
            s[0] = n[1]+n[2]; s[1] = n[2]-n[0]; s[2] = -(n[0]+n[1]);
          }

        Dune::FieldVector<T2,dim> t; // compute cross product s * n
        t[0] = n[1]*s[2] - n[2]*s[1];
        t[1] = n[2]*s[0] - n[0]*s[2];
        t[2] = n[0]*s[1] - n[1]*s[0];

        RT[0][0] =  0;  RT[0][1] =   s[0];  RT[0][2] =   s[1];  RT[0][3] =   s[2];
        RT[1][0] =  0;  RT[1][1] =   t[0];  RT[1][2] =   t[1];  RT[1][3] =   t[2];
        RT[2][0] =  1;  RT[2][1] = c*n[0];  RT[2][2] = c*n[1];  RT[2][3] = c*n[2];
        RT[3][0] = -1;  RT[3][1] = c*n[0];  RT[3][2] = c*n[1];  RT[3][3] = c*n[2];
      }

      template<typename T1, typename T2, typename T3>
      static void eigenvectors (T1 c, const Dune::FieldVector<T2,dim>& n, Dune::FieldMatrix<T3,dim+1,dim+1>& RT)
      {
        Dune::FieldVector<T2,dim> s;
        s[0] = n[1]-n[2]; s[1] = n[2]-n[0]; s[2] = n[0]-n[1];
        if (s.two_norm()<1e-5)
          {
            s[0] = n[1]+n[2]; s[1] = n[2]-n[0]; s[2] = -(n[0]+n[1]);
          }

        Dune::FieldVector<T2,dim> t; // compute cross product s * n
        t[0] = n[1]*s[2] - n[2]*s[1];
        t[1] = n[2]*s[0] - n[0]*s[2];
        t[2] = n[0]*s[1] - n[1]*s[0];

        RT[0][0] =  0;  RT[1][0] =   s[0];  RT[2][0] =   s[1];  RT[3][0] =   s[2];
        RT[0][1] =  0;  RT[1][1] =   t[0];  RT[2][1] =   t[1];  RT[3][1] =   t[2];
        RT[0][2] =  1;  RT[1][2] = c*n[0];  RT[2][2] = c*n[1];  RT[3][2] = c*n[2];
        RT[0][3] = -1;  RT[1][3] = c*n[0];  RT[2][3] = c*n[1];  RT[3][3] = c*n[2];
      }
    };

    /** Spatial local operator for discontinuous Galerkin method for the equations
        of linear acoustics in conservative form:

        \nabla \cdot \{ w \}  = 0 in \Omega
        \nabla \{ c^2 \rho \} = 0 in \Omega
        A^-(x) u = A^-(x) g on \partial\Omega

        Where u = (\rho,w) is the solution with dim+1 components, w = \bar\rho v is the momentum.

        - Assumes that the local function space is a power space
        with dim+1 identical components.
        - Assumes Galerkin method, i.e. U=V

        \tparam T parameter class
        \tparam FEM Finite Element Map needed to select the cache
    */
    template<typename T, typename FEM>
    class DGLinearAcousticsSpatialOperator :
      public NumericalJacobianApplyVolume<DGLinearAcousticsSpatialOperator<T,FEM> >,
      public NumericalJacobianVolume<DGLinearAcousticsSpatialOperator<T,FEM> >,
      public NumericalJacobianApplySkeleton<DGLinearAcousticsSpatialOperator<T,FEM> >,
      public NumericalJacobianSkeleton<DGLinearAcousticsSpatialOperator<T,FEM> >,
      public NumericalJacobianApplyBoundary<DGLinearAcousticsSpatialOperator<T,FEM> >,
      public NumericalJacobianBoundary<DGLinearAcousticsSpatialOperator<T,FEM> >,
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
      DGLinearAcousticsSpatialOperator (T& param_, int overintegration_=0)
        : param(param_), overintegration(overintegration_), cache(20)
      {
      }

      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
        // Get types
        using namespace Indices;
        using DGSpace = TypeTree::Child<LFSV,_0>;
        using RF = typename DGSpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType;
        using size_type = typename DGSpace::Traits::SizeType;

        // paranoia check number of number of components
        if (TypeTree::degree(lfsv)!=dim+1)
          DUNE_THROW(Dune::Exception,"need exactly dim+1 components!");

        // get local function space that is identical for all components
        const auto& dgspace = child(lfsv,_0);

        // Reference to cell
        const auto& cell = eg.entity();

        // Get geometry
        auto geo = eg.geometry();

        // evaluate speed of sound (assumed constant per element)
        auto ref_el = referenceElement(geo);
        auto localcenter = ref_el.position(0,0);
        auto c2 = param.c(cell,localcenter);
        c2 = c2*c2; // square it

        // Transformation
        typename EG::Geometry::JacobianInverseTransposed jac;

        // Initialize vectors outside for loop
        Dune::FieldVector<RF,dim+1> u(0.0);
        std::vector<Dune::FieldVector<RF,dim> > gradphi(dgspace.size());

        // loop over quadrature points
        const int order = dgspace.finiteElement().localBasis().order();
        const int intorder = overintegration+2*order;
        for (const auto& ip : quadratureRule(geo,intorder))
          {
            // evaluate basis functions
            auto& phi = cache[order].evaluateFunction(ip.position(),dgspace.finiteElement().localBasis());

            // evaluate u
            u = 0.0;
            for (size_type k=0; k<=dim; k++) // for all components
              for (size_type j=0; j<dgspace.size(); j++) // for all basis functions
                u[k] += x(lfsv.child(k),j)*phi[j];
            // std::cout << "  u at " << ip.position() << " : " << u << std::endl;

            // evaluate gradient of basis functions (we assume Galerkin method lfsu=lfsv)
            auto& js = cache[order].evaluateJacobian(ip.position(),dgspace.finiteElement().localBasis());

            // compute global gradients
            jac = geo.jacobianInverseTransposed(ip.position());
            for (size_type i=0; i<dgspace.size(); i++)
              jac.mv(js[i][0],gradphi[i]);

            // integrate
            auto factor = ip.weight() * geo.integrationElement(ip.position());
            for (size_type k=0; k<dgspace.size(); k++) // loop over all vector-valued (!) basis functions (with identical components)
              {
                // component i=0
                for (size_type j=0; j<dim; j++)
                  r.accumulate(lfsv.child(0),k, - u[j+1]*gradphi[k][j]*factor);
                // components i=1...d
                for (size_type i=1; i<=dim; i++)
                  r.accumulate(lfsv.child(i),k, - c2*u[0]*gradphi[k][i-1]*factor);
              }
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
        // Get types
        using namespace Indices;
        using DGSpace = TypeTree::Child<LFSV,_0>;
        using DF = typename DGSpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType;
        using RF = typename DGSpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType;
        using size_type = typename DGSpace::Traits::SizeType;

        // Get local function space that is identical for all components
        const auto& dgspace_s = child(lfsv_s,_0);
        const auto& dgspace_n = child(lfsv_n,_0);

        // Normal: assume faces are planar
        const Dune::FieldVector<DF,dim> n_F = ig.centerUnitOuterNormal();

        // References to inside and outside cells
        const auto& cell_inside = ig.inside();
        const auto& cell_outside = ig.outside();

        // Get geometries
        auto geo = ig.geometry();
        auto geo_inside = cell_inside.geometry();
        auto geo_outside = cell_outside.geometry();

        // Get geometry of intersection in local coordinates of cell_inside and cell_outside
        auto geo_in_inside = ig.geometryInInside();
        auto geo_in_outside = ig.geometryInOutside();

        // Evaluate speed of sound (assumed constant per element)
        auto ref_el_inside = referenceElement(geo_inside);
        auto ref_el_outside = referenceElement(geo_outside);
        auto local_inside = ref_el_inside.position(0,0);
        auto local_outside = ref_el_outside.position(0,0);
        auto c_s = param.c(cell_inside,local_inside);
        auto c_n = param.c(cell_outside,local_outside);

        // Compute A+ (outgoing waves)
        Dune::FieldMatrix<DF,dim+1,dim+1> RT;
        LinearAcousticsEigenvectors<dim>::eigenvectors_transposed(c_s,n_F,RT);
        Dune::FieldVector<DF,dim+1> alpha;
        for (int i=0; i<=dim; i++) alpha[i] = RT[dim-1][i]; // row dim-1 corresponds to eigenvalue +c
        Dune::FieldVector<DF,dim+1> unit(0.0);
        unit[dim-1] = 1.0;
        Dune::FieldVector<DF,dim+1> beta;
        RT.solve(beta,unit);
        Dune::FieldMatrix<DF,dim+1,dim+1> A_plus_s;
        for (int i=0; i<=dim; i++)
          for (int j=0; j<=dim; j++)
            A_plus_s[i][j] = c_s*alpha[i]*beta[j];

        // Compute A- (incoming waves)
        LinearAcousticsEigenvectors<dim>::eigenvectors_transposed(c_n,n_F,RT);
        for (int i=0; i<=dim; i++) alpha[i] = RT[dim][i]; // row dim corresponds to eigenvalue -c
        unit = 0.0;
        unit[dim] = 1.0;
        RT.solve(beta,unit);
        Dune::FieldMatrix<DF,dim+1,dim+1> A_minus_n;
        for (int i=0; i<=dim; i++)
          for (int j=0; j<=dim; j++)
            A_minus_n[i][j] = -c_n*alpha[i]*beta[j];

        // Initialize vectors outside for loop
        Dune::FieldVector<RF,dim+1> u_s(0.0);
        Dune::FieldVector<RF,dim+1> u_n(0.0);
        Dune::FieldVector<RF,dim+1> f(0.0);

        // Loop over quadrature points
        const int order_s = dgspace_s.finiteElement().localBasis().order();
        const int order_n = dgspace_n.finiteElement().localBasis().order();
        const int intorder = overintegration+1+2*std::max(order_s,order_n);
        for (const auto& ip : quadratureRule(geo,intorder))
          {
            // Position of quadrature point in local coordinates of elements
            auto iplocal_s = geo_in_inside.global(ip.position());
            auto iplocal_n = geo_in_outside.global(ip.position());

            // Evaluate basis functions
            auto& phi_s = cache[order_s].evaluateFunction(iplocal_s,dgspace_s.finiteElement().localBasis());
            auto& phi_n = cache[order_n].evaluateFunction(iplocal_n,dgspace_n.finiteElement().localBasis());

            // Evaluate u from inside and outside
            u_s = 0.0;
            for (size_type i=0; i<=dim; i++) // for all components
              for (size_type k=0; k<dgspace_s.size(); k++) // for all basis functions
                u_s[i] += x_s(lfsv_s.child(i),k)*phi_s[k];
            u_n = 0.0;
            for (size_type i=0; i<=dim; i++) // for all components
              for (size_type k=0; k<dgspace_n.size(); k++) // for all basis functions
                u_n[i] += x_n(lfsv_n.child(i),k)*phi_n[k];

            // Compute numerical flux at integration point
            f = 0.0;
            A_plus_s.umv(u_s,f);
            // std::cout << "  after A_plus*u_s  " << f << std::endl;
            A_minus_n.umv(u_n,f);
            // std::cout << "  after A_minus*u_n " << f << std::endl;

            // Integrate
            auto factor = ip.weight() * geo.integrationElement(ip.position());
            for (size_type k=0; k<dgspace_s.size(); k++) // loop over all vector-valued (!) basis functions (with identical components)
              for (size_type i=0; i<=dim; i++) // loop over all components
                r_s.accumulate(lfsv_s.child(i),k, f[i]*phi_s[k]*factor);
            for (size_type k=0; k<dgspace_n.size(); k++) // loop over all vector-valued (!) basis functions (with identical components)
              for (size_type i=0; i<=dim; i++) // loop over all components
                r_n.accumulate(lfsv_n.child(i),k, - f[i]*phi_n[k]*factor);
          }

        // std::cout << "  residual_s: ";
        // for (size_type i=0; i<r_s.size(); i++) std::cout << r_s[i] << " ";
        // std::cout << std::endl;
        // std::cout << "  residual_n: ";
        // for (size_type i=0; i<r_n.size(); i++) std::cout << r_n[i] << " ";
        // std::cout << std::endl;

      }

      // Skeleton integral depending on test and ansatz functions
      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_boundary (const IG& ig,
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           R& r_s) const
      {
        // Get types
        using namespace Indices;
        using DGSpace = TypeTree::Child<LFSV,_0>;
        using DF = typename DGSpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType;
        using RF = typename DGSpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType;
        using size_type = typename DGSpace::Traits::SizeType;

        // Get local function space that is identical for all components
        const auto& dgspace_s = child(lfsv_s,_0);

        // Normal: assume faces are planar
        const Dune::FieldVector<DF,dim> n_F = ig.centerUnitOuterNormal();

        // Reference to inside cell
        const auto& cell_inside = ig.inside();

        // Get geometries
        auto geo = ig.geometry();
        auto geo_inside = cell_inside.geometry();

        // Get geometry of intersection in local coordinates of cell_inside
        auto geo_in_inside = ig.geometryInInside();

        // Evaluate speed of sound (assumed constant per element)
        auto ref_el_inside = referenceElement(geo_inside);
        auto local_inside = ref_el_inside.position(0,0);
        auto c_s = param.c(cell_inside,local_inside);

        // Compute A+ (outgoing waves)
        Dune::FieldMatrix<DF,dim+1,dim+1> RT;
        LinearAcousticsEigenvectors<dim>::eigenvectors_transposed(c_s,n_F,RT);
        Dune::FieldVector<DF,dim+1> alpha;
        for (int i=0; i<=dim; i++) alpha[i] = RT[dim-1][i]; // row dim-1 corresponds to eigenvalue +c
        Dune::FieldVector<DF,dim+1> unit(0.0);
        unit[dim-1] = 1.0;
        Dune::FieldVector<DF,dim+1> beta;
        RT.solve(beta,unit);
        Dune::FieldMatrix<DF,dim+1,dim+1> A_plus_s;
        for (int i=0; i<=dim; i++)
          for (int j=0; j<=dim; j++)
            A_plus_s[i][j] = c_s*alpha[i]*beta[j];

        // Compute A- (incoming waves)
        LinearAcousticsEigenvectors<dim>::eigenvectors_transposed(c_s,n_F,RT);
        for (int i=0; i<=dim; i++) alpha[i] = RT[dim][i]; // row dim corresponds to eigenvalue -c
        unit = 0.0;
        unit[dim] = 1.0;
        RT.solve(beta,unit);
        Dune::FieldMatrix<DF,dim+1,dim+1> A_minus_n;
        for (int i=0; i<=dim; i++)
          for (int j=0; j<=dim; j++)
            A_minus_n[i][j] = -c_s*alpha[i]*beta[j];

        // Initialize vectors outside for loop
        Dune::FieldVector<RF,dim+1> u_s(0.0);
        Dune::FieldVector<RF,dim+1> f(0.0);

        // Loop over quadrature points
        const int order_s = dgspace_s.finiteElement().localBasis().order();
        const int intorder = overintegration+1+2*order_s;
        for (const auto& ip : quadratureRule(geo,intorder))
          {
            // Position of quadrature point in local coordinates of elements
            auto iplocal_s = geo_in_inside.global(ip.position());

            // Evaluate basis functions
            auto& phi_s = cache[order_s].evaluateFunction(iplocal_s,dgspace_s.finiteElement().localBasis());

            // Evaluate u from inside and outside
            u_s = 0.0;
            for (size_type i=0; i<=dim; i++) // for all components
              for (size_type k=0; k<dgspace_s.size(); k++) // for all basis functions
                u_s[i] += x_s(lfsv_s.child(i),k)*phi_s[k];
            // std::cout << "  u_s " << u_s << std::endl;

            // Evaluate boundary condition
            Dune::FieldVector<RF,dim+1> u_n(param.g(ig.intersection(),ip.position(),u_s));
            // std::cout << "  u_n " << u_n << " bc: " << param.g(ig.intersection(),ip.position(),u_s) << std::endl;

            // Compute numerical flux at integration point
            f = 0.0;
            A_plus_s.umv(u_s,f);
            // std::cout << "  after A_plus*u_s  " << f << std::endl;
            A_minus_n.umv(u_n,f);
            // std::cout << "  after A_minus*u_n " << f << std::endl;

            // Integrate
            auto factor = ip.weight() * geo.integrationElement(ip.position());
            for (size_type k=0; k<dgspace_s.size(); k++) // loop over all vector-valued (!) basis functions (with identical components)
              for (size_type i=0; i<=dim; i++) // loop over all components
                r_s.accumulate(lfsv_s.child(i),k, f[i]*phi_s[k]*factor);
          }

        // std::cout << "  residual_s: ";
        // for (size_type i=0; i<r_s.size(); i++) std::cout << r_s[i] << " ";
        // std::cout << std::endl;
      }

      // Volume integral depending only on test functions
      template<typename EG, typename LFSV, typename R>
      void lambda_volume (const EG& eg, const LFSV& lfsv, R& r) const
      {
        // Get types
        using namespace Indices;
        using DGSpace = TypeTree::Child<LFSV,_0>;
        using size_type = typename DGSpace::Traits::SizeType;

        // Get local function space that is identical for all components
        const DGSpace& dgspace = child(lfsv,_0);

        // Reference to cell
        const auto& cell = eg.entity();

        // Get geometries
        auto geo = eg.geometry();

        // Loop over quadrature points
        const int order_s = dgspace.finiteElement().localBasis().order();
        const int intorder = overintegration+2*order_s;
        for (const auto& ip : quadratureRule(geo,intorder))
          {
            // Evaluate right hand side
            auto q(param.q(cell,ip.position()));

            // Evaluate basis functions
            auto& phi = cache[order_s].evaluateFunction(ip.position(),dgspace.finiteElement().localBasis());

            // Integrate
            auto factor = ip.weight() * geo.integrationElement(ip.position());
            for (size_type k=0; k<=dim; k++) // for all components
              for (size_type i=0; i<dgspace.size(); i++) // for all test functions of this component
                r.accumulate(lfsv.child(k),i, - q[k]*phi[i]*factor);
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
      using LocalBasisType = typename FEM::Traits::FiniteElementType::Traits::LocalBasisType;
      using Cache = Dune::PDELab::LocalBasisCache<LocalBasisType>;
      std::vector<Cache> cache;
    };



    /** a local operator for the mass operator of a vector valued lfs (L_2 integral)
     *
     * \f{align*}{
     \int_\Omega uv dx
     * \f}
     */
    template<typename T, typename FEM>
    class DGLinearAcousticsTemporalOperator :
      public NumericalJacobianApplyVolume<DGLinearAcousticsTemporalOperator<T,FEM> >,
        public LocalOperatorDefaultFlags,
        public InstationaryLocalOperatorDefaultMethods<typename T::Traits::RangeFieldType>
    {
      enum { dim = T::Traits::GridViewType::dimension };
    public:
      // pattern assembly flags
      enum { doPatternVolume = true };

      // residual assembly flags
      enum { doAlphaVolume = true };

      DGLinearAcousticsTemporalOperator (T& param_, int overintegration_=0)
        : param(param_), overintegration(overintegration_), cache(20)
      {}

      // define sparsity pattern of operator representation
      template<typename LFSU, typename LFSV, typename LocalPattern>
      void pattern_volume (const LFSU& lfsu, const LFSV& lfsv,
                           LocalPattern& pattern) const
      {
        // paranoia check number of number of components
        if (TypeTree::degree(lfsu)!=TypeTree::degree(lfsv))
          DUNE_THROW(Dune::Exception,"need U=V!");
        if (TypeTree::degree(lfsv)!=dim+1)
          DUNE_THROW(Dune::Exception,"need exactly dim+1 components!");

        for (size_t k=0; k<TypeTree::degree(lfsv); k++)
          for (size_t i=0; i<lfsv.child(k).size(); ++i)
            for (size_t j=0; j<lfsu.child(k).size(); ++j)
              pattern.addLink(lfsv.child(k),i,lfsu.child(k),j);
      }

      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
        // get types
        using namespace Indices;
        using DGSpace = TypeTree::Child<LFSV,_0>;
        using RF = typename DGSpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType;
        using size_type = typename DGSpace::Traits::SizeType;

        // get local function space that is identical for all components
        const auto& dgspace = child(lfsv,_0);

        // get geometry
        auto geo = eg.geometry();

        // Initialize vectors outside for loop
        Dune::FieldVector<RF,dim+1> u(0.0);

        // loop over quadrature points
        const int order = dgspace.finiteElement().localBasis().order();
        const int intorder = overintegration+2*order;
        for (const auto& ip : quadratureRule(geo,intorder))
          {
            // evaluate basis functions
            auto& phi = cache[order].evaluateFunction(ip.position(),dgspace.finiteElement().localBasis());

            // evaluate u
            u = 0.0;
            for (size_type k=0; k<=dim; k++) // for all components
              for (size_type j=0; j<dgspace.size(); j++) // for all basis functions
                u[k] += x(lfsv.child(k),j)*phi[j];

            // integrate
            auto factor = ip.weight() * geo.integrationElement(ip.position());
            for (size_type k=0; k<=dim; k++) // for all components
              for (size_type i=0; i<dgspace.size(); i++) // for all test functions of this component
                r.accumulate(lfsv.child(k),i, u[k]*phi[i]*factor);
          }
      }

      // jacobian of volume term
      template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                            M & mat) const
      {
        // get types
        using namespace Indices;
        using DGSpace = TypeTree::Child<LFSV,_0>;
        using size_type = typename DGSpace::Traits::SizeType;

        // get local function space that is identical for all components
        const auto& dgspace = child(lfsv,_0);

        // get geometry
        auto geo = eg.geometry();

        // loop over quadrature points
        const int order = dgspace.finiteElement().localBasis().order();
        const int intorder = overintegration+2*order;
        for (const auto& ip : quadratureRule(geo,intorder))
          {
            // evaluate basis functions
            auto& phi = cache[order].evaluateFunction(ip.position(),dgspace.finiteElement().localBasis());

            // integrate
            auto factor = ip.weight() * geo.integrationElement(ip.position());
            for (size_type k=0; k<=dim; k++) // for all components
              for (size_type i=0; i<dgspace.size(); i++) // for all test functions of this component
                for (size_type j=0; j<dgspace.size(); j++) // for all ansatz functions of this component
                  mat.accumulate(lfsv.child(k),i,lfsu.child(k),j, phi[j]*phi[i]*factor);
          }
      }

    private:
      T& param;
      int overintegration;
      using LocalBasisType = typename FEM::Traits::FiniteElementType::Traits::LocalBasisType;
      using Cache = Dune::PDELab::LocalBasisCache<LocalBasisType>;
      std::vector<Cache> cache;
    };

  }
}

#endif // DUNE_PDELAB_LOCALOPERATOR_LINEARACOUSTICSDG_HH
