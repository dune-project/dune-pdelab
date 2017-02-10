// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_LOCALOPERATOR_ERRORINDICATORDG_HH
#define DUNE_PDELAB_LOCALOPERATOR_ERRORINDICATORDG_HH

#include <algorithm>
#include <vector>

#include <dune/common/fvector.hh>

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/pdelab/localoperator/convectiondiffusiondg.hh>
#include <dune/pdelab/localoperator/convectiondiffusionparameter.hh>
#include <dune/pdelab/localoperator/eval.hh>
#include <dune/pdelab/localoperator/flags.hh>


// Note:
// The residual-based error estimator implemented here (for h-refinement only!)
// is taken from the PhD thesis
// "Robust A Posteriori Error Estimation for Discontinuous Galerkin Methods for Convection Diffusion Problems"
// by Liang Zhu (2010)
//

namespace Dune {
  namespace PDELab {

    /**
     * \brief a local operator for residual-based error estimation
     *
     * A call to residual() of a grid operator space will assemble
     * the quantity \f$\eta_T^2\f$ for each cell. Note that the squares
     * of the cell indicator \f$\eta_T\f$ is stored. To compute the global
     * error estimate sum up all values and take the square root.
     *
     * Assumptions and limitations:
     * - Assumes that LFSU is \f$P_1\f$/\f$Q_1\f$ finite element space
     *   and LFSV is a \f$P_0\f$ finite element space (one value per cell).
     * - Convection term is ignored (but reaction term is included)
     *
     * \tparam T model of ConvectionDiffusionParameterInterface
     */
    template<typename T>
    class ConvectionDiffusionDG_ErrorIndicator
      : public Dune::PDELab::LocalOperatorDefaultFlags
    {
      enum { dim = T::Traits::GridViewType::dimension };

      using Real = typename T::Traits::RangeFieldType;
      using BCType = typename ConvectionDiffusionBoundaryConditions::Type;

    public:
      // pattern assembly flags
      enum { doPatternVolume = false };
      enum { doPatternSkeleton = false };

      // residual assembly flags
      enum { doAlphaVolume  = true };
      enum { doAlphaSkeleton  = true };
      enum { doAlphaBoundary  = true };

      //! constructor: pass parameter object
      ConvectionDiffusionDG_ErrorIndicator ( T& param_,
                                             ConvectionDiffusionDGMethod::Type method_=ConvectionDiffusionDGMethod::NIPG,
                                             ConvectionDiffusionDGWeights::Type weights_=ConvectionDiffusionDGWeights::weightsOff,
                                             Real gamma_=0.0
                                             )
        : param(param_),
          method(method_),
          weights(weights_),
          gamma(gamma_) // interior penalty parameter, same as alpha in ConvectionDiffusionDG
      {}

      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
        // define types
        using RF = typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType;
        using RangeType = typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType;
        using size_type = typename LFSU::Traits::SizeType;

        // dimensions
        const int dim = EG::Geometry::mydimension;

        // Reference to cell
        const auto& cell = eg.entity();

        // Get geometry
        auto geo  = eg.geometry();

        // evaluate diffusion tensor at cell center, assume it is constant over elements
        auto ref_el = referenceElement(geo);
        auto localcenter = ref_el.position(0,0);
        auto A = param.A(cell,localcenter);

        static_assert(dim == 2 || dim == 3,
                      "The computation of epsilon looks very "
                      "much like it will only work in 2D or 3D.  If you think "
                      "otherwise, replace this static assert with a comment "
                      "that explains why.  --Jö");
        auto epsilon = std::min( A[0][0], A[1][1]);
        if( dim>2 ) epsilon = std::min( A[2][2], epsilon );

        // select quadrature rule
        // pOrder is constant on all grid elements (h-adaptive scheme).
        const int pOrder = lfsu.finiteElement().localBasis().order();

        // Initialize vectors outside for loop
        std::vector<RangeType> phi(lfsu.size());
        Dune::FieldVector<RF,dim> gradu(0.0);

        // loop over quadrature points
        RF sum(0.0);
        const int intorder = 2 * pOrder;
        for (const auto &qp : quadratureRule(geo,intorder))
          {
            // evaluate basis functions
            lfsu.finiteElement().localBasis().evaluateFunction(qp.position(),phi);

            // evaluate u
            RF u=0.0;
            for (size_type i=0; i<lfsu.size(); i++)
              u += x(lfsu,i)*phi[i];

            // evaluate reaction term
            auto c = param.c(cell,qp.position());

            // evaluate right hand side parameter function
            auto f = param.f(cell,qp.position());

            // evaluate convection term
            auto beta = param.b(cell,qp.position());


            /**********************/
            /* Evaluate Gradients */
            /**********************/
            gradu = 0.0;
            evalGradient( qp.position(), cell, lfsu, x, gradu );


            // integrate f^2
            auto factor = qp.weight() * geo.integrationElement(qp.position());

            auto square = f - (beta*gradu) - c*u; // + eps * Laplacian_u (TODO for pMax=2)
            square *= square;
            sum += square * factor;
          }

        // accumulate cell indicator
        auto h_T = diameter(geo);

        // L.Zhu: First term, interior residual squared
        auto eta_RK = h_T * h_T / pOrder / pOrder / epsilon * sum;

        // add contributions
        r.accumulate( lfsv, 0, eta_RK );
      }


      // skeleton integral depending on test and ansatz functions
      // each face is only visited ONCE!
      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_skeleton (const IG& ig,
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                           R& r_s, R& r_n) const
      {
        // Define types
        using RF = typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType;

        // Dimension
        const int dim = IG::dimension;

        // Refererences to inside and outside cells
        const auto& cell_inside = ig.inside();
        const auto& cell_outside = ig.outside();

        // Get geometries
        auto geo = ig.geometry();
        auto geo_inside = cell_inside.geometry();
        auto geo_outside = cell_outside.geometry();

        // Get geometry of intersection in local coordinates of inside_cell and outside_cell
        auto geo_in_inside = ig.geometryInInside();
        auto geo_in_outside = ig.geometryInOutside();


        // Evaluate permeability tensors
        auto ref_el_inside = referenceElement(geo_inside);
        auto ref_el_outside = referenceElement(geo_outside);
        auto local_inside = ref_el_inside.position(0,0);
        auto local_outside = ref_el_outside.position(0,0);
        auto A_s = param.A(cell_inside,local_inside);
        auto A_n = param.A(cell_outside,local_outside);

        static_assert(dim == 2 || dim == 3,
                      "The computation of epsilon_s and epsilon_n looks very "
                      "much like it will only work in 2D or 3D.  If you think "
                      "otherwise, replace this static assert with a comment "
                      "that explains why.  --Jö");

        auto epsilon_s = std::min( A_s[0][0], A_s[1][1]);
        if( dim>2 ) epsilon_s = std::min( A_s[2][2], epsilon_s );

        auto epsilon_n = std::min( A_n[0][0], A_n[1][1]);
        if( dim>2 ) epsilon_n = std::min( A_n[2][2], epsilon_n );

        const int pOrder_s = lfsu_s.finiteElement().localBasis().order();

        auto n_F = ig.centerUnitOuterNormal();

        RF flux_jump_L2normSquare(0.0);
        RF uh_jump_L2normSquare(0.0);

        // Declare vectors outside for loop
        Dune::FieldVector<RF,dim> An_F_s;
        Dune::FieldVector<RF,dim> An_F_n;
        Dune::FieldVector<RF,dim> gradu_s;
        Dune::FieldVector<RF,dim> gradu_n;

        // loop over quadrature points and integrate normal flux
        const int intorder = 2*pOrder_s;
        for (const auto &qp : quadratureRule(geo,intorder))
          {
            // position of quadrature point in local coordinates of elements
            const auto &iplocal_s = geo_in_inside .global(qp.position());
            const auto &iplocal_n = geo_in_outside.global(qp.position());

            // Diffusion tensor at quadrature point
            A_s = param.A( cell_inside,  iplocal_s );
            A_n = param.A( cell_outside, iplocal_n );

            A_s.mv(n_F,An_F_s);
            A_n.mv(n_F,An_F_n);

            /**********************/
            /* Evaluate Functions */
            /**********************/

            // evaluate uDG_s, uDG_n
            RF uDG_s=0.0;
            evalFunction( iplocal_s, lfsu_s, x_s, uDG_s );

            RF uDG_n=0.0;
            evalFunction( iplocal_n, lfsu_n, x_n, uDG_n );


            /**********************/
            /* Evaluate Gradients */
            /**********************/
            gradu_s = 0.0;
            evalGradient( iplocal_s, cell_inside, lfsu_s, x_s, gradu_s );
            gradu_n = 0.0;
            evalGradient( iplocal_n, cell_outside, lfsu_n, x_n, gradu_n );


            // integrate
            auto factor = qp.weight() * geo.integrationElement(qp.position());

            // evaluate flux jump term
            auto flux_jump = (An_F_s*gradu_s)-(An_F_n*gradu_n);
            flux_jump_L2normSquare += flux_jump * flux_jump * factor;

            // evaluate jump term
            auto jump_uDG = (uDG_s - uDG_n);
            uh_jump_L2normSquare += jump_uDG * jump_uDG * factor ;
          }

        // accumulate indicator
        auto h_face = diameter(geo);

        // L.Zhu: second term, edge residual
        auto eta_Ek_s = 0.5 * h_face * flux_jump_L2normSquare;
        auto eta_Ek_n = eta_Ek_s; // equal on both sides of the intersection!

        // L.Zhu: third term, edge jumps
        auto eta_Jk_s = 0.5 * ( gamma / h_face + h_face / epsilon_s) * uh_jump_L2normSquare;
        auto eta_Jk_n = 0.5 * ( gamma / h_face + h_face / epsilon_n) * uh_jump_L2normSquare;

        // add contributions from both sides of the intersection
        r_s.accumulate( lfsv_s, 0, eta_Ek_s + eta_Jk_s );
        r_n.accumulate( lfsv_n, 0, eta_Ek_n + eta_Jk_n );
      }


      // boundary integral depending on test and ansatz functions
      // We put the Dirchlet evaluation also in the alpha term to save some geometry evaluations
      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_boundary (const IG& ig,
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           R& r_s) const
      {
        // define types
        using RF = typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType;

        // dimensions
        const int dim = IG::dimension;

        // Reference to inside cell
        const auto& cell_inside = ig.inside();

        // Get geometries
        auto geo = ig.geometry();
        auto geo_inside = cell_inside.geometry();

        // Get geometry of intersection in local coordinates of inside_cell
        auto geo_in_inside = ig.geometryInInside();

        // reference elements
        auto ref_el = referenceElement(geo);
        auto ref_el_inside = referenceElement(geo_inside);

        // evaluate permeability tensors
        auto local_inside = ref_el_inside.position(0,0);
        auto A_s = param.A(cell_inside,local_inside);

        static_assert(dim == 2 || dim == 3,
                      "The computation of epsilon_s looks very "
                      "much like it will only work in 2D or 3D.  If you think "
                      "otherwise, replace this static assert with a comment "
                      "that explains why.  --Jö");

        auto epsilon_s = std::min( A_s[0][0], A_s[1][1]);
        if( dim>2 ) epsilon_s = std::min( A_s[2][2], epsilon_s );

        // select quadrature rule
        const int pOrder_s = lfsu_s.finiteElement().localBasis().order();
        const int intorder = 2*pOrder_s;

        // evaluate boundary condition
        auto face_local = ref_el.position(0,0);
        auto bctype = param.bctype(ig.intersection(),face_local);
        if (bctype != ConvectionDiffusionBoundaryConditions::Dirichlet)
          return;

        RF uh_jump_L2normSquare(0.0);

        // loop over quadrature points and integrate normal flux
        for (const auto &qp : quadratureRule(geo,intorder))
          {
            // position of quadrature point in local coordinates of elements
            const auto &iplocal_s = geo_in_inside.global(qp.position());

            // evaluate Dirichlet boundary condition
            auto gDirichlet = param.g( cell_inside, iplocal_s );

            /**********************/
            /* Evaluate Functions */
            /**********************/

            // evaluate uDG_s
            RF uDG_s=0.0;
            evalFunction( iplocal_s, lfsu_s, x_s, uDG_s );

            // integrate
            auto factor = qp.weight() * geo.integrationElement(qp.position());

            // evaluate jump term
            auto jump_uDG = uDG_s - gDirichlet;
            uh_jump_L2normSquare += jump_uDG * jump_uDG * factor;
          }

        // accumulate indicator
        auto h_face = diameter(geo);

        // L.Zhu: third term, edge jumps on the Dirichlet boundary
        auto eta_Jk_s = (gamma / h_face + h_face / epsilon_s) * uh_jump_L2normSquare;

        // add contributions
        r_s.accumulate( lfsv_s, 0, eta_Jk_s );  // boundary edge
      }

    private:
      T& param;  // two phase parameter class
      ConvectionDiffusionDGMethod::Type method;
      ConvectionDiffusionDGWeights::Type weights;
      Real gamma; // interior penalty parameter, same as alpha in class TransportOperatorDG

      template<class GEO>
      typename GEO::ctype diameter (const GEO& geo) const
      {
        using DF = typename GEO::ctype;
        DF hmax = -1.0E00;
        const int dim = GEO::coorddimension;
        for (int i=0; i<geo.corners(); i++)
          {
            Dune::FieldVector<DF,dim> xi = geo.corner(i);
            for (int j=i+1; j<geo.corners(); j++)
              {
                Dune::FieldVector<DF,dim> xj = geo.corner(j);
                xj -= xi;
                hmax = std::max(hmax,xj.two_norm());
              }
          }
        return hmax;
      }

    };


  }
}
#endif // DUNE_PDELAB_LOCALOPERATOR_ERRORINDICATORDG_HH
