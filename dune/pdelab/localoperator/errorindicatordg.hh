// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_ERRORINDICATORDG_HH
#define DUNE_PDELAB_ERRORINDICATORDG_HH

#include<vector>

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/geometry/type.hh>
#include<dune/geometry/referenceelements.hh>
#include<dune/geometry/quadraturerules.hh>

#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/idefault.hh>
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/finiteelement/localbasiscache.hh>

#include "convectiondiffusionparameter.hh"
#include "convectiondiffusiondg.hh"
#include "eval.hh"

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

      typedef typename T::Traits::RangeFieldType Real;
      typedef typename ConvectionDiffusionBoundaryConditions::Type BCType;

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
        // domain and range field type
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RangeType;
        typedef typename LFSU::Traits::SizeType size_type;

        // dimensions
        const int dim = EG::Geometry::dimension;

        // extract objects
        auto cell = eg.entity();
        auto geo  = eg.geometry();

        const auto &gt = geo.type();

        const auto &ref = Dune::ReferenceElements<DF,dim>::general(gt);

        const auto &localcenter = ref.position(0,0);

        // Diffusion tensor at cell center
        auto A = param.A(cell,localcenter);

        static_assert(dim == 2 || dim == 3,
                      "The computation of epsilon looks very "
                      "much like it will only work in 2D or 3D.  If you think "
                      "otherwise, replace this static assert with a comment "
                      "that explains why.  --Jö");
        RF epsilon = std::min( A[0][0], A[1][1]);
        if( dim>2 ) epsilon = std::min( A[2][2], epsilon );

        // select quadrature rule
        // pOrder is constant on all grid elements (h-adaptive scheme).
        const int pOrder = lfsu.finiteElement().localBasis().order();
        const int intorder = 2 * pOrder;
        const auto rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

        RF sum(0.0);
        std::vector<RangeType> phi(lfsu.size());

        // loop over quadrature points
        for (const auto &qp : rule)
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

            Dune::FieldVector<RF,dim> gradu(0.0);
            evalGradient( qp.position(), cell, lfsu, x, gradu );


            // integrate f^2
            RF factor = qp.weight() * geo.integrationElement(qp.position());

            RF square = f - (beta*gradu) - c*u; // + eps * Laplacian_u (TODO for pMax=2)
            square *= square;
            sum += square * factor;
          }

        // accumulate cell indicator
        DF h_T = diameter(geo);

        // L.Zhu: First term, interior residual squared
        RF eta_RK = h_T * h_T / pOrder / pOrder / epsilon * sum;

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
        // domain and range field type
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;

        // dimensions
        const int dim = IG::dimension;

        // extract objects
        auto cell_inside = ig.inside();
        auto cell_outside = ig.outside();

        auto geo            = ig.geometry();
        auto geo_in_inside  = ig.geometryInInside();
        auto geo_in_outside = ig.geometryInOutside();

        const auto &gtface = geo.type();

        const auto &insideRef =
          Dune::ReferenceElements<DF,dim>::general(cell_inside.type());
        const auto &outsideRef =
          Dune::ReferenceElements<DF,dim>::general(cell_outside.type());

        const auto &inside_local = insideRef.position(0,0);
        const auto &outside_local = outsideRef.position(0,0);

        // evaluate permeability tensors
        auto A_s = param.A(cell_inside,inside_local);
        auto A_n = param.A(cell_outside,outside_local);

        static_assert(dim == 2 || dim == 3,
                      "The computation of epsilon_s and epsilon_n looks very "
                      "much like it will only work in 2D or 3D.  If you think "
                      "otherwise, replace this static assert with a comment "
                      "that explains why.  --Jö");
        RF epsilon_s = std::min( A_s[0][0], A_s[1][1]);
        if( dim>2 ) epsilon_s = std::min( A_s[2][2], epsilon_s );

        RF epsilon_n = std::min( A_n[0][0], A_n[1][1]);
        if( dim>2 ) epsilon_n = std::min( A_n[2][2], epsilon_n );

        // select quadrature rule
        const int pOrder_s = lfsu_s.finiteElement().localBasis().order();
        const int intorder = 2*pOrder_s;
        const auto& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,intorder);

        const auto &n_F = ig.centerUnitOuterNormal();

        RF flux_jump_L2normSquare(0.0);
        RF uh_jump_L2normSquare(0.0);

        // loop over quadrature points and integrate normal flux
        for (const auto &qp : rule)
          {
            // position of quadrature point in local coordinates of elements
            const auto &iplocal_s = geo_in_inside .global(qp.position());
            const auto &iplocal_n = geo_in_outside.global(qp.position());

            // Diffusion tensor at quadrature point
            A_s = param.A( cell_inside,  iplocal_s );
            A_n = param.A( cell_outside, iplocal_n );

            Dune::FieldVector<RF,dim> An_F_s;
            A_s.mv(n_F,An_F_s);
            Dune::FieldVector<RF,dim> An_F_n;
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

            Dune::FieldVector<RF,dim> gradu_s(0.0);
            evalGradient( iplocal_s, cell_inside, lfsu_s, x_s, gradu_s );

            Dune::FieldVector<RF,dim> gradu_n(0.0);
            evalGradient( iplocal_n, cell_outside, lfsu_n, x_n, gradu_n );


            // integrate
            RF factor = qp.weight() * geo.integrationElement(qp.position());

            // evaluate flux jump term
            RF flux_jump = (An_F_s*gradu_s)-(An_F_n*gradu_n);
            flux_jump_L2normSquare += flux_jump * flux_jump * factor;

            // evaluate jump term
            RF jump_uDG = (uDG_s - uDG_n);
            uh_jump_L2normSquare += jump_uDG * jump_uDG * factor ;
          }

        // accumulate indicator
        DF h_face = diameter(ig.geometry());

        // L.Zhu: second term, edge residual
        RF eta_Ek_s = 0.5 * h_face * flux_jump_L2normSquare;
        RF eta_Ek_n = eta_Ek_s; // equal on both sides of the intersection!

        // L.Zhu: third term, edge jumps
        RF eta_Jk_s = 0.5 * ( gamma / h_face + h_face / epsilon_s) * uh_jump_L2normSquare;
        RF eta_Jk_n = 0.5 * ( gamma / h_face + h_face / epsilon_n) * uh_jump_L2normSquare;

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
        // domain and range field type
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;

        // dimensions
        const int dim = IG::dimension;

        // extract objects
        auto cell_inside = ig.inside();

        auto geo            = ig.geometry();
        auto geo_in_inside  = ig.geometryInInside();

        const auto &gtface = geo.type();

        const auto &ref =
          Dune::ReferenceElements<DF,dim-1>::general(gtface);
        const auto &insideRef =
          Dune::ReferenceElements<DF,dim>::general(cell_inside.type());

        const auto &inside_local = insideRef.position(0,0);

        // evaluate permeability tensors
        auto A_s = param.A(cell_inside,inside_local);

        static_assert(dim == 2 || dim == 3,
                      "The computation of epsilon_s looks very "
                      "much like it will only work in 2D or 3D.  If you think "
                      "otherwise, replace this static assert with a comment "
                      "that explains why.  --Jö");
        RF epsilon_s = std::min( A_s[0][0], A_s[1][1]);
        if( dim>2 ) epsilon_s = std::min( A_s[2][2], epsilon_s );

        // select quadrature rule
        const int pOrder_s = lfsu_s.finiteElement().localBasis().order();
        const int intorder = 2*pOrder_s;
        const auto& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,intorder);

        // evaluate boundary condition
        const auto &face_local = ref.position(0,0);
        BCType bctype = param.bctype(ig.intersection(),face_local);
        if (bctype != ConvectionDiffusionBoundaryConditions::Dirichlet)
          return;

        RF uh_jump_L2normSquare(0.0);

        // loop over quadrature points and integrate normal flux
        for (const auto &qp : rule)
          {
            // position of quadrature point in local coordinates of elements
            const auto &iplocal_s = geo_in_inside .global(qp.position());

            // evaluate Dirichlet boundary condition
            RF gDirichlet = param.g( cell_inside, iplocal_s );

            /**********************/
            /* Evaluate Functions */
            /**********************/

            // evaluate uDG_s
            RF uDG_s=0.0;
            evalFunction( iplocal_s, lfsu_s, x_s, uDG_s );

            // integrate
            RF factor = qp.weight() * geo.integrationElement(qp.position());

            // evaluate jump term
            RF jump_uDG = uDG_s - gDirichlet;
            uh_jump_L2normSquare += jump_uDG * jump_uDG * factor;
          }

        // accumulate indicator
        DF h_face = diameter(ig.geometry());

        // L.Zhu: third term, edge jumps on the Dirichlet boundary
        RF eta_Jk_s = (gamma / h_face + h_face / epsilon_s) * uh_jump_L2normSquare;

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
        typedef typename GEO::ctype DF;
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
#endif
