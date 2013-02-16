// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_ERRORINDICATORDG_HH
#define DUNE_PDELAB_ERRORINDICATORDG_HH

#include<vector>

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/static_assert.hh>
#include<dune/geometry/type.hh>
#include<dune/geometry/referenceelements.hh>
#include<dune/geometry/quadraturerules.hh>

#include<dune/pdelab/common/geometrywrapper.hh>
#include<dune/pdelab/gridoperatorspace/gridoperatorspace.hh>
#include<dune/pdelab/gridoperatorspace/gridoperatorspaceutilities.hh>
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

    /** a local operator for residual-based error estimation
     *  
     * A call to residual() of a grid operator space will assemble
     * the quantity \eta_T^2 for each cell. Note that the squares
     * of the cell indicator \eta_T is stored. To compute the global
     * error estimate sum up all values and take the square root.
     *
     * Assumptions and limitations:
     * - Assumes that LFSU is P_1/Q_1 finite element space
     *   and LFSV is a P_0 finite element space (one value per cell).
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

        // pOrder is constant on all grid elements (h-adaptive scheme).
        const int pOrder = lfsu.finiteElement().localBasis().order();
        const int intorder = 2 * pOrder;
        
        Dune::GeometryType gt = eg.geometry().type();
        // Diffusion tensor at cell center
        typename T::Traits::PermTensorType A;
        Dune::FieldVector<DF,dim> localcenter = Dune::GenericReferenceElements<DF,dim>::general(gt).position(0,0);
        A = param.A(eg.entity(),localcenter);
        RF epsilon = std::min( A[0][0], A[1][1]);
        if( dim>2 ) epsilon = std::min( A[2][2], epsilon );

        // select quadrature rule
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

        // loop over quadrature points
        RF sum(0.0);
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // evaluate basis functions
            std::vector<RangeType> phi(lfsu.size());
            lfsu.finiteElement().localBasis().evaluateFunction(it->position(),phi);

            // evaluate u
            RF u=0.0;
            for (size_type i=0; i<lfsu.size(); i++)
              u += x(lfsu,i)*phi[i];

            // evaluate reaction term
            typename T::Traits::RangeFieldType c = param.c(eg.entity(),it->position());

            // evaluate right hand side parameter function
            typename T::Traits::RangeFieldType f = param.f(eg.entity(),it->position());

            // integrate f^2
            RF factor = it->weight() * eg.geometry().integrationElement(it->position());
      
            // evaluate convection term
            typename T::Traits::RangeType beta
              = param.b(eg.entity(),it->position());

            // compute gradient of u_h
            Dune::FieldVector<RF,dim> gradu(0.0);
            evalGradient( it->position(), eg.entity(), lfsu, x, gradu );

            RF square = f - (beta*gradu) - c*u; // + eps * Laplacian_u (TODO for pMax=2)
            square *= square;
            sum += square * factor;
          }
        
        // accumulate cell indicator 
        DF h_T = diameter(eg.geometry());

        // L.Zhu: First term, interior residual squared
        RF eta_RK = h_T * h_T / pOrder / pOrder / epsilon * sum;
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
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RangeType;
        
        // dimensions
        const int dim = IG::dimension;
        
        // evaluate permeability tensors
        const Dune::FieldVector<DF,dim>&
          inside_local = Dune::GenericReferenceElements<DF,dim>::general(ig.inside()->type()).position(0,0);
        const Dune::FieldVector<DF,dim>&
          outside_local = Dune::GenericReferenceElements<DF,dim>::general(ig.outside()->type()).position(0,0);
        typename T::Traits::PermTensorType A_s, A_n;
        A_s = param.A(*(ig.inside()),inside_local);
        A_n = param.A(*(ig.outside()),outside_local);

        RF epsilon_s = std::min( A_s[0][0], A_s[1][1]);
        if( dim>2 ) epsilon_s = std::min( A_s[2][2], epsilon_s );

        RF epsilon_n = std::min( A_n[0][0], A_n[1][1]);
        if( dim>2 ) epsilon_n = std::min( A_n[2][2], epsilon_n );

        // select quadrature rule
        const int pOrder_s = lfsu_s.finiteElement().localBasis().order();
        const int intorder = 2*pOrder_s;
        Dune::GeometryType gtface = ig.geometryInInside().type();
        const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,intorder);

        const Dune::FieldVector<DF,dim> n_F = ig.centerUnitOuterNormal();

        RF flux_jump_L2normSquare(0.0);
        RF uh_jump_L2normSquare(0.0);

        // loop over quadrature points and integrate normal flux
        for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // position of quadrature point in local coordinates of elements 
            Dune::FieldVector<DF,dim> iplocal_s = ig.geometryInInside().global(it->position());
            Dune::FieldVector<DF,dim> iplocal_n = ig.geometryInOutside().global(it->position());

            // Diffusion tensor at quadrature point
            typename T::Traits::PermTensorType A_s = param.A( *(ig.inside()), iplocal_s );
            typename T::Traits::PermTensorType A_n = param.A( *(ig.outside()), iplocal_n);

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
            evalGradient( iplocal_s, *(ig.inside()), lfsu_s, x_s, gradu_s );

            Dune::FieldVector<RF,dim> gradu_n(0.0);
            evalGradient( iplocal_n, *(ig.outside()), lfsu_n, x_n, gradu_n );


            // integrate
            RF factor = it->weight() * ig.geometry().integrationElement(it->position());

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
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RangeType;
        
        // dimensions
        const int dim = IG::dimension;
        
        // evaluate permeability tensors
        const Dune::FieldVector<DF,dim>&
          inside_local = Dune::GenericReferenceElements<DF,dim>::general(ig.inside()->type()).position(0,0);
        typename T::Traits::PermTensorType A_s;
        A_s = param.A(*(ig.inside()),inside_local);

        RF epsilon_s = std::min( A_s[0][0], A_s[1][1]);
        if( dim>2 ) epsilon_s = std::min( A_s[2][2], epsilon_s );
        
        const Dune::FieldVector<DF,dim> n_F = ig.centerUnitOuterNormal();

        // select quadrature rule
        const int pOrder_s = lfsu_s.finiteElement().localBasis().order();
        const int intorder = 2 * pOrder_s;

        Dune::GeometryType gtface = ig.geometryInInside().type();
        const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,intorder);

        // evaluate boundary condition
        const Dune::FieldVector<DF,dim-1>
          face_local = Dune::GenericReferenceElements<DF,dim-1>::general(gtface).position(0,0);
        BCType bctype = param.bctype(ig.intersection(),face_local);
        if (bctype != ConvectionDiffusionBoundaryConditions::Dirichlet)
          return;

        // loop over quadrature points and integrate normal flux
        RF uh_jump_L2normSquare(0.0);
        for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // position of quadrature point in local coordinates of elements 
            Dune::FieldVector<DF,dim> iplocal_s = ig.geometryInInside().global(it->position());

            // evaluate Dirichlet boundary condition
            RF gDirichlet = param.g( *(ig.inside()), iplocal_s );

            /**********************/
            /* Evaluate Functions */
            /**********************/

            // evaluate uDG_s
            RF uDG_s=0.0;
            evalFunction( iplocal_s, lfsu_s, x_s, uDG_s );
            RF jump_uDG = uDG_s - gDirichlet;
                
            // integrate
            RF factor = it->weight() * ig.geometry().integrationElement(it->position());

            uh_jump_L2normSquare += jump_uDG * jump_uDG * factor;

          }

        // accumulate indicator
        DF h_face = diameter(ig.geometry());
        
        // L.Zhu: third term, edge jumps on the Dirichlet boundary
        RF eta_Jk_s = (gamma / h_face + h_face / epsilon_s) * uh_jump_L2normSquare;
        
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
