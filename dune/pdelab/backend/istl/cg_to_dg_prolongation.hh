// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_BACKEND_ISTL_CG_TO_DG_PROLONGATION_HH
#define DUNE_PDELAB_BACKEND_ISTL_CG_TO_DG_PROLONGATION_HH

// this is here for backwards compatibility and deprecation warnings, remove after 2.5.0
#include "ensureistlinclude.hh"

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/dynmatrix.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/localfunctions/common/interfaceswitch.hh>

#include <dune/typetree/pairtraversal.hh>
#include <dune/typetree/transformation.hh>
#include <dune/typetree/visitor.hh>

#include <dune/pdelab/common/geometrywrapper.hh>
#include <dune/pdelab/gridoperator/common/localmatrix.hh>
#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/localoperator/idefault.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup LocalOperator
    //! \ingroup PDELab
    //! \{

    // a local operator to compute DG shift matrix needed for some AMG variants
    class CG2DGProlongation :
      public FullVolumePattern,
      public LocalOperatorDefaultFlags,
      public InstationaryLocalOperatorDefaultMethods<double>
    {
    public:
      // pattern assembly flags
      enum { doPatternVolume = true };

      // residual assembly flags
      enum { doAlphaVolume  = true };

      CG2DGProlongation() = default;

      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
        // Switches between local and global interface
        typedef FiniteElementInterfaceSwitch<
          typename LFSU::Traits::FiniteElementType
          > FESwitchU;
        typedef FiniteElementInterfaceSwitch<
          typename LFSV::Traits::FiniteElementType
          > FESwitchV;
        typedef BasisInterfaceSwitch<
          typename FESwitchU::Basis
          > BasisSwitch;

        // domain and range field type
        typedef typename BasisSwitch::DomainField DF;
        typedef typename BasisSwitch::RangeField RF;
        typedef typename BasisSwitch::Range RangeType;
        typedef typename LFSU::Traits::SizeType size_type;

        // dimensions
        const int dim = EG::Geometry::mydimension;

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        int intorder = 2+2*std::max(lfsu.finiteElement().localBasis().order(),lfsv.finiteElement().localBasis().order());

        // matrices
        //std::cout << "n_DG=" << n_DG << " n_CG=" << n_CG << std::endl;
        size_type n_DG = lfsv.size();
        size_type n_CG = lfsu.size();
        Dune::DynamicMatrix<RF> MDG(n_DG,n_DG,0.0); // DG mass matrix in element
        Dune::DynamicVector<RF> b(n_DG,0.0); // right hand side vector

        // loop over quadrature points
        std::vector<RangeType> phi(n_DG); // values of DG basis functions
        std::vector<RangeType> psi(n_CG); // values of CG basis functions

        for(const auto& ip : Dune::QuadratureRules<DF,dim>::rule(gt,intorder)) {
          // evaluate basis functions
          FESwitchU::basis(lfsu.finiteElement()).evaluateFunction(ip.position(),psi); // eval CG basis ("Ansatz functions")
          FESwitchV::basis(lfsv.finiteElement()).evaluateFunction(ip.position(),phi); // eval DG basis ("Test functions")

          RF factor = ip.weight() * eg.geometry().integrationElement(ip.position());

          // accumulate to local mass matrix
          for (size_type i=0; i<n_DG; i++)
            for (size_type j=0; j<n_DG; j++)
              MDG[i][j] += phi[j]*phi[i]*factor;

          RF u = 0.0;
          for(size_type i=0; i<n_CG; i++)
            u += x(lfsu,i) * psi[i];

          // accumulate to RHS vector
          for(size_type i=0; i<n_DG; i++)
            b[i] += u * phi[i] * factor;
        }

        // invert mass matrix
        MDG.invert();

        // compute MDG^-1 b which is the prolongation
        for(size_type i=0; i<n_DG; i++) {
          RF res_entry = 0.0;
          for(size_type j=0; j<n_DG; j++)
            res_entry += MDG[i][j] * b[j];
          r.accumulate(lfsv,i, res_entry);
        }
      }

      // jacobian of volume term
      template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_volume (const EG& eg, const LFSU& lfsu, const X&, const LFSV& lfsv,
                            M & mat) const
      {
        // Switches between local and global interface
        typedef FiniteElementInterfaceSwitch<
          typename LFSU::Traits::FiniteElementType
          > FESwitchU;
        typedef FiniteElementInterfaceSwitch<
          typename LFSV::Traits::FiniteElementType
          > FESwitchV;
        typedef BasisInterfaceSwitch<
          typename FESwitchU::Basis
          > BasisSwitch;

        // domain and range field type
        typedef typename BasisSwitch::DomainField DF;
        typedef typename BasisSwitch::RangeField RF;
        typedef typename BasisSwitch::Range RangeType;
        typedef typename LFSU::Traits::SizeType size_type;

        // dimensions
        const int dim = EG::Geometry::mydimension;

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        int intorder = 2+2*std::max(lfsu.finiteElement().localBasis().order(),lfsv.finiteElement().localBasis().order());

        // matrices
        //std::cout << "n_DG=" << n_DG << " n_CG=" << n_CG << std::endl;
        size_type n_DG = lfsv.size();
        size_type n_CG = lfsu.size();
        typedef Dune::DynamicMatrix<RF> MATRIX;
        MATRIX MDG(n_DG,n_DG,0.0); // DG mass matrix in element
        MATRIX B(n_DG,n_CG,0.0); // right hand side matrix

        // loop over quadrature points
        std::vector<RangeType> phi(n_DG); // values of DG basis functions
        std::vector<RangeType> psi(n_CG); // values of CG basis functions

        for(const auto& ip : Dune::QuadratureRules<DF,dim>::rule(gt,intorder)) {
          // evaluate basis functions
          FESwitchU::basis(lfsu.finiteElement()).evaluateFunction(ip.position(),psi); // eval CG basis ("Ansatz functions")
          FESwitchV::basis(lfsv.finiteElement()).evaluateFunction(ip.position(),phi); // eval DG basis ("Test functions")

          RF factor = ip.weight() * eg.geometry().integrationElement(ip.position());

          // accumulate to local mass matrix
          for (size_type i=0; i<n_DG; i++)
            for (size_type j=0; j<n_DG; j++)
              MDG[i][j] += phi[j]*phi[i]*factor;

          // accumulate to rhs matrix
          for (size_type i=0; i<n_DG; i++)
            for (size_type j=0; j<n_CG; j++)
              B[i][j] += psi[j]*phi[i]*factor;
        }

        // invert mass matrix
        MDG.invert();

        // compute MDG^-1 B which is the prolongation matrix
        for(size_type i=0; i<n_DG; i++)
          for(size_type j=0; j<n_CG; j++) {
            RF mat_entry = 0.0;
            for(size_type k=0; k<n_DG; k++)
              mat_entry += MDG[i][k] * B[k][j];
            mat.accumulate(lfsv,i,lfsu,j, mat_entry);
          }
      }
    };

    //! \} group LocalOperator
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_BACKEND_ISTL_CG_TO_DG_PROLONGATION_HH
