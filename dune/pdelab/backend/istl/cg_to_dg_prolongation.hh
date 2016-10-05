// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_BACKEND_ISTL_CG_TO_DG_PROLONGATION_HH
#define DUNE_PDELAB_BACKEND_ISTL_CG_TO_DG_PROLONGATION_HH

// this is here for backwards compatibility and deprecation warnings, remove after 2.5.0
#include "ensureistlinclude.hh"

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

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

    namespace CG2DGHelper { // hide some TMP code

      template <typename Imp>
      struct WrappedLocalShapeFunctionTraits {
        typedef typename Imp::Traits::FiniteElementType::
        Traits::LocalBasisType::Traits::RangeType RangeType;
        typedef typename Imp::Traits::FiniteElementType::
        Traits::LocalBasisType::Traits::DomainType DomainType;
      };

      // evaluate a localfunction as a function on a different element
      template<typename Imp>
      class WrappedLocalShapeFunction
      {
        const Imp & _imp;
        const int _comp;

        typedef typename Imp::Traits::FiniteElementType FEM;
        typedef FiniteElementInterfaceSwitch<FEM> FESwitch;
        typedef BasisInterfaceSwitch<typename FESwitch::Basis > BasisSwitch;
        typedef typename BasisSwitch::DomainField DF;
        typedef typename BasisSwitch::Range RT;
        enum { dim = BasisSwitch::dimDomainLocal };
      public:
        typedef WrappedLocalShapeFunction<Imp> Traits;
        WrappedLocalShapeFunction (const Imp& imp, int comp) :
          _imp(imp), _comp(comp) {}

        void evaluate(const Dune::FieldVector<DF,dim> & x,
          Dune::FieldVector<DF,1> & y) const
        {
          std::vector<RT> v;
          _imp.finiteElement().localBasis().evaluateFunction(x,v);
          y = v[_comp];
        }
      };

      template <typename R>
      class ComputeCG2DGVisitor :
        public TypeTree::DefaultPairVisitor,
        public TypeTree::DynamicTraversal,
        public TypeTree::VisitTree
      {
        LocalMatrix<R>& _mat;

      public:
        ComputeCG2DGVisitor(LocalMatrix<R>& mat) :
          _mat(mat)
        {}

        template<typename LFSU, typename LFSV, typename TreePath>
        void leaf(const LFSU& lfsu, const LFSV& lfsv, TreePath treePath) const
        {
          // map from CG (lfsu) 2 DG (lfsv)
          typedef typename LFSV::Traits::FiniteElementType DG_FEM;
          typedef FiniteElementInterfaceSwitch<DG_FEM> FESwitch;
          typedef BasisInterfaceSwitch<typename FESwitch::Basis > BasisSwitch;
          typedef typename BasisSwitch::DomainField DF;
          std::vector<DF> v;
          for (unsigned int i=0; i<lfsu.size(); i++)
          {
            // create function f, which wraps a CG shape function
            WrappedLocalShapeFunction<LFSU> f(lfsu, i);
            // interpolate f into DG space
            FESwitch::interpolation(lfsv.finiteElement()).
              interpolate(f, v);
            // store coefficients
            for (unsigned int j=0; j<lfsv.size(); j++)
            {
              _mat(lfsv,j,lfsu,i) = v[j];
            }
          }
        }
      };

    } // end namespace CG2DGHelper

    // a local operator to compute DG shift matrix needed for some AMG variants
    class CG2DGProlongation :
      public FullVolumePattern,
      public LocalOperatorDefaultFlags,
      public InstationaryLocalOperatorDefaultMethods<double>
    {
      template<typename LFSU, typename LFSV, typename R>
      void computeCG2DG(const LFSU & lfsu, const LFSV & lfsv,
        LocalMatrix<R>& mat) const
      {
        // lfsu: CG
        // lfsv: DG
        CG2DGHelper::ComputeCG2DGVisitor<R> cg2dg(mat);
        TypeTree::applyToTreePair(lfsu, lfsv, cg2dg);
      }
    public:
      // pattern assembly flags
      enum { doPatternVolume = true };

      // residual assembly flags
      enum { doAlphaVolume  = true };

      CG2DGProlongation  () {}

      // alpha_volume:
      // not implemented, as it should never be used. We just miss-use the assembler to
      // assemble the shift-matrix

      // jacobian of skeleton term
      template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_volume (const EG&, const LFSU& lfsu, const X&, const LFSV& lfsv,
                            M & mat) const
      {
        computeCG2DG(lfsu, lfsv, mat.container());
      }
    };

    //! \} group LocalOperator
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_BACKEND_ISTL_CG_TO_DG_PROLONGATION_HH
