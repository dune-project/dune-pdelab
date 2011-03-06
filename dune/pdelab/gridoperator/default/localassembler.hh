#ifndef DUNE_PDELAB_DEFAULT_LOCAL_ASSEMBLER_HH
#define DUNE_PDELAB_DEFAULT_LOCAL_ASSEMBLER_HH

#include <dune/pdelab/gridoperator/default/residualengine.hh>
#include <dune/pdelab/gridoperator/default/patternengine.hh>
#include <dune/pdelab/gridoperator/default/jacobianengine.hh>
#include <dune/pdelab/gridoperator/common/assemblerutilities.hh>
#include <dune/pdelab/common/typetree.hh>

namespace Dune{
  namespace PDELab{

    /**
       \brief The local assembler for DUNE grids

       \tparam GFSU GridFunctionSpace for ansatz functions
       \tparam GFSV GridFunctionSpace for test functions
       \tparam X The solution vector representation type
       \tparam R The residual vector representation type
       \tparam A The jacobian matrix representation type
       \tparam B The matrix backend
       \tparam P The matrix pattern representation type
       \tparam CU   Constraints maps for the individual dofs (trial space)
       \tparam CV   Constraints maps for the individual dofs (test space)

    */
    template<typename GFSU, typename GFSV, typename LOP, 
             typename X, typename R, typename A, typename B, typename P,
             typename CU=Dune::PDELab::EmptyTransformation,
             typename CV=Dune::PDELab::EmptyTransformation>
    class LocalAssembler : public Dune::PDELab::LocalAssemblerBase<B,CU,CV>{
    public:

      //! The base class of this local assembler
      typedef Dune::PDELab::LocalAssemblerBase<B,CU,CV> Base;

      typedef typename Base::Traits Traits;

      //! The current grid view type
      typedef typename GFSU::Traits::GridViewType GridView;

      //! The local operator 
      typedef LOP LocalOperator;

      //! Types related to the finte element spaces
      //! @{
      typedef FiniteElementInterfaceSwitch<typename GFSU::Traits::FiniteElementType> UFiniteElement;
      typedef Dune::BasisInterfaceSwitch<typename UFiniteElement::Basis> UBasis;
      typedef FiniteElementInterfaceSwitch<typename GFSV::Traits::FiniteElementType> VFiniteElement;
      typedef Dune::BasisInterfaceSwitch<typename VFiniteElement::Basis> VBasis;
      typedef typename UBasis::Range URange;
      typedef typename VBasis::Range VRange;
      typedef typename UBasis::RangeField RangeField;
      //! @}

      //! The local function spaces
      //! @{
      // Types of local function spaces
      typedef Dune::PDELab::LocalFunctionSpace<GFSU, Dune::PDELab::TrialSpaceTag> LFSU;
      typedef Dune::PDELab::LocalFunctionSpace<GFSV, Dune::PDELab::TestSpaceTag> LFSV;
      //! @}

      //! The local operators type for real numbers e.g. time
      typedef typename LOP::RealType Real;

      //! The residual representation type
      typedef R Residual;

      //! The solution representation type
      typedef X Solution;

      //! The jacobian representation type
      typedef A Jacobian;

      //! The matrix pattern representation type
      typedef P Pattern;

      //! The local assembler engines
      //! @{
      typedef LocalPatternAssemblerEngine<LocalAssembler> LocalPatternAssemblerEngine;
      typedef LocalResidualAssemblerEngine<LocalAssembler> LocalResidualAssemblerEngine;
      typedef LocalJacobianAssemblerEngine<LocalAssembler> LocalJacobianAssemblerEngine;

      friend class LocalPatternAssemblerEngine<LocalAssembler>;
      friend class LocalResidualAssemblerEngine<LocalAssembler>;
      friend class LocalJacobianAssemblerEngine<LocalAssembler>;
      //! @}

      //! Constructor with empty constraints
      LocalAssembler (LOP & lop_) 
        : lop(lop_),  weight(1.0), 
          pattern_engine(*this), residual_engine(*this), jacobian_engine(*this)
      {}

      //! Constructor for non trivial constraints
      LocalAssembler (LOP & lop_, const CU& cu_, const CV& cv_) 
        : Base(cu_, cv_), 
          lop(lop_),  weight(1.0),
          pattern_engine(*this), residual_engine(*this), jacobian_engine(*this)
      {}

      //! Notifies the local assembler about the current time of
      //! assembling. Should be called before assembling if the local
      //! operator has time dependencies.
      void setTime(Real time_){
        lop.setTime(time_);
      }

      //! Notifies the assembler about the current weight of assembling.
      void setWeight(RangeField weight_){
        weight = weight_;
      }

      //! Time stepping interface
      //! @{
      void preStage (Real time_, int r_) { lop.preStage(time_,r_); }
      void preStep (Real time_, Real dt_, std::size_t stages_){ lop.preStep(time_,dt_,stages_); }
      void postStep (){ lop.postStep(); }
      void postStage (){ lop.postStage(); }
      Real suggestTimestep (Real dt) const{return lop.suggestTimestep(dt); }

      //! @}

      //! Access methods which provid "ready to use" engines
      //! @{

      //! Returns a reference to the requested engine. This engine is
      //! completely configured and ready to use.
      LocalPatternAssemblerEngine & localPatternAssemblerEngine
      (Pattern & p)
      {
        pattern_engine.setPattern(p);
        return pattern_engine;
      }

      //! Returns a reference to the requested engine. This engine is
      //! completely configured and ready to use.
      LocalResidualAssemblerEngine & localResidualAssemblerEngine
      (Residual & r, const Solution & x)
      {
        residual_engine.setResidual(r);
        residual_engine.setSolution(x);
        return residual_engine;
      }

      //! Returns a reference to the requested engine. This engine is
      //! completely configured and ready to use.
      LocalJacobianAssemblerEngine & localJacobianAssemblerEngine
      (Jacobian & a, const Solution & x)
      {
        jacobian_engine.setJacobian(a);
        jacobian_engine.setSolution(x);
        return jacobian_engine;
      }

      //! @}

      //! \brief Query methods for the assembler engines. Theses methods
      //! do not belong to the assembler interface, but simplify the
      //! implementations of query methods in the engines;
      //! @{
      static bool doAlphaVolume() { return LOP::doAlphaVolume; }
      static bool doLambdaVolume() { return LOP::doLambdaVolume; }
      static bool doAlphaSkeleton() { return LOP::doAlphaSkeleton; }
      static bool doLambdaSkeleton() { return LOP::doLambdaSkeleton; }
      static bool doAlphaBoundary()  { return LOP::doAlphaBoundary; }
      static bool doLambdaBoundary() { return LOP::doLambdaBoundary; }
      static bool doAlphaVolumePostSkeleton()  { return LOP::doAlphaVolumePostSkeleton; }
      static bool doLambdaVolumePostSkeleton() { return LOP::doLambdaVolumePostSkeleton; }
      static bool doSkeletonTwoSided()  { assert(!LOP::doSkeletonTwoSided); return false; }
      static bool doPatternVolume()  { return LOP::doPatternVolume; }
      static bool doPatternSkeleton()  { return LOP::doPatternSkeleton; }
      static bool doPatternBoundary()  { return LOP::doPatternBoundary; }
      static bool doPatternVolumePostSkeleton()  { return LOP::doPatternVolumePostSkeleton; }
      //! @}

    private:

      //! The local operator
      LOP & lop;

      //! The current weight of assembling
      RangeField weight;

      //! The engine member objects
      //! @{
      LocalPatternAssemblerEngine  pattern_engine;
      LocalResidualAssemblerEngine residual_engine;
      LocalJacobianAssemblerEngine jacobian_engine;
      //! @}
    };

  };
};
#endif
