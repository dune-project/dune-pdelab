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
    template<typename GO, typename LOP, bool nonoverlapping_mode = false>
    class DefaultLocalAssembler : 
      public Dune::PDELab::LocalAssemblerBase<typename GO::Traits::MatrixBackend,
                                              typename GO::Traits::TrialGridFunctionSpaceConstraints,
                                              typename GO::Traits::TestGridFunctionSpaceConstraints>
    {

    public:

      //! The traits class
      typedef Dune::PDELab::LocalAssemblerTraits<GO> Traits;

      //! The local operators type for real numbers e.g. time
      typedef typename Traits::Residual::ElementType RangeField;
      typedef RangeField Real;

      typedef typename Traits::TrialGridFunctionSpace GFSU;
      typedef typename Traits::TestGridFunctionSpace GFSV;

      typedef typename Traits::TrialGridFunctionSpaceConstraints CU;
      typedef typename Traits::TestGridFunctionSpaceConstraints CV;

      //! The base class of this local assembler
      typedef Dune::PDELab::LocalAssemblerBase<typename Traits::MatrixBackend,CU,CV> Base;

      //! The current grid view type
      typedef typename GFSU::Traits::GridViewType GridView;

      //! The local operator
      typedef LOP LocalOperator;

      static const bool isNonOverlapping = nonoverlapping_mode;

      //! The local function spaces
      //! @{
      // Types of local function spaces
      typedef Dune::PDELab::LocalFunctionSpace<GFSU, Dune::PDELab::TrialSpaceTag> LFSU;
      typedef Dune::PDELab::LocalFunctionSpace<GFSV, Dune::PDELab::TestSpaceTag> LFSV;
      //! @}

      //! The local assembler engines
      //! @{
      typedef DefaultLocalPatternAssemblerEngine<DefaultLocalAssembler> LocalPatternAssemblerEngine;
      typedef DefaultLocalResidualAssemblerEngine<DefaultLocalAssembler> LocalResidualAssemblerEngine;
      typedef DefaultLocalJacobianAssemblerEngine<DefaultLocalAssembler> LocalJacobianAssemblerEngine;

      friend class DefaultLocalPatternAssemblerEngine<DefaultLocalAssembler>;
      friend class DefaultLocalResidualAssemblerEngine<DefaultLocalAssembler>;
      friend class DefaultLocalJacobianAssemblerEngine<DefaultLocalAssembler>;
      //! @}

      //! Constructor with empty constraints
      DefaultLocalAssembler (LOP & lop_)
        : lop(lop_),  weight(1.0), doConstraintsPostProcessing(true),
          pattern_engine(*this), residual_engine(*this), jacobian_engine(*this)
      {}

      //! Constructor for non trivial constraints
      DefaultLocalAssembler (LOP & lop_, const CU& cu_, const CV& cv_)
        : Base(cu_, cv_),
          lop(lop_),  weight(1.0), doConstraintsPostProcessing(true),
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
      (typename Traits::MatrixPattern & p)
      {
        pattern_engine.setPattern(p);
        return pattern_engine;
      }

      //! Returns a reference to the requested engine. This engine is
      //! completely configured and ready to use.
      LocalResidualAssemblerEngine & localResidualAssemblerEngine
      (typename Traits::Residual & r, const typename Traits::Solution & x)
      {
        residual_engine.setResidual(r);
        residual_engine.setSolution(x);
        return residual_engine;
      }

      //! Returns a reference to the requested engine. This engine is
      //! completely configured and ready to use.
      LocalJacobianAssemblerEngine & localJacobianAssemblerEngine
      (typename Traits::Jacobian & a, const typename Traits::Solution & x)
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
      static bool doSkeletonTwoSided()  { return LOP::doSkeletonTwoSided; }
      static bool doPatternVolume()  { return LOP::doPatternVolume; }
      static bool doPatternSkeleton()  { return LOP::doPatternSkeleton; }
      static bool doPatternBoundary()  { return LOP::doPatternBoundary; }
      static bool doPatternVolumePostSkeleton()  { return LOP::doPatternVolumePostSkeleton; }
      //! @}

      //! This method allows to set the behavior with regard to the
      //! post processing of the constraints. It is called by the
      //! setupGridOperators() method of the GridOperator and should
      //! not be called directly.
      void constraintsPostProcessing(bool v){
        doConstraintsPostProcessing = v;
      }

    private:

      //! The local operator
      LOP & lop;

      //! The current weight of assembling
      RangeField weight;

      //! Indicates whether this local operator has to perform post
      //! processing with regard to the constraints
      bool doConstraintsPostProcessing;

      //! The engine member objects
      //! @{
      LocalPatternAssemblerEngine  pattern_engine;
      LocalResidualAssemblerEngine residual_engine;
      LocalJacobianAssemblerEngine jacobian_engine;
      //! @}

    };

  }
}
#endif
