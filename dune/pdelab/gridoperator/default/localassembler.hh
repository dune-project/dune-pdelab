#ifndef DUNE_PDELAB_DEFAULT_LOCAL_ASSEMBLER_HH
#define DUNE_PDELAB_DEFAULT_LOCAL_ASSEMBLER_HH

#include <dune/typetree/typetree.hh>

#include <dune/pdelab/gridoperator/default/residualengine.hh>
#include <dune/pdelab/gridoperator/default/patternengine.hh>
#include <dune/pdelab/gridoperator/default/jacobianengine.hh>
#include <dune/pdelab/gridoperator/default/jacobianapplyengine.hh>
#include <dune/pdelab/gridoperator/common/assemblerutilities.hh>
#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>

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

      // The GridOperator has to be a friend to modify the do{Pre,Post}Processing flags
      template<typename, typename, typename,
               typename, typename, typename, typename,
               typename, typename, bool>
      friend class GridOperator;

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
      typedef LFSIndexCache<LFSU,CU> LFSUCache;
      typedef LFSIndexCache<LFSV,CV> LFSVCache;

      typedef LFSUCache NoConstraintsLFSUCache DUNE_DEPRECATED_MSG("NoConstraintsLFSUCache is deprecated, use LFSUCache instead and use the runtime interface to disable constraints caching");
      typedef LFSVCache NoConstraintsLFSVCache DUNE_DEPRECATED_MSG("NoConstraintsLFSVCache is deprecated, use LFSVCache instead and use the runtime interface to disable constraints caching");

      //! @}

      //! The local assembler engines
      //! @{
      typedef DefaultLocalPatternAssemblerEngine<DefaultLocalAssembler> LocalPatternAssemblerEngine;
      typedef DefaultLocalResidualAssemblerEngine<DefaultLocalAssembler> LocalResidualAssemblerEngine;
      typedef DefaultLocalJacobianAssemblerEngine<DefaultLocalAssembler> LocalJacobianAssemblerEngine;
      typedef DefaultLocalJacobianApplyAssemblerEngine<DefaultLocalAssembler> LocalJacobianApplyAssemblerEngine;

      friend class DefaultLocalPatternAssemblerEngine<DefaultLocalAssembler>;
      friend class DefaultLocalResidualAssemblerEngine<DefaultLocalAssembler>;
      friend class DefaultLocalJacobianAssemblerEngine<DefaultLocalAssembler>;
      friend class DefaultLocalJacobianApplyAssemblerEngine<DefaultLocalAssembler>;
      //! @}

      //! Constructor with empty constraints
      DefaultLocalAssembler (LOP & lop_, shared_ptr<typename GO::BorderDOFExchanger> border_dof_exchanger)
        : lop(lop_),  weight(1.0), doPreProcessing(true), doPostProcessing(true),
          pattern_engine(*this,border_dof_exchanger), residual_engine(*this), jacobian_engine(*this), jacobian_apply_engine(*this)
        , _reconstruct_border_entries(isNonOverlapping)
      {}

      //! Constructor for non trivial constraints
      DefaultLocalAssembler (LOP & lop_, const CU& cu_, const CV& cv_,
                             shared_ptr<typename GO::BorderDOFExchanger> border_dof_exchanger)
        : Base(cu_, cv_),
          lop(lop_),  weight(1.0), doPreProcessing(true), doPostProcessing(true),
          pattern_engine(*this,border_dof_exchanger), residual_engine(*this), jacobian_engine(*this), jacobian_apply_engine(*this)
        , _reconstruct_border_entries(isNonOverlapping)
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

      bool reconstructBorderEntries() const
      {
        return _reconstruct_border_entries;
      }

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

      //! Returns a reference to the requested engine. This engine is
      //! completely configured and ready to use.
      LocalJacobianApplyAssemblerEngine & localJacobianApplyAssemblerEngine
      (typename Traits::Residual & r, const typename Traits::Solution & x)
      {
        jacobian_apply_engine.setResidual(r);
        jacobian_apply_engine.setSolution(x);
        return jacobian_apply_engine;
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

      //! This method allows to set the behavior with regard to any
      //! preprocessing within the engines. It is called by the
      //! setupGridOperators() method of the GridOperator and should
      //! not be called directly.
      void preProcessing(bool v)
      {
        doPreProcessing = v;
      }

      //! This method allows to set the behavior with regard to any
      //! postprocessing within the engines. It is called by the
      //! setupGridOperators() method of the GridOperator and should
      //! not be called directly.
      void postProcessing(bool v)
      {
        doPostProcessing = v;
      }

    private:

      //! The local operator
      LOP & lop;

      //! The current weight of assembling
      RangeField weight;

      //! Indicates whether this local operator has to perform pre
      //! processing
      bool doPreProcessing;

      //! Indicates whether this local operator has to perform post
      //! processing
      bool doPostProcessing;

      //! The engine member objects
      //! @{
      LocalPatternAssemblerEngine  pattern_engine;
      LocalResidualAssemblerEngine residual_engine;
      LocalJacobianAssemblerEngine jacobian_engine;
      LocalJacobianApplyAssemblerEngine jacobian_apply_engine;
      //! @}

      bool _reconstruct_border_entries;

    };

  }
}
#endif
