#ifndef DUNE_PDELAB_GRIDOPERATOR_DEFAULT_LOCALASSEMBLER_HH
#define DUNE_PDELAB_GRIDOPERATOR_DEFAULT_LOCALASSEMBLER_HH

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

      //! @}

      //! The local assembler engines
      //! @{
      typedef DefaultLocalPatternAssemblerEngine<DefaultLocalAssembler> LocalPatternAssemblerEngine;
      typedef DefaultLocalResidualAssemblerEngine<DefaultLocalAssembler> LocalResidualAssemblerEngine;
      typedef DefaultLocalJacobianAssemblerEngine<DefaultLocalAssembler> LocalJacobianAssemblerEngine;
      typedef DefaultLocalJacobianApplyAssemblerEngine<DefaultLocalAssembler> LocalJacobianApplyAssemblerEngine;

      // friend declarations such that engines are able to call scatter_jacobian() and add_entry() from base class
      friend class DefaultLocalPatternAssemblerEngine<DefaultLocalAssembler>;
      friend class DefaultLocalJacobianAssemblerEngine<DefaultLocalAssembler>;
      //! @}

      //! Constructor with empty constraints
      DefaultLocalAssembler (LOP & lop, std::shared_ptr<typename GO::BorderDOFExchanger> border_dof_exchanger)
        : lop_(lop),  weight_(1.0), doPreProcessing_(true), doPostProcessing_(true),
          pattern_engine(*this,border_dof_exchanger), residual_engine(*this), jacobian_engine(*this)
        , jacobian_apply_engine(*this)
        , _reconstruct_border_entries(isNonOverlapping)
      {}

      //! Constructor for non trivial constraints
      DefaultLocalAssembler (LOP & lop, const CU& cu, const CV& cv,
                             std::shared_ptr<typename GO::BorderDOFExchanger> border_dof_exchanger)
        : Base(cu, cv),
          lop_(lop),  weight_(1.0), doPreProcessing_(true), doPostProcessing_(true),
          pattern_engine(*this,border_dof_exchanger), residual_engine(*this), jacobian_engine(*this)
        , jacobian_apply_engine(*this)
        , _reconstruct_border_entries(isNonOverlapping)
      {}

      //! get a reference to the local operator
      LOP &localOperator()
      {
        return lop_;
      }

      //! get a reference to the local operator
      const LOP &localOperator() const
      {
        return lop_;
      }

      //! Notifies the local assembler about the current time of
      //! assembling. Should be called before assembling if the local
      //! operator has time dependencies.
      void setTime(Real time_)
      {
        lop_.setTime(time_);
      }

      //! Obtain the weight that was set last
      RangeField weight() const
      {
        return weight_;
      }

      //! Notifies the assembler about the current weight of assembling.
      void setWeight(RangeField weight){
        weight_ = weight;
      }

      //! Time stepping interface
      //! @{
      void preStage (Real time_, int r_) { lop_.preStage(time_,r_); }
      void preStep (Real time_, Real dt_, std::size_t stages_){ lop_.preStep(time_,dt_,stages_); }
      void postStep (){ lop_.postStep(); }
      void postStage (){ lop_.postStage(); }
      Real suggestTimestep (Real dt) const{return lop_.suggestTimestep(dt); }
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
      (const typename Traits::Domain & update, typename Traits::Range & result)
      {
        jacobian_apply_engine.setUpdate(update);
        jacobian_apply_engine.setResult(result);
        return jacobian_apply_engine;
      }

      //! Returns a reference to the requested engine. This engine is
      //! completely configured and ready to use.
      LocalJacobianApplyAssemblerEngine & localJacobianApplyAssemblerEngine
      (const typename Traits::Domain & solution, const typename Traits::Domain & update, typename Traits::Range & result)
      {
        jacobian_apply_engine.setSolution(solution);
        jacobian_apply_engine.setUpdate(update);
        jacobian_apply_engine.setResult(result);
        return jacobian_apply_engine;
      }

      //! @}

      //! \brief Query methods for the assembler engines. Theses methods
      //! do not belong to the assembler interface, but simplify the
      //! implementations of query methods in the engines;
      //! @{
      static constexpr bool doAlphaVolume() { return LOP::doAlphaVolume; }
      static constexpr bool doLambdaVolume() { return LOP::doLambdaVolume; }
      static constexpr bool doAlphaSkeleton() { return LOP::doAlphaSkeleton; }
      static constexpr bool doLambdaSkeleton() { return LOP::doLambdaSkeleton; }
      static constexpr bool doAlphaBoundary()  { return LOP::doAlphaBoundary; }
      static constexpr bool doLambdaBoundary() { return LOP::doLambdaBoundary; }
      static constexpr bool doAlphaVolumePostSkeleton()  { return LOP::doAlphaVolumePostSkeleton; }
      static constexpr bool doLambdaVolumePostSkeleton() { return LOP::doLambdaVolumePostSkeleton; }
      static constexpr bool doSkeletonTwoSided()  { return LOP::doSkeletonTwoSided; }
      static constexpr bool doPatternVolume()  { return LOP::doPatternVolume; }
      static constexpr bool doPatternSkeleton()  { return LOP::doPatternSkeleton; }
      static constexpr bool doPatternBoundary()  { return LOP::doPatternBoundary; }
      static constexpr bool doPatternVolumePostSkeleton()  { return LOP::doPatternVolumePostSkeleton; }
      static constexpr bool isLinear() { return LOP::isLinear;}
      //! @}

      //! Query whether to do preprocessing in the engines
      /**
       * This method is used by the engines.
       */
      bool doPreProcessing() const { return doPreProcessing_; }

      //! This method allows to set the behavior with regard to any
      //! preprocessing within the engines. It is called by the
      //! setupGridOperators() method of the GridOperator and should
      //! not be called directly.
      void preProcessing(bool v)
      {
        doPreProcessing_ = v;
      }

      //! Query whether to do postprocessing in the engines
      /**
       * This method is used by the engines.
       */
      bool doPostProcessing() const { return doPostProcessing_; }

      //! This method allows to set the behavior with regard to any
      //! postprocessing within the engines. It is called by the
      //! setupGridOperators() method of the GridOperator and should
      //! not be called directly.
      void postProcessing(bool v)
      {
        doPostProcessing_ = v;
      }

    private:

      //! The local operator
      LOP & lop_;

      //! The current weight of assembling
      RangeField weight_;

      //! Indicates whether this local operator has to perform pre
      //! processing
      bool doPreProcessing_;

      //! Indicates whether this local operator has to perform post
      //! processing
      bool doPostProcessing_;

      //! The engine member objects
      //! @{
      LocalPatternAssemblerEngine  pattern_engine;
      LocalResidualAssemblerEngine residual_engine;
      LocalJacobianAssemblerEngine jacobian_engine;
      LocalJacobianApplyAssemblerEngine jacobian_apply_engine;
      //! @}

      bool _reconstruct_border_entries;
    };

  } // end namespace PDELab
} // end namespace Dune
#endif
