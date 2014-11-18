#ifndef DUNE_PDELAB_TBB_LOCALASSEMBLER_HH
#define DUNE_PDELAB_TBB_LOCALASSEMBLER_HH

#include <dune/typetree/typetree.hh>

#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>
#include <dune/pdelab/gridoperator/common/assemblerutilities.hh>
#include <dune/pdelab/gridoperator/default/jacobianapplyengine.hh>
#include <dune/pdelab/gridoperator/default/jacobianengine.hh>
#include <dune/pdelab/gridoperator/default/patternengine.hh>
#include <dune/pdelab/gridoperator/tbb/jacobianengine.hh>
#include <dune/pdelab/gridoperator/tbb/residualengine.hh>

namespace Dune{
  namespace PDELab{

    /**
       \brief The local assembler for DUNE grids, with support for tbb

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
    class TBBLocalAssembler :
      public Dune::PDELab::LocalAssemblerBase<typename GO::Traits::MatrixBackend,
                                              typename GO::Traits::TrialGridFunctionSpaceConstraints,
                                              typename GO::Traits::TestGridFunctionSpaceConstraints>
    {

      // The GridOperator has to be a friend to modify the do{Pre,Post}Processing flags
      template<typename, typename, typename, typename, typename,
               typename, typename, typename, typename,
               typename, typename, bool>
      friend class TBBGridOperator;

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

      typedef LFSIndexCache<LFSU,EmptyTransformation> NoConstraintsLFSUCache;
      typedef LFSIndexCache<LFSV,EmptyTransformation> NoConstraintsLFSVCache;

      typedef typename GO::LockManager LockManager;

      //! @}

      //! The local assembler engines
      //! @{
      typedef DefaultLocalPatternAssemblerEngine<TBBLocalAssembler>
        LocalPatternAssemblerEngine;
      typedef TBBLocalResidualAssemblerEngine<TBBLocalAssembler>
        LocalResidualAssemblerEngine;
      typedef TBBLocalJacobianAssemblerEngine<TBBLocalAssembler>
        LocalJacobianAssemblerEngine;
      typedef DefaultLocalJacobianApplyAssemblerEngine<TBBLocalAssembler>
        LocalJacobianApplyAssemblerEngine;

      friend class DefaultLocalPatternAssemblerEngine<TBBLocalAssembler>;
      friend class TBBLocalResidualAssemblerEngine<TBBLocalAssembler>;
      friend class TBBLocalJacobianAssemblerEngine<TBBLocalAssembler>;
      friend class DefaultLocalJacobianApplyAssemblerEngine<TBBLocalAssembler>;
      //! @}

      //! Constructor with empty constraints
      TBBLocalAssembler
      ( LOP & lop_,
        shared_ptr<typename GO::BorderDOFExchanger> border_dof_exchanger,
        const std::shared_ptr<LockManager> &lockManager)
        : lop(lop_),  weight(1.0), doPreProcessing(true), doPostProcessing(true),
          pattern_engine(*this,border_dof_exchanger), residual_engine(*this), jacobian_engine(*this), jacobian_apply_engine(*this)
        , _reconstruct_border_entries(isNonOverlapping)
        , lockManager_(lockManager)
      {}

      //! Constructor for non trivial constraints
      TBBLocalAssembler
      ( LOP & lop_, const CU& cu_, const CV& cv_,
        shared_ptr<typename GO::BorderDOFExchanger> border_dof_exchanger,
        const std::shared_ptr<LockManager> &lockManager)
        : Base(cu_, cv_),
          lop(lop_),  weight(1.0), doPreProcessing(true), doPostProcessing(true),
          pattern_engine(*this,border_dof_exchanger), residual_engine(*this), jacobian_engine(*this), jacobian_apply_engine(*this)
        , _reconstruct_border_entries(isNonOverlapping)
        , lockManager_(lockManager)
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
        residual_engine.setLockManager(*lockManager_);
        return residual_engine;
      }

      //! Returns a reference to the requested engine. This engine is
      //! completely configured and ready to use.
      LocalJacobianAssemblerEngine & localJacobianAssemblerEngine
      (typename Traits::Jacobian & a, const typename Traits::Solution & x)
      {
        jacobian_engine.setJacobian(a);
        jacobian_engine.setSolution(x);
        jacobian_engine.setLockManager(*lockManager_);
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

      //! The lock manager to use for the local engines.
      std::shared_ptr<LockManager> lockManager_;
    };

    /**
       \brief The local assembler for DUNE grids, with support for threading

       This version is for uses with threaded assemblers.  It collects the
       to-be-updated entries in a buffer and writes them once the commit limit
       is reached, obtaining a lock on the vector to prevent races.  At least
       for the engines that support that.

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
    class BatchedTBBLocalAssembler :
      public Dune::PDELab::LocalAssemblerBase<typename GO::Traits::MatrixBackend,
                                              typename GO::Traits::TrialGridFunctionSpaceConstraints,
                                              typename GO::Traits::TestGridFunctionSpaceConstraints>
    {

      // The GridOperator has to be a friend to modify the do{Pre,Post}Processing flags
      template<typename, typename, typename,
               typename, typename, typename, typename,
               typename, typename, bool>
      friend class GridOperator;

      typedef typename GO::Mutex Mutex;

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

      typedef LFSIndexCache<LFSU,EmptyTransformation> NoConstraintsLFSUCache;
      typedef LFSIndexCache<LFSV,EmptyTransformation> NoConstraintsLFSVCache;

      //! @}

      //! The local assembler engines
      //! @{
      typedef DefaultLocalPatternAssemblerEngine<BatchedTBBLocalAssembler>
        LocalPatternAssemblerEngine;
      typedef BatchedTBBLocalResidualAssemblerEngine<BatchedTBBLocalAssembler>
        LocalResidualAssemblerEngine;
      typedef BatchedTBBLocalJacobianAssemblerEngine<BatchedTBBLocalAssembler>
        LocalJacobianAssemblerEngine;
      typedef DefaultLocalJacobianApplyAssemblerEngine
        <BatchedTBBLocalAssembler> LocalJacobianApplyAssemblerEngine;

      friend class DefaultLocalPatternAssemblerEngine
        <BatchedTBBLocalAssembler>;
      friend class BatchedTBBLocalResidualAssemblerEngine
        <BatchedTBBLocalAssembler>;
      friend class BatchedTBBLocalJacobianAssemblerEngine
        <BatchedTBBLocalAssembler>;
      friend class DefaultLocalJacobianApplyAssemblerEngine
        <BatchedTBBLocalAssembler>;
      //! @}

      //! Constructor with empty constraints
      BatchedTBBLocalAssembler
      ( LOP & lop_,
        shared_ptr<typename GO::BorderDOFExchanger> border_dof_exchanger)
        : lop(lop_),  weight(1.0), doPreProcessing(true), doPostProcessing(true),
          pattern_engine(*this,border_dof_exchanger), residual_engine(*this), jacobian_engine(*this), jacobian_apply_engine(*this)
        , _reconstruct_border_entries(isNonOverlapping)
      {}

      //! Constructor for non trivial constraints
      BatchedTBBLocalAssembler
      ( LOP & lop_, const CU& cu_, const CV& cv_,
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
        residual_engine.setResidual(r, residual_mutex);
        residual_engine.setSolution(x);
        return residual_engine;
      }

      //! Returns a reference to the requested engine. This engine is
      //! completely configured and ready to use.
      LocalJacobianAssemblerEngine & localJacobianAssemblerEngine
      (typename Traits::Jacobian & a, const typename Traits::Solution & x)
      {
        jacobian_engine.setJacobian(a, jacobian_mutex);
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

      // must be aquired by the engines that support mutithreading to update
      // the residual vector.  Since this is used by the engines, it is
      // important that it is initialized before them
      Mutex residual_mutex;
      // must be aquired by the engines that support mutithreading to update
      // the jacobian matrix.  Since this is used by the engines, it is
      // important that it is initialized before them
      Mutex jacobian_mutex;

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

    //! The local assembler for DUNE grids, with support for tbb and coloring
    template<typename GO, typename LOP, bool nonoverlapping_mode = false>
    class ColoredTBBLocalAssembler :
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

      typedef LFSIndexCache<LFSU,EmptyTransformation> NoConstraintsLFSUCache;
      typedef LFSIndexCache<LFSV,EmptyTransformation> NoConstraintsLFSVCache;
      //! @}

      //! The local assembler engines
      //! @{
      typedef DefaultLocalPatternAssemblerEngine<ColoredTBBLocalAssembler>
        LocalPatternAssemblerEngine;
      typedef ColoredTBBLocalResidualAssemblerEngine<ColoredTBBLocalAssembler>
        LocalResidualAssemblerEngine;
      typedef ColoredTBBLocalJacobianAssemblerEngine<ColoredTBBLocalAssembler>
        LocalJacobianAssemblerEngine;
      typedef DefaultLocalJacobianApplyAssemblerEngine
        <ColoredTBBLocalAssembler> LocalJacobianApplyAssemblerEngine;

      friend class DefaultLocalPatternAssemblerEngine
        <ColoredTBBLocalAssembler>;
      friend class ColoredTBBLocalResidualAssemblerEngine
        <ColoredTBBLocalAssembler>;
      friend class DefaultLocalResidualAssemblerEngine
        <ColoredTBBLocalAssembler>;
      friend class ColoredTBBLocalJacobianAssemblerEngine
        <ColoredTBBLocalAssembler>;
      friend class DefaultLocalJacobianAssemblerEngine
        <ColoredTBBLocalAssembler>;
      friend class DefaultLocalJacobianApplyAssemblerEngine
        <ColoredTBBLocalAssembler>;
      //! @}

      //! Constructor with empty constraints
      ColoredTBBLocalAssembler
      ( LOP & lop_,
        shared_ptr<typename GO::BorderDOFExchanger> border_dof_exchanger)
        : lop(lop_),  weight(1.0), doPreProcessing(true), doPostProcessing(true),
          pattern_engine(*this,border_dof_exchanger), residual_engine(*this), jacobian_engine(*this), jacobian_apply_engine(*this)
        , _reconstruct_border_entries(isNonOverlapping)
      {}

      //! Constructor for non trivial constraints
      ColoredTBBLocalAssembler
      ( LOP & lop_, const CU& cu_, const CV& cv_,
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
#endif //  DUNE_PDELAB_TBB_LOCALASSEMBLER_HH
