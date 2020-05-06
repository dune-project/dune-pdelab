#ifndef DUNE_PDELAB_GRIDOPERATOR_FASTDG_JACOBIANENGINE_HH
#define DUNE_PDELAB_GRIDOPERATOR_FASTDG_JACOBIANENGINE_HH

#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/gridfunctionspace/localvector.hh>
#include <dune/pdelab/gridoperator/common/localmatrix.hh>
#include <dune/pdelab/gridoperator/common/diagonallocalmatrix.hh>
#include <dune/pdelab/gridoperator/common/gridoperatorutilities.hh>
#include <dune/pdelab/gridoperator/common/assemblerutilities.hh>
#include <dune/pdelab/gridoperator/common/localassemblerenginebase.hh>
#include <dune/pdelab/localoperator/callswitch.hh>
#include <dune/pdelab/localoperator/flags.hh>

namespace Dune{
  namespace PDELab{

    /**
       \brief The fast DG local assembler engine for DUNE grids which
       assembles the jacobian matrix

       \tparam LA The local assembler

    */
    template<typename LA>
    class FastDGLocalJacobianAssemblerEngine
      : public LocalAssemblerEngineBase
    {
    public:

      template<typename TrialConstraintsContainer, typename TestConstraintsContainer>
      bool needsConstraintsCaching(const TrialConstraintsContainer& cu, const TestConstraintsContainer& cv)
      {
        return cu.containsNonDirichletConstraints() || cv.containsNonDirichletConstraints();
      }

      //! The type of the wrapping local assembler
      typedef LA LocalAssembler;

      //! The type of the local operator
      typedef typename LA::LocalOperator LOP;

      //! The local function spaces
      typedef typename LA::LFSU LFSU;
      typedef typename LA::LFSUCache LFSUCache;
      typedef typename LFSU::Traits::GridFunctionSpace GFSU;
      typedef typename LA::LFSV LFSV;
      typedef typename LA::LFSVCache LFSVCache;
      typedef typename LFSV::Traits::GridFunctionSpace GFSV;

      //! The type of the jacobian matrix
      typedef typename LA::Traits::Jacobian Jacobian;
      typedef typename Jacobian::ElementType JacobianElement;
      typedef typename Jacobian::template AliasedLocalView<LFSVCache,LFSUCache> JacobianView;

      //! The type of the solution vector
      typedef typename LA::Traits::Solution Solution;
      typedef typename Solution::ElementType SolutionElement;
      typedef typename Solution::template ConstAliasedLocalView<LFSUCache> SolutionView;

      /**
         \brief Constructor

         \param [in] local_assembler_ The local assembler object which
         creates this engine
      */
      FastDGLocalJacobianAssemblerEngine(const LocalAssembler & local_assembler_)
        : local_assembler(local_assembler_),
          lop(local_assembler_.localOperator())
      {}

      //! copy contructor
      /**
       * \note This does not create an exact copy.  Instead it copies the
       *       global views, such that they point to the same global vector.
       *       Local matrices/vectors are constructed freshly without copying
       *       the content.  Views into the local matrices/vectors are
       *       constructed freshly so they reference the local matrices/vector
       *       in the new object, and are given unit weight.  This essentially
       *       creates an engine object that is not currently bound to any
       *       entity, but otherwise behaves like the object it was contructed
       *       from.
       * \note This constructor is needed to implement splitting constructors
       *       in derived classes.
       */
      FastDGLocalJacobianAssemblerEngine
      (const FastDGLocalJacobianAssemblerEngine &other) :
        local_assembler(other.local_assembler), lop(other.lop),
        global_s_s_view(other.global_s_s_view),
        global_s_n_view(other.global_s_n_view),
        global_a_ss_view(other.global_a_ss_view),
        global_a_sn_view(other.global_a_sn_view),
        global_a_ns_view(other.global_a_ns_view),
        global_a_nn_view(other.global_a_nn_view)
      { }

      //! Query methods for the global grid assembler
      //! @{
      bool requireSkeleton() const
      { return local_assembler.doAlphaSkeleton(); }
      bool requireSkeletonTwoSided() const
      { return local_assembler.doSkeletonTwoSided(); }
      bool requireUVVolume() const
      { return local_assembler.doAlphaVolume(); }
      bool requireUVSkeleton() const
      { return local_assembler.doAlphaSkeleton(); }
      bool requireUVBoundary() const
      { return local_assembler.doAlphaBoundary(); }
      bool requireUVVolumePostSkeleton() const
      { return local_assembler.doAlphaVolumePostSkeleton(); }
      //! @}

      //! Public access to the wrapping local assembler
      const LocalAssembler & localAssembler() const
      {
        return local_assembler;
      }

      //! Trial space constraints
      const typename LocalAssembler::Traits::TrialGridFunctionSpaceConstraints& trialConstraints() const
      {
        return localAssembler().trialConstraints();
      }

      //! Test space constraints
      const typename LocalAssembler::Traits::TestGridFunctionSpaceConstraints& testConstraints() const
      {
        return localAssembler().testConstraints();
      }

      //! Set current residual vector. Should be called prior to
      //! assembling.
      void setJacobian(Jacobian & jacobian_)
      {
        global_a_ss_view.attach(jacobian_);
        global_a_sn_view.attach(jacobian_);
        global_a_ns_view.attach(jacobian_);
        global_a_nn_view.attach(jacobian_);
      }

      //! Set current solution vector. Should be called prior to
      //! assembling.
      void setSolution(const Solution & solution_)
      {
        global_s_s_view.attach(solution_);
        global_s_n_view.attach(solution_);
      }

      //! Called immediately after binding of local function space in
      //! global assembler.
      //! @{
      template<typename EG, typename LFSUC, typename LFSVC>
      void onBindLFSUV(const EG & eg, const LFSUC & lfsu_cache, const LFSVC & lfsv_cache)
      {
        global_s_s_view.bind(lfsu_cache);
        global_a_ss_view.bind(lfsv_cache,lfsu_cache);
      }

      template<typename IG, typename LFSUC, typename LFSVC>
      void onBindLFSUVOutside(const IG & ig,
                              const LFSUC & lfsu_s_cache, const LFSVC & lfsv_s_cache,
                              const LFSUC & lfsu_n_cache, const LFSVC & lfsv_n_cache)
      {
        global_s_n_view.bind(lfsu_n_cache);
        global_a_sn_view.bind(lfsv_s_cache,lfsu_n_cache);
        global_a_ns_view.bind(lfsv_n_cache,lfsu_s_cache);
        global_a_nn_view.bind(lfsv_n_cache,lfsu_n_cache);
      }

      //! @}

      //! Called when the local function space is about to be rebound or
      //! discarded
      //! @{
      template<typename EG, typename LFSUC, typename LFSVC>
      void onUnbindLFSUV(const EG & eg, const LFSUC & lfsu_cache, const LFSVC & lfsv_cache)
      {
        global_a_ss_view.unbind();
      }

      template<typename IG, typename LFSUC, typename LFSVC>
      void onUnbindLFSUVOutside(const IG & ig,
                                const LFSUC & lfsu_s_cache, const LFSVC & lfsv_s_cache,
                                const LFSUC & lfsu_n_cache, const LFSVC & lfsv_n_cache)
      {
        global_a_sn_view.unbind();
        global_a_ns_view.unbind();
        global_a_nn_view.unbind();
      }

      //! @}

      //! Methods for loading of the local function's coefficients
      //! @{
      template<typename LFSUC>
      void loadCoefficientsLFSUInside(const LFSUC & lfsu_cache)
      {}
      template<typename LFSUC>
      void loadCoefficientsLFSUOutside(const LFSUC & lfsu_n_cache)
      {}
      template<typename LFSUC>
      void loadCoefficientsLFSUCoupling(const LFSUC & lfsu_c_cache)
      {
        DUNE_THROW(Dune::NotImplemented,"No coupling lfsu_cache available for ");
      }
      //! @}

      //! Notifier functions, called immediately before and after assembling
      //! @{
      void postAssembly(const GFSU& gfsu, const GFSV& gfsv)
      {
        Jacobian& jacobian = global_a_ss_view.container();
        global_s_s_view.detach();
        global_s_n_view.detach();
        global_a_ss_view.detach();
        global_a_sn_view.detach();
        global_a_ns_view.detach();
        global_a_nn_view.detach();

        if(local_assembler.doPostProcessing())
          local_assembler.handle_dirichlet_constraints(gfsv,jacobian);
      }
      //! @}

      //! Assembling methods
      //! @{

      /** Assemble on a given cell without function spaces.

          \return If true, the assembling for this cell is assumed to
          be complete and the assembler continues with the next grid
          cell.
       */
      template<typename EG>
      bool assembleCell(const EG & eg)
      {
        return LocalAssembler::isNonOverlapping && eg.entity().partitionType() != Dune::InteriorEntity;
      }

      template<typename EG, typename LFSUC, typename LFSVC>
      void assembleUVVolume(const EG & eg, const LFSUC & lfsu_cache, const LFSVC & lfsv_cache)
      {
        global_a_ss_view.setWeight(local_assembler.weight());
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doAlphaVolume>::
          jacobian_volume(lop,eg,lfsu_cache.localFunctionSpace(),global_s_s_view,lfsv_cache.localFunctionSpace(),global_a_ss_view);
      }

      template<typename IG, typename LFSUC, typename LFSVC>
      void assembleUVSkeleton(const IG & ig, const LFSUC & lfsu_s_cache, const LFSVC & lfsv_s_cache,
                              const LFSUC & lfsu_n_cache, const LFSVC & lfsv_n_cache)
      {
        global_a_ss_view.setWeight(local_assembler.weight());
        global_a_sn_view.setWeight(local_assembler.weight());
        global_a_ns_view.setWeight(local_assembler.weight());
        global_a_nn_view.setWeight(local_assembler.weight());
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doAlphaSkeleton>::
          jacobian_skeleton(lop,ig,
                            lfsu_s_cache.localFunctionSpace(),global_s_s_view,lfsv_s_cache.localFunctionSpace(),
                            lfsu_n_cache.localFunctionSpace(),global_s_n_view,lfsv_n_cache.localFunctionSpace(),
                            global_a_ss_view, global_a_sn_view,
                            global_a_ns_view, global_a_nn_view);
      }

      template<typename IG, typename LFSUC, typename LFSVC>
      void assembleUVBoundary(const IG & ig, const LFSUC & lfsu_s_cache, const LFSVC & lfsv_s_cache)
      {
        global_a_ss_view.setWeight(local_assembler.weight());
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doAlphaBoundary>::
          jacobian_boundary(lop,ig,lfsu_s_cache.localFunctionSpace(),global_s_s_view,lfsv_s_cache.localFunctionSpace(),global_a_ss_view);
      }

      template<typename IG, typename LFSUC, typename LFSVC>
      static void assembleUVEnrichedCoupling(const IG & ig,
                                             const LFSUC & lfsu_s_cache, const LFSVC & lfsv_s_cache,
                                             const LFSUC & lfsu_n_cache, const LFSVC & lfsv_n_cache,
                                             const LFSUC & lfsu_coupling_cache, const LFSVC & lfsv_coupling_cache)
      {
        DUNE_THROW(Dune::NotImplemented,"Assembling of coupling spaces is not implemented for ");
      }

      template<typename IG, typename LFSVC>
      static void assembleVEnrichedCoupling(const IG & ig,
                                            const LFSVC & lfsv_s_cache,
                                            const LFSVC & lfsv_n_cache,
                                            const LFSVC & lfsv_coupling_cache)
      {
        DUNE_THROW(Dune::NotImplemented,"Assembling of coupling spaces is not implemented for ");
      }

      template<typename EG, typename LFSUC, typename LFSVC>
      void assembleUVVolumePostSkeleton(const EG & eg, const LFSUC & lfsu_cache, const LFSVC & lfsv_cache)
      {
        global_a_ss_view.setWeight(local_assembler.weight());
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doAlphaVolumePostSkeleton>::
          jacobian_volume_post_skeleton(lop,eg,lfsu_cache.localFunctionSpace(),global_s_s_view,lfsv_cache.localFunctionSpace(),global_a_ss_view);
      }

      //! @}

    private:
      //! Reference to the wrapping local assembler object which
      //! constructed this engine
      const LocalAssembler & local_assembler;

      //! Reference to the local operator
      const LOP & lop;

      //! Pointer to the current solution vector for which to assemble
      SolutionView global_s_s_view;
      SolutionView global_s_n_view;

      //! Pointer to the current residual vector in which to assemble
      JacobianView global_a_ss_view;
      JacobianView global_a_sn_view;
      JacobianView global_a_ns_view;
      JacobianView global_a_nn_view;

      //! The local vectors and matrices as required for assembling
      //! @{
      typedef Dune::PDELab::TrialSpaceTag LocalTrialSpaceTag;
      typedef Dune::PDELab::TestSpaceTag LocalTestSpaceTag;

      typedef Dune::PDELab::LocalVector<SolutionElement, LocalTrialSpaceTag> SolutionVector;
      typedef typename std::conditional<
        std::is_base_of<
          lop::DiagonalJacobian,
          LOP
          >::value,
        Dune::PDELab::DiagonalLocalMatrix<JacobianElement>,
        Dune::PDELab::LocalMatrix<JacobianElement>
        >::type JacobianMatrix;

      //! @}

    }; // End of class FastDGLocalJacobianAssemblerEngine

  }
}
#endif
