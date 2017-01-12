// -*- tab-width: 2; indent-tabs-mode: nil -*-
// vi: set et ts=2 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDOPERATOR_ONESTEP_JACOBIANAPPLYENGINE_HH
#define DUNE_PDELAB_GRIDOPERATOR_ONESTEP_JACOBIANAPPLYENGINE_HH

#include <dune/pdelab/gridoperator/onestep/enginebase.hh>
#include <cmath>

namespace Dune{
  namespace PDELab{

    /**
       \brief The local assembler engine for one-step methods which
       applies the jacobian without explicitly assembling it.

       \tparam OSLA The one-step local assembler.

    */
    template<typename OSLA>
    class OneStepLocalJacobianApplyAssemblerEngine
      : public OneStepLocalAssemblerEngineBase<OSLA,
                                               typename OSLA::LocalAssemblerDT0::LocalJacobianApplyAssemblerEngine,
                                               typename OSLA::LocalAssemblerDT1::LocalJacobianApplyAssemblerEngine
                                               >
    {
      typedef OneStepLocalAssemblerEngineBase<OSLA,
                                              typename OSLA::LocalAssemblerDT0::LocalJacobianApplyAssemblerEngine,
                                              typename OSLA::LocalAssemblerDT1::LocalJacobianApplyAssemblerEngine
                                              > BaseT;

      using BaseT::la;
      using BaseT::lae0;
      using BaseT::lae1;
      using BaseT::implicit;
      using BaseT::setLocalAssemblerEngineDT0;
      using BaseT::setLocalAssemblerEngineDT1;
    public:
      //! The type of the wrapping local assembler
      typedef OSLA LocalAssembler;

      typedef typename OSLA::LocalAssemblerDT0 LocalAssemblerDT0;
      typedef typename OSLA::LocalAssemblerDT1 LocalAssemblerDT1;

      typedef typename LocalAssemblerDT0::LocalJacobianApplyAssemblerEngine JacobianEngineDT0;
      typedef typename LocalAssemblerDT1::LocalJacobianApplyAssemblerEngine JacobianEngineDT1;

      //! The type of the residual vector
      typedef typename OSLA::Traits::Residual Residual;

      //! The type of the solution vector
      typedef typename OSLA::Traits::Solution Solution;

      //! The type for real numbers
      typedef typename OSLA::Real Real;

      /**
         \brief Constructor

         \param[in] local_assembler_ The local assembler object which creates this engine.
      */
      OneStepLocalJacobianApplyAssemblerEngine(LocalAssembler& local_assembler_)
        : BaseT(local_assembler_)
        , invalid_residual(nullptr)
        , invalid_solution(nullptr)
        , residual(invalid_residual), solution(invalid_solution)
      {}

      //! Set current solution vector. Must be called before
      //! setResidual(). Should be called prior to assembling.
      void setSolution(const Solution& solution_)
      {
        solution = &solution_;
      }

      //! Set current residual vector. Should be called prior to
      //! assembling.
      void setResidual(Residual& residual_)
      {
        residual = &residual_;

        // Initialize the engines of the two wrapped local assemblers
        assert(solution != invalid_solution);
        setLocalAssemblerEngineDT0
          (la.child0().localJacobianApplyAssemblerEngine(*residual,*solution));
        setLocalAssemblerEngineDT1
          (la.child1().localJacobianApplyAssemblerEngine(*residual,*solution));
      }

      //! When multiple engines are combined in one assembling
      //! procedure, this method allows to reset the weights which may
      //! have been changed by the other engines.
      void setWeights()
      {
        la.child0().setWeight(b_rr * la.dt_factor0());
        la.child1().setWeight(la.dt_factor1());
      }

      //! Notifier functions, called immediately before and after assembling
      //! @{
      void preAssembly()
      {
        lae0->preAssembly();
        lae1->preAssembly();

        // Extract the coefficients of the time step scheme
        b_rr = la.method().b(la.stage(),la.stage());
        d_r = la.method().d(la.stage());

        // Here we only want to know whether this stage is implicit
        using std::abs;
        implicit = abs(b_rr) > 1e-6;

        // prepare local operators for stage
        la.child0().setTime(la.timeAtStage());
        la.child1().setTime(la.timeAtStage());

        setWeights();
      }

      template<typename GFSU, typename GFSV>
      void postAssembly(const GFSU& gfsu, const GFSV& gfsv)
      {
        lae0->postAssembly(gfsu,gfsv);
        lae1->postAssembly(gfsu,gfsv);
      }
      //! @}

      //! @name Multithreading support
      //! @{

      //! initialize another engine from this engine.
      /**
       * This is called instead of \c other.preAssembly().  It is an
       * oppertunity to do any setup that is needed at the begin of a thread.
       * It can also be used to copy or split data from \c *this to \c other.
       */
      void split(OneStepLocalJacobianApplyAssemblerEngine &other)
      {
        BaseT::split(other);
        b_rr = other.b_rr;
        d_r = other.d_r;
      }

      //! @}

    private:

      //! Default value indicating an invalid residual pointer
      Residual * const invalid_residual;

      //! Default value indicating an invalid solution pointer
      Solution * const invalid_solution;

      //! Pointer to the current constant part residual vector in
      //! which to assemble
      Residual * residual;

      //! Pointer to the current residual vector in which to assemble
      const Solution * solution;

      //! Coefficients of time stepping scheme
      Real b_rr, d_r;

    }; // end class OneStepLocalJacobianApplyAssemblerEngine

  } // end namespace PDELab
} // end Dune
#endif
