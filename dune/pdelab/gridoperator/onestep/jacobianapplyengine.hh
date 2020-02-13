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

      //! The type of the result vector
      typedef typename OSLA::Traits::Range Range;

      //! The type of the solution vector
      typedef typename OSLA::Traits::Domain Domain;

      //! The type for real numbers
      typedef typename OSLA::Real Real;

      /**
         \brief Constructor

         \param[in] local_assembler_ The local assembler object which creates this engine.
      */
      OneStepLocalJacobianApplyAssemblerEngine(LocalAssembler& local_assembler_)
        : BaseT(local_assembler_)
        , invalid_result(nullptr)
        , invalid_solution(nullptr)
        , invalid_update(nullptr)
        , result(invalid_result)
        , solution(invalid_solution)
        , update(invalid_update)
      {}

      //! Set current solution vector. Must be called before
      //! setResult(). Should be called prior to assembling.
      void setSolution(const Domain& solution_)
      {
        solution = &solution_;
      }

      //! Set current update vector. Must be called before
      //! setResult(). Should be called prior to assembling.
      void setUpdate(const Domain& update_)
      {
        update = &update_;
      }

      //! Set current result vector. Should be called prior to
      //! assembling.
      void setResult(Range& result_)
      {
        result = &result_;

        // Initialize the engines of the two wrapped local assemblers
        assert(update != invalid_update);
        setLocalAssemblerEngineDT0
          (la.la0.localJacobianApplyAssemblerEngine(*solution,*update,*result));
        setLocalAssemblerEngineDT1
          (la.la1.localJacobianApplyAssemblerEngine(*solution,*update,*result));
      }

      //! When multiple engines are combined in one assembling
      //! procedure, this method allows to reset the weights which may
      //! have been changed by the other engines.
      void setWeights()
      {
        la.la0.setWeight(b_rr * la.dt_factor0);
        la.la1.setWeight(la.dt_factor1);
      }

      //! Notifier functions, called immediately before and after assembling
      //! @{
      void preAssembly()
      {
        lae0->preAssembly();
        lae1->preAssembly();

        // Extract the coefficients of the time step scheme
        b_rr = la.osp_method->b(la.stage,la.stage);
        d_r = la.osp_method->d(la.stage);

        // Here we only want to know whether this stage is implicit
        using std::abs;
        implicit = abs(b_rr) > 1e-6;

        // prepare local operators for stage
        la.la0.setTime(la.time + d_r * la.dt);
        la.la1.setTime(la.time + d_r * la.dt);

        setWeights();
      }

      template<typename GFSU, typename GFSV>
      void postAssembly(const GFSU& gfsu, const GFSV& gfsv)
      {
        lae0->postAssembly(gfsu,gfsv);
        lae1->postAssembly(gfsu,gfsv);
      }
      //! @}

    private:

      //! Default value indicating an invalid result pointer
      Range * const invalid_result;

      //! Default value indicating an invalid solution pointer
      Domain * const invalid_solution;

      //! Default value indicating an invalid update pointer
      Domain * const invalid_update;

      //! Pointer to the current constant part result vector in
      //! which to assemble
      Range * result;

      //! Pointer to the current result vector in which to assemble
      const Domain * solution;

      //! Pointer to the current update vector
      const Domain * update;

      //! Coefficients of time stepping scheme
      Real b_rr, d_r;

    }; // end class OneStepLocalJacobianApplyAssemblerEngine

  } // end namespace PDELab
} // end Dune
#endif
