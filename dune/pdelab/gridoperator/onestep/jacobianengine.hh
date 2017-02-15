#ifndef DUNE_PDELAB_GRIDOPERATOR_ONESTEP_JACOBIANENGINE_HH
#define DUNE_PDELAB_GRIDOPERATOR_ONESTEP_JACOBIANENGINE_HH

#include <dune/pdelab/gridoperator/onestep/enginebase.hh>
#include <cmath>

namespace Dune{
  namespace PDELab{

    /**
       \brief The local assembler engine for one step methods which
       assembles the residual vector

       \tparam LA The local one step assembler

    */
    template<typename OSLA>
    class OneStepLocalJacobianAssemblerEngine
      : public OneStepLocalAssemblerEngineBase<OSLA,
                                               typename OSLA::LocalAssemblerDT0::LocalJacobianAssemblerEngine,
                                               typename OSLA::LocalAssemblerDT1::LocalJacobianAssemblerEngine
                                               >
    {

      typedef OneStepLocalAssemblerEngineBase<OSLA,
                                              typename OSLA::LocalAssemblerDT0::LocalJacobianAssemblerEngine,
                                              typename OSLA::LocalAssemblerDT1::LocalJacobianAssemblerEngine
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

      typedef typename LocalAssemblerDT0::LocalJacobianAssemblerEngine JacobianEngineDT0;
      typedef typename LocalAssemblerDT1::LocalJacobianAssemblerEngine JacobianEngineDT1;

      //! The type of the residual vector
      typedef typename OSLA::Traits::Jacobian Jacobian;

      //! The type of the solution vector
      typedef typename OSLA::Traits::Solution Solution;

      //! The type for real numbers
      typedef typename OSLA::Real Real;

      /**
         \brief Constructor

         \param [in] local_assembler_ The local assembler object which
         creates this engine
      */
      OneStepLocalJacobianAssemblerEngine(const LocalAssembler & local_assembler_)
        : BaseT(local_assembler_),
          invalid_jacobian(nullptr),
          invalid_solution(nullptr),
          jacobian(invalid_jacobian), solution(invalid_solution)
      {}


      //! Set current solution vector. Must be called before
      //! setResidual(). Should be called prior to assembling.
      void setSolution(const Solution & solution_)
      {
        solution = &solution_;
      }

      //! Set current residual vector. Should be called prior to
      //! assembling.
      void setJacobian(Jacobian & jacobian_)
      {
        jacobian = &jacobian_;

        assert(solution != invalid_solution);

        // Initialize the engines of the two wrapped local assemblers
        setLocalAssemblerEngineDT0(la.la0.localJacobianAssemblerEngine(*jacobian,*solution));
        setLocalAssemblerEngineDT1(la.la1.localJacobianAssemblerEngine(*jacobian,*solution));
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

      //! Default value indicating an invalid residual pointer
      Jacobian * const invalid_jacobian;

      //! Default value indicating an invalid solution pointer
      Solution * const invalid_solution;

      //! Pointer to the current constant part residual vector in
      //! which to assemble
      Jacobian * jacobian;

      //! Pointer to the current residual vector in which to assemble
      const Solution * solution;

      //! Coefficients of time stepping scheme
      Real b_rr, d_r;

    }; // End of class OneStepLocalJacobianAssemblerEngine

  }
}

#endif // DUNE_PDELAB_GRIDOPERATOR_ONESTEP_JACOBIANENGINE_HH
