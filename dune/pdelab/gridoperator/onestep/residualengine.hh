#ifndef DUNE_PDELAB_GRIDOPERATOR_ONESTEP_RESIDUALENGINE_HH
#define DUNE_PDELAB_GRIDOPERATOR_ONESTEP_RESIDUALENGINE_HH

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
    class OneStepLocalResidualAssemblerEngine
      : public OneStepLocalAssemblerEngineBase<OSLA,
                                               typename OSLA::LocalAssemblerDT0::LocalResidualAssemblerEngine,
                                               typename OSLA::LocalAssemblerDT1::LocalResidualAssemblerEngine
                                               >
    {

      typedef OneStepLocalAssemblerEngineBase<OSLA,
                                              typename OSLA::LocalAssemblerDT0::LocalResidualAssemblerEngine,
                                              typename OSLA::LocalAssemblerDT1::LocalResidualAssemblerEngine
                                              > BaseT;

      using BaseT::la;
      using BaseT::lae0;
      using BaseT::lae1;
      using BaseT::implicit;
      using BaseT::setLocalAssemblerEngineDT0;
      using BaseT::setLocalAssemblerEngineDT1;

    public:
      //! The type of the wrapping local assembler
      typedef OSLA OneStepLocalAssembler;

      //! Types of the subordinate assemblers and engines
      //! @{
      typedef typename OSLA::LocalAssemblerDT0 LocalAssemblerDT0;
      typedef typename OSLA::LocalAssemblerDT1 LocalAssemblerDT1;

      typedef typename LocalAssemblerDT0::LocalResidualAssemblerEngine ResidualEngineDT0;
      typedef typename LocalAssemblerDT1::LocalResidualAssemblerEngine ResidualEngineDT1;
      //! @}

      //! The type of the residual vector
      typedef typename OSLA::Traits::Residual Residual;

      //! The type of the solution vector
      typedef typename OSLA::Traits::Solution Solution;

      //! The type for real numbers
      typedef typename OSLA::Real Real;

      typedef OSLA LocalAssembler;

      /**
         \brief Constructor

         \param [in] local_assembler_ The local assembler object which
         creates this engine
      */
      OneStepLocalResidualAssemblerEngine(const LocalAssembler & local_assembler_)
        : BaseT(local_assembler_)
        , invalid_residual(nullptr)
        , invalid_solution(nullptr)
        , residual_0(invalid_residual)
        , residual_1(invalid_residual)
        , const_residual_0(invalid_residual)
        , const_residual_1(invalid_residual)
        , solution(invalid_solution)
      {}

      //! Set current solution vector. Must be called before
      //! setResidual(). Should be called prior to assembling.
      void setSolution(const Solution & solution_)
      {
        solution = &solution_;
      }

      //! Set current const residual vector. Must be called before
      //! setResidual(). Should be called prior to assembling.
      void setConstResidual(const Residual &const_residual_)
      {
        const_residual_0 = &const_residual_;
        const_residual_1 = &const_residual_;
      }

      //! Set current const residual vector. Should be called prior to
      //! assembling.
      void setResidual(Residual & residual_)
      {
        residual_0 = &residual_;
        residual_1 = &residual_;

        // Initialize the engines of the two wrapped local assemblers
        assert(solution != invalid_solution);
        setLocalAssemblerEngineDT0(la.la0.localResidualAssemblerEngine(*residual_0,*solution));
        setLocalAssemblerEngineDT1(la.la1.localResidualAssemblerEngine(*residual_1,*solution));
      }

      //! Set current const residual vectors. Must be called before
      //! setResidual(). Should be called prior to assembling. Here,
      //! separate vectors are used for the operators corresponding to
      //! the time dervatives of order one and zero.
      void setConstResiduals(const Residual &const_residual_0_, const Residual &const_residual_1_)
      {
        const_residual_0 = &const_residual_0_;
        const_residual_1 = &const_residual_1_;
      }

      //! Set current const residual vectors. Should be called prior
      //! to assembling. Here, separate vectors are used for the
      //! operators corresponding to the time dervatives of order one
      //! and zero.
      void setResiduals(Residual & residual_0_, Residual & residual_1_)
      {
        residual_0 = &residual_0_;
        residual_1 = &residual_1_;

        // Initialize the engines of the two wrapped local assemblers
        assert(solution != invalid_solution);
        setLocalAssemblerEngineDT0(la.la0.localResidualAssemblerEngine(*residual_0,*solution));
        setLocalAssemblerEngineDT1(la.la1.localResidualAssemblerEngine(*residual_1,*solution));
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

        // Update residual vectors with constant part
        assert(const_residual_0 != invalid_residual);
        assert(const_residual_1 != invalid_residual);
        *residual_0 += *const_residual_0;
        if(residual_0 != residual_1){
          assert(const_residual_0 != const_residual_1);
          *residual_1 += *const_residual_1;
        }

        lae0->postAssembly(gfsu,gfsv);
        lae1->postAssembly(gfsu,gfsv);
      }
      //! @}


    private:

      //! Default value indicating an invalid residual pointer
      Residual * const invalid_residual;

      //! Default value indicating an invalid solution pointer
      Solution * const invalid_solution;

      //! Pointer to the current constant part residual vector in
      //! which to assemble the residual corresponding to the operator
      //! representing the time derivative of order zero and one.
      //! @{
      Residual * residual_0;
      Residual * residual_1;
      //! @}

      //! Pointer to the current constant part residual vectors in
      //! which to assemble the residual corresponding to the operator
      //! representing the time derivative of order zero and one.
      //! @{
      const Residual * const_residual_0;
      const Residual * const_residual_1;
      //! @}

      //! Pointer to the current residual vector in which to assemble
      const Solution * solution;

      //! Coefficients of time stepping scheme
      Real b_rr, d_r;

    }; // End of class OneStepLocalResidualAssemblerEngine

  }
}

#endif // DUNE_PDELAB_GRIDOPERATOR_ONESTEP_RESIDUALENGINE_HH
