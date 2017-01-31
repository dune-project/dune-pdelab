#ifndef DUNE_PDELAB_GRIDOPERATOR_ONESTEP_PATTERNENGINE_HH
#define DUNE_PDELAB_GRIDOPERATOR_ONESTEP_PATTERNENGINE_HH

#include <dune/pdelab/gridoperator/common/assemblerutilities.hh>
#include <dune/pdelab/gridoperator/onestep/enginebase.hh>

namespace Dune{
  namespace PDELab{

    /**
       \brief The local assembler engine for OneStep sub triangulations which
       creates the matrix pattern

       \tparam LA The local udg assembler

    */
    template<typename OSLA>
    class OneStepLocalPatternAssemblerEngine
      : public OneStepLocalAssemblerEngineBase<OSLA,
                                               typename OSLA::LocalAssemblerDT0::LocalPatternAssemblerEngine,
                                               typename OSLA::LocalAssemblerDT1::LocalPatternAssemblerEngine
                                               >
    {

      typedef OneStepLocalAssemblerEngineBase<OSLA,
                                              typename OSLA::LocalAssemblerDT0::LocalPatternAssemblerEngine,
                                              typename OSLA::LocalAssemblerDT1::LocalPatternAssemblerEngine
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

      //! The type of the matrix pattern container
      typedef typename LocalAssembler::Traits::MatrixPattern Pattern;
      typedef Dune::PDELab::LocalSparsityPattern LocalPattern;

      /**
         \brief Constructor

         \param [in] la_ The local assembler object which
         creates this engine
      */
      OneStepLocalPatternAssemblerEngine(const LocalAssembler & la_)
        : BaseT(la_)
        , invalid_pattern(nullptr)
        , pattern(invalid_pattern)
      {}

      //! Set current residual vector. Should be called prior to
      //! assembling.
      void setPattern(Pattern & pattern_)
      {

        // Set pointer to global pattern
        pattern = &pattern_;

        // Initialize the engines of the two wrapped local assemblers
        setLocalAssemblerEngineDT0(la.la0.localPatternAssemblerEngine(pattern_));
        setLocalAssemblerEngineDT1(la.la1.localPatternAssemblerEngine(pattern_));
      }


      //! @name Notification functions
      //! @{
      void preAssembly()
      {
        implicit = la.osp_method->implicit();

        lae0->preAssembly();
        lae1->preAssembly();
      }

      template<typename GFSU, typename GFSV>
      void postAssembly(const GFSU& gfsu, const GFSV& gfsv)
      {
        lae0->postAssembly(gfsu,gfsv);
        lae1->postAssembly(gfsu,gfsv);
      }
      //! @}

    private:

      //! Default value indicating an invalid solution pointer
      Pattern * const invalid_pattern;

      //! Pointer to the current matrix pattern container
      Pattern * pattern;

    }; // End of class OneStepLocalPatternAssemblerEngine

  }
}

#endif // DUNE_PDELAB_GRIDOPERATOR_ONESTEP_PATTERNENGINE_HH
