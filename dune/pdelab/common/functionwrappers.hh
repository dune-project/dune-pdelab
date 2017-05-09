// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_COMMON_FUNCTIONWRAPPERS_HH
#define DUNE_PDELAB_COMMON_FUNCTIONWRAPPERS_HH

#include <vector>
#include <tuple>

#include <dune/typetree/nodetags.hh>

#include "function.hh"

namespace Dune {
  namespace PDELab {

    //! \ingroup PDELab_FunctionAdapters
    //! \{

    //////////////////////////////////////////////////////////////////////
    //
    //  PointwiseGridFunctionAdapter
    //

    //! A function wrapper which can map a set of gridfunctions through an
    //! PointwiseAdapterEngine
    /**
     * \tparam Engine The type of the engine
     * \tparam FN     The types of the functions.  Currently, N up to 9 is
     *                supported.
     */
    template<typename Engine, typename F0, typename... Functions>
    class PointwiseGridFunctionAdapter :
      public GridFunctionBase<
      typename F0::Traits,
      PointwiseGridFunctionAdapter<Engine, F0, Functions...> >
    {
    public:
      typedef typename F0::Traits Traits;

    private:
      const Engine& engine;
      std::tuple<const F0*, const Functions*...> storage;

      //! \brief evaluate a particular function from storage
      /**
       * \tparam I    integral constant which we use to select the correct function
       *
       * \param I     the integral constant, i.e. index
       * \param e     The element at which to evaluate
       * \param x     The local coordinates inside x where to evaluate.
       * \param y     Where to write the result of the evaluation.
       */
      template<unsigned int I, unsigned int N>
      void evaluate(
        const typename Traits::ElementType& e,
        const typename Traits::DomainType& x,
        std::vector<typename Traits::RangeType>& y,
        std::integral_constant<unsigned int, I>,
        std::integral_constant<unsigned int, N>,
        std::integral_constant<bool, true>) const {
        std::get<I>(storage)->evaluate(e, x, y[I]);
        evaluate(e,x,y,
          std::integral_constant<unsigned int, I+1>(),
          std::integral_constant<unsigned int, N>(),
          std::integral_constant<bool, (I+1<N)>()
          );
      }

      //! end of evaluate recursion
      template<unsigned int I>
      void evaluate(
        const typename Traits::ElementType& e,
        const typename Traits::DomainType& x,
        std::vector<typename Traits::RangeType>& y,
        std::integral_constant<unsigned int, I>,
        std::integral_constant<unsigned int, I>,
        std::integral_constant<bool, false>) const
      {}

    public:
      //! construct a PointwiseGridFunctionAdapter
      /**
       * \param engine_ A reference to the engine.  The referenced engine
       *                object must live at least as long as this adapter
       *                object is used for evaluation.
       * \param f0_     Reference to the first function.
       * \param ...     References to the other functions. These referenced function
       *                objects must live at least as long as this adapter
       *                object is used for evaluation. Currently, up to 9 functions
       *                are supported.
       */
      template<typename... F>
      PointwiseGridFunctionAdapter(const Engine& engine_, const F&... functions)
        : engine(engine_), storage(&functions...)
      {}

	  inline void evaluate (const typename Traits::ElementType& e,
							const typename Traits::DomainType& x,
							typename Traits::RangeType& y) const
	  {
        static const unsigned int N = sizeof...(Functions)+1;
        std::vector<typename Traits::RangeType> in(N);
        evaluate(e, x, in,
          std::integral_constant<unsigned int, 0>(),
          std::integral_constant<unsigned int, N>(),
          std::integral_constant<bool, true>()
          );
        engine.evaluate(y, in);
	  }

	  inline const typename Traits::GridViewType& getGridView () const
	  {
        return std::get<0>(storage)->getGridView();
	  }

    };

    //! syntactic sugar for easy creation of PointwiseGridFunctionAdapter
    //! objects
    template<typename Engine, typename... Functions>
    PointwiseGridFunctionAdapter<Engine, Functions...>
    makePointwiseGridFunctionAdapter(const Engine& engine, const Functions&... f)
    {
      return PointwiseGridFunctionAdapter
        <Engine,Functions...>
        (engine,f...);
    }

    //////////////////////////////////////////////////////////////////////
    //
    //  AdapterEngines
    //

    //! Interface of a pointwise adapter engine
    class PointwiseAdapterEngine
    {
    public:

      //! calculate the adapted value from a set of input values
      /**
       * \tparam Domain type of input value
       * \tparam Range  type of output value
       *
       * \param out Where to store the output value
       * \param in  The list of input values
       */
      template<typename Domain, typename Range>
      void evaluate(Range& out,
                    const std::vector<Domain>& in) const;
    };

    //! Scale the output value
    /**
     * \tparam S type of the scaling factor
     *
     * \implements PointwiseAdapterEngine
     */
    template<typename S>
    class PointwiseScaleAdapterEngine
    {
      S scale;

    public:

      //! create a PointwiseScaleAdapterEngine
      /**
       * \param scale_ The scaling factor
       */
      PointwiseScaleAdapterEngine(const S scale_)
        : scale(scale_)
      {}

      //! \copydoc PointwiseAdapterEngine::evaluate
      template<typename Domain, typename Range>
      void evaluate(Range& out,
                    const std::vector<Domain>& in) const {
        assert(in.size() == 1);
        out = in[0];
        out *= scale;
      }
    };
    //! syntactic sugar to create a PointwiseScaleAdapterEngine
    /**
     * \relates PointwiseScaleAdapterEngine
     */
    template<typename S>
    PointwiseScaleAdapterEngine<S>
    makePointwiseScaleAdapterEngine(const S scale) {
      return PointwiseScaleAdapterEngine<S>(scale);
    }

    //! Sum all terms
    /**
     * \tparam S type of the scaling factor
     *
     * \implements PointwiseAdapterEngine
     */
    class PointwiseSumAdapterEngine
    {
    public:

      //! \copydoc PointwiseAdapterEngine::evaluate
      template<typename Domain, typename Range>
      void evaluate(Range& out,
                    const std::vector<Domain>& in) const {
        out = 0;
        for(unsigned i = 0; i < in.size(); ++i)
          out += in[i];
      }
    };

    //! \} Function

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_COMMON_FUNCTIONWRAPPERS_HH
