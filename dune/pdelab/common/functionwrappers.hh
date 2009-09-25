// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_FUNCTIONWRAPPERS_HH
#define DUNE_PDELAB_FUNCTIONWRAPPERS_HH

#include <vector>

#include "function.hh"

namespace Dune {
  namespace PDELab {

    //! \addtogroup PDELab_Function Function
    //! \ingroup PDELab
    //! \{

    //////////////////////////////////////////////////////////////////////
    //
    //  PointwiseGridFunctionAdapter
    //

    struct Nil {};

    //! A function wrapper which can map a set of gridfunctions through an
    //! PointwiseAdapterEngine
    /**
     * \tparam Engine The type of the engine
     * \tparam FN     The types of the functions.  Currently, N up to 9 is
     *                supported.
     */
    template<typename Engine, typename F0,
             typename F1 = Nil, typename F2 = Nil, typename F3 = Nil,
             typename F4 = Nil, typename F5 = Nil, typename F6 = Nil,
             typename F7 = Nil, typename F8 = Nil, typename F9 = Nil>
    class PointwiseGridFunctionAdapter :
      public GridFunctionInterface<
      typename F0::Traits,
        PointwiseGridFunctionAdapter<
          Engine,
          F0, F1, F2, F3, F4, F5, F6, F7, F8, F9> >
    {
    public:
      typedef typename F0::Traits Traits;

    private:
      const Engine& engine;
      const F0& f0;
      const F1& f1;
      const F2& f2;
      const F3& f3;
      const F4& f4;
      const F5& f5;
      const F6& f6;
      const F7& f7;
      const F8& f8;
      const F9& f9;
      //! hold the number of non-Nil functions
      unsigned size;

      //! increment a sum (but only if type(*f) != Nil via overloading)
      /**
       * \tparam F The type of *f
       *
       * \param sum The variable to increment
       * \param f   Dummy parameter, only its type is of importance to select
       *            the correct method from the overloaded set
       */
      template<typename F>
      static void inc(unsigned &sum, const F* f)
      { ++sum; }
      //! increment a sum (well, don't, since type(*f) == Nil)
      /**
       * \param sum The variable to increment (well, not to imcrement)
       * \param f   Dummy parameter, only its type is of importance to select
       *            the correct method from the overloaded set
       */
      static void inc(unsigned &sum, const Nil* f)
      { }
      //! count the number of non-Nil functions
      /**
       * uses the inc()-methods
       */
      static unsigned count() {
        unsigned sum = 0;
        inc(sum, (F0*)0);
        inc(sum, (F1*)0);
        inc(sum, (F2*)0);
        inc(sum, (F3*)0);
        inc(sum, (F4*)0);
        inc(sum, (F5*)0);
        inc(sum, (F6*)0);
        inc(sum, (F7*)0);
        inc(sum, (F8*)0);
        inc(sum, (F9*)0);
        return sum;
      }

      //! evaluate a function (but only if type(f) != Nil via overloading)
      /**
       * \tparam F The type of f
       *
       * \param f     The function to evaluate
       * \param index The index to write to in the std::vector y.  This index
       *              is incremented after use.
       * \param e     The element at which to evaluate
       * \param x     The local coordinates inside x where to evaluate.
       * \param y     Where to write the result of the evaluation.
       */
      template<typename F>
      static void evaluate(const F &f, unsigned &index,
                           const typename Traits::ElementType& e,
                           const typename Traits::DomainType& x,
                           std::vector<typename Traits::RangeType>& y) {
        f.evaluate(e, x, y[index]);
        ++index;
      }
      //! evaluate a function (well, don't, since type(f) == Nil)
      /**
       * \param f     Dummy argument, only its type Nil is of importance to
       *              select the correct method from the overload set.
       * \param index The index to write to in the std::vector y.  Since this
       *              dummy methid doesn't actually write anythin to y, the
       *              index is not incremented afterwards.
       * \param e     The element at which to evaluate
       * \param x     The local coordinates inside x where to evaluate.
       * \param y     Where to write the result of the evaluation.  Well, its
       *              not actually used in this dummy method
       */
      static void evaluate(const Nil &f, unsigned &index,
                           const typename Traits::ElementType& e,
                           const typename Traits::DomainType& x,
                           std::vector<typename Traits::RangeType>& y)
      { }

      //! hide non-const getGridView() since we cannot implement it
      DUNE_DEPRECATED
	  inline const typename Traits::GridViewType& getGridView ();

    public:
      //! construct a PointwiseGridFunctionAdapter
      /**
       * \param engine_ A reference to the engine.  The referenced engine
       *                object must live at least as long as this adapter
       *                object is used for evaluation.
       * \param fN_     References to the functions.  These referenced funtion
       *                objects must live at least as long as this adapter
       *                object is used for evaluation.  Currently, N up to 9
       *                is supported.
       */
      PointwiseGridFunctionAdapter(const Engine& engine_, const F0& f0_,
                                   const F1& f1_ = Nil(), const F2& f2_ = Nil(), const F3& f3_ = Nil(),
                                   const F4& f4_ = Nil(), const F5& f5_ = Nil(), const F6& f6_ = Nil(),
                                   const F7& f7_ = Nil(), const F8& f8_ = Nil(), const F9& f9_ = Nil())
        : engine(engine_), f0(f0_)
        , f1(f1_), f2(f2_), f3(f3_)
        , f4(f4_), f5(f5_), f6(f6_)
        , f7(f7_), f8(f8_), f9(f9_)
        , size(count())
      {}

	  inline void evaluate (const typename Traits::ElementType& e,
							const typename Traits::DomainType& x,
							typename Traits::RangeType& y) const
	  {
        std::vector<typename Traits::RangeType> in(size);
        unsigned index = 0;
        evaluate(f0, index, e, y, in);
        evaluate(f1, index, e, y, in);
        evaluate(f2, index, e, y, in);
        evaluate(f3, index, e, y, in);
        evaluate(f4, index, e, y, in);
        evaluate(f5, index, e, y, in);
        evaluate(f6, index, e, y, in);
        evaluate(f7, index, e, y, in);
        evaluate(f8, index, e, y, in);
        evaluate(f9, index, e, y, in);
        engine.evaluate(y, in);
	  }

	  inline const typename Traits::GridViewType& getGridView () const
	  {
		return f0.getGridView();
	  }
    };
      
    //! syntactic sugar for easy creation of PointwiseGridFunctionAdapter
    //! objects
    template<typename Engine, typename F0,
             typename F1 = Nil, typename F2 = Nil, typename F3 = Nil,
             typename F4 = Nil, typename F5 = Nil, typename F6 = Nil,
             typename F7 = Nil, typename F8 = Nil, typename F9 = Nil>
    PointwiseGridFunctionAdapter<Engine, F0, F1, F2, F3, F4, F5, F6, F7, F8, F9>
    makePointwiseGridFunctionAdapter(const Engine& engine, const F0& f0,
                                     const F1& f1 = Nil(), const F2& f2 = Nil(), const F3& f3 = Nil(),
                                     const F4& f4 = Nil(), const F5& f5 = Nil(), const F6& f6 = Nil(),
                                     const F7& f7 = Nil(), const F8& f8 = Nil(), const F9& f9 = Nil())
    {
      return PointwiseGridFunctionAdapter
        <Engine,F0,F1,F2,F3,F4,F5,F6,F7,F8,F9>
        (engine,f0,f1,f2,f3,f4,f5,f6,f7,f8,f9);
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
        out = in[1];
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

#endif // DUNE_PDELAB_FUNCTIONWRAPPERS_HH
