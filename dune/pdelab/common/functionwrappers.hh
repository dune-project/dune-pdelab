// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_FUNCTIONWRAPPERS_HH
#define DUNE_PDELAB_FUNCTIONWRAPPERS_HH

#include <vector>

#include <dune/pdelab/common/typetree/nodetags.hh>

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
    template<typename Engine, typename F0,
             typename F1 = TypeTree::EmptyNode,
             typename F2 = TypeTree::EmptyNode,
             typename F3 = TypeTree::EmptyNode,
             typename F4 = TypeTree::EmptyNode,
             typename F5 = TypeTree::EmptyNode,
             typename F6 = TypeTree::EmptyNode,
             typename F7 = TypeTree::EmptyNode,
             typename F8 = TypeTree::EmptyNode,
             typename F9 = TypeTree::EmptyNode>
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
      //! hold the number of non-EmptyNode functions
      unsigned size;

      //! increment a sum (but only if type(*f) != EmptyNode via overloading)
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
      //! increment a sum (well, don't, since type(*f) == EmptyNode)
      /**
       * \param sum The variable to increment (well, not to imcrement)
       * \param f   Dummy parameter, only its type is of importance to select
       *            the correct method from the overloaded set
       */
      static void inc(unsigned &sum, const TypeTree::EmptyNode* f)
      { }
      //! count the number of non-EmptyNode functions
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

      //! \brief evaluate a function (but only if type(f) != EmptyNode via
      //!        overloading)
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
      //! evaluate a function (well, don't, since type(f) == EmptyNode)
      /**
       * \param f     Dummy argument, only its type EmptyNode is of importance
       *              to select the correct method from the overload set.
       * \param index The index to write to in the std::vector y.  Since this
       *              dummy methid doesn't actually write anythin to y, the
       *              index is not incremented afterwards.
       * \param e     The element at which to evaluate
       * \param x     The local coordinates inside x where to evaluate.
       * \param y     Where to write the result of the evaluation.  Well, its
       *              not actually used in this dummy method
       */
      static void evaluate(const TypeTree::EmptyNode &f, unsigned &index,
                           const typename Traits::ElementType& e,
                           const typename Traits::DomainType& x,
                           std::vector<typename Traits::RangeType>& y)
      { }

      //! hide non-const getGridView() since we cannot implement it
      DUNE_DEPRECATED
	  inline const typename Traits::GridViewType& getGridView ();

    public:
#ifdef DOXYGEN
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
      PointwiseGridFunctionAdapter(const Engine& engine_, const F0& f0_, ...)
      {}
#else
      PointwiseGridFunctionAdapter(const Engine& engine_, const F0& f0_,
                                   const F1& f1_ = TypeTree::EmptyNode(),
                                   const F2& f2_ = TypeTree::EmptyNode(),
                                   const F3& f3_ = TypeTree::EmptyNode(),
                                   const F4& f4_ = TypeTree::EmptyNode(),
                                   const F5& f5_ = TypeTree::EmptyNode(),
                                   const F6& f6_ = TypeTree::EmptyNode(),
                                   const F7& f7_ = TypeTree::EmptyNode(),
                                   const F8& f8_ = TypeTree::EmptyNode(),
                                   const F9& f9_ = TypeTree::EmptyNode())
        : engine(engine_), f0(f0_)
        , f1(f1_), f2(f2_), f3(f3_)
        , f4(f4_), f5(f5_), f6(f6_)
        , f7(f7_), f8(f8_), f9(f9_)
        , size(count())
      {}
#endif // DOXYGEN

	  inline void evaluate (const typename Traits::ElementType& e,
							const typename Traits::DomainType& x,
							typename Traits::RangeType& y) const
	  {
        std::vector<typename Traits::RangeType> in(size);
        unsigned index = 0;
        evaluate(f0, index, e, x, in);
        evaluate(f1, index, e, x, in);
        evaluate(f2, index, e, x, in);
        evaluate(f3, index, e, x, in);
        evaluate(f4, index, e, x, in);
        evaluate(f5, index, e, x, in);
        evaluate(f6, index, e, x, in);
        evaluate(f7, index, e, x, in);
        evaluate(f8, index, e, x, in);
        evaluate(f9, index, e, x, in);
        engine.evaluate(y, in);
	  }

	  inline const typename Traits::GridViewType& getGridView () const
	  {
		return f0.getGridView();
	  }
    };

    //! syntactic sugar for easy creation of PointwiseGridFunctionAdapter
    //! objects
    template<typename Engine, typename F0>
    PointwiseGridFunctionAdapter<Engine, F0>
    makePointwiseGridFunctionAdapter(const Engine& engine, const F0& f0)
    {
      return PointwiseGridFunctionAdapter
        <Engine,F0>
        (engine,f0);
    }

    //! syntactic sugar for easy creation of PointwiseGridFunctionAdapter
    //! objects
    template<typename Engine, typename F0, typename F1>
    PointwiseGridFunctionAdapter<Engine, F0, F1>
    makePointwiseGridFunctionAdapter(const Engine& engine, const F0& f0,
                                     const F1& f1)
    {
      return PointwiseGridFunctionAdapter
        <Engine,F0,F1>
        (engine,f0,f1);
    }

    //! syntactic sugar for easy creation of PointwiseGridFunctionAdapter
    //! objects
    template<typename Engine, typename F0, typename F1, typename F2>
    PointwiseGridFunctionAdapter<Engine, F0, F1, F2>
    makePointwiseGridFunctionAdapter(const Engine& engine, const F0& f0,
                                     const F1& f1, const F2& f2)
    {
      return PointwiseGridFunctionAdapter
        <Engine,F0,F1,F2>
        (engine,f0,f1,f2);
    }

    //! syntactic sugar for easy creation of PointwiseGridFunctionAdapter
    //! objects
    template<typename Engine, typename F0, typename F1, typename F2,
             typename F3>
    PointwiseGridFunctionAdapter<Engine, F0, F1, F2, F3>
    makePointwiseGridFunctionAdapter(const Engine& engine, const F0& f0,
                                     const F1& f1, const F2& f2, const F3& f3)
    {
      return PointwiseGridFunctionAdapter
        <Engine,F0,F1,F2,F3>
        (engine,f0,f1,f2,f3);
    }

    //! syntactic sugar for easy creation of PointwiseGridFunctionAdapter
    //! objects
    template<typename Engine, typename F0, typename F1, typename F2,
             typename F3, typename F4>
    PointwiseGridFunctionAdapter<Engine, F0, F1, F2, F3, F4>
    makePointwiseGridFunctionAdapter(const Engine& engine, const F0& f0,
                                     const F1& f1, const F2& f2, const F3& f3,
                                     const F4& f4)
    {
      return PointwiseGridFunctionAdapter
        <Engine,F0,F1,F2,F3,F4>
        (engine,f0,f1,f2,f3,f4);
    }

    //! syntactic sugar for easy creation of PointwiseGridFunctionAdapter
    //! objects
    template<typename Engine, typename F0, typename F1, typename F2,
             typename F3, typename F4, typename F5>
    PointwiseGridFunctionAdapter<Engine, F0, F1, F2, F3, F4, F5>
    makePointwiseGridFunctionAdapter(const Engine& engine, const F0& f0,
                                     const F1& f1, const F2& f2, const F3& f3,
                                     const F4& f4, const F5& f5)
    {
      return PointwiseGridFunctionAdapter
        <Engine,F0,F1,F2,F3,F4,F5>
        (engine,f0,f1,f2,f3,f4,f5);
    }

    //! syntactic sugar for easy creation of PointwiseGridFunctionAdapter
    //! objects
    template<typename Engine, typename F0, typename F1, typename F2,
             typename F3, typename F4, typename F5, typename F6>
    PointwiseGridFunctionAdapter<Engine, F0, F1, F2, F3, F4, F5, F6>
    makePointwiseGridFunctionAdapter(const Engine& engine, const F0& f0,
                                     const F1& f1, const F2& f2, const F3& f3,
                                     const F4& f4, const F5& f5, const F6& f6)
    {
      return PointwiseGridFunctionAdapter
        <Engine,F0,F1,F2,F3,F4,F5,F6>
        (engine,f0,f1,f2,f3,f4,f5,f6);
    }

    //! syntactic sugar for easy creation of PointwiseGridFunctionAdapter
    //! objects
    template<typename Engine, typename F0, typename F1, typename F2,
             typename F3, typename F4, typename F5, typename F6, typename F7>
    PointwiseGridFunctionAdapter<Engine, F0, F1, F2, F3, F4, F5, F6, F7>
    makePointwiseGridFunctionAdapter(const Engine& engine, const F0& f0,
                                     const F1& f1, const F2& f2, const F3& f3,
                                     const F4& f4, const F5& f5, const F6& f6,
                                     const F7& f7)
    {
      return PointwiseGridFunctionAdapter
        <Engine,F0,F1,F2,F3,F4,F5,F6,F7>
        (engine,f0,f1,f2,f3,f4,f5,f6,f7);
    }

    //! syntactic sugar for easy creation of PointwiseGridFunctionAdapter
    //! objects
    template<typename Engine, typename F0, typename F1, typename F2,
             typename F3, typename F4, typename F5, typename F6, typename F7,
             typename F8>
    PointwiseGridFunctionAdapter<Engine, F0, F1, F2, F3, F4, F5, F6, F7, F8>
    makePointwiseGridFunctionAdapter(const Engine& engine, const F0& f0,
                                     const F1& f1, const F2& f2, const F3& f3,
                                     const F4& f4, const F5& f5, const F6& f6,
                                     const F7& f7, const F8& f8)
    {
      return PointwiseGridFunctionAdapter
        <Engine,F0,F1,F2,F3,F4,F5,F6,F7,F8>
        (engine,f0,f1,f2,f3,f4,f5,f6,f7,f8);
    }

    //! syntactic sugar for easy creation of PointwiseGridFunctionAdapter
    //! objects
    template<typename Engine, typename F0, typename F1, typename F2,
             typename F3, typename F4, typename F5, typename F6, typename F7,
             typename F8, typename F9>
    PointwiseGridFunctionAdapter<Engine, F0, F1, F2, F3, F4, F5, F6, F7, F8,
                                 F9>
    makePointwiseGridFunctionAdapter(const Engine& engine, const F0& f0,
                                     const F1& f1, const F2& f2, const F3& f3,
                                     const F4& f4, const F5& f5, const F6& f6,
                                     const F7& f7, const F8& f8, const F9& f9)
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

#endif // DUNE_PDELAB_FUNCTIONWRAPPERS_HH
