// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_FUNCTION_CONST_HH
#define DUNE_PDELAB_FUNCTION_CONST_HH

#include <dune/pdelab/common/function.hh>

namespace Dune {
  namespace PDELab {

    //! GridFunction returning a constant value everywhere
    /**
     * \tparam GV       The type of the GridView
     * \tparam RF       The type of the range field
     * \tparam dimRange The dimension of the Range
     */
    template<typename GV, typename RF, unsigned dimR = 1>
    class ConstGridFunction
      : public AnalyticGridFunctionBase<
          AnalyticGridFunctionTraits<GV,RF,dimR>,
          ConstGridFunction<GV,RF,dimR>
        >,
        public InstationaryFunctionDefaults
    {
    public:
      typedef AnalyticGridFunctionTraits<GV,RF,dimR> Traits;
      typedef AnalyticGridFunctionBase<
        Traits,
        ConstGridFunction<GV,RF, dimR> > BaseT;

      //! Contruct a Const GridFunction
      /**
       * \param gv   The GridView to use.  It is passed as a reference to
       *             AnalyticGridFunctionBase (look there for the requirements
       *             of this argument).
       * \param val_ The value tu return on evaluation.  This class stores a
       *             copy of that value.
       */
      ConstGridFunction(const GV& gv,
                        const typename Traits::RangeType& val_ = 1)
        : BaseT(gv)
        , val(val_)
      {}

      //! evaluate the function globally
      /**
       * \param x  Position in global coordinates where to evaluate.
       * \param y  The resulting value.
       */
      inline void
      evaluateGlobal(const typename Traits::DomainType& x,
                     typename Traits::RangeType& y) const
      {
        y = val;
      }

    private:
      typename Traits::RangeType val;
    };

    //! BoundaryGridFunction returning a constant value everywhere
    /**
     * \tparam GV       The type of the GridView
     * \tparam RF       The type of the range field
     * \tparam dimRange The dimension of the Range
     */
    template<typename GV, typename RF, unsigned dimR = 1>
    class ConstBoundaryGridFunction
      : public BoundaryGridFunctionBase<
          BoundaryGridFunctionTraits<
            GV,
            RF,dimR,FieldVector<RF,dimR> >,
          ConstBoundaryGridFunction<GV, RF, dimR> >
    {
    public:
      //! export Traits class
      typedef BoundaryGridFunctionTraits<
        GV,
        RF,dimR,FieldVector<RF,dimR> > Traits;

    private:
      typedef BoundaryGridFunctionBase<
        Traits,
        ConstBoundaryGridFunction<GV,RF,dimR> > BaseT;

    public:
      //! Contruct a ConstBoundaryGridFunction
      /**
       * \param gv_  The GridView to use.  This class stores a copy of it.
       * \param val_ The value tu return on evaluation.  This class stores a
       *             copy of that value.
       */
      ConstBoundaryGridFunction(const GV& gv_,
                                const typename Traits::RangeType& val_ = 1)
        : gv(gv_)
        , val(val_)
      {}

      //! evaluate the function
      /**
       * \tparam I Type of intersection
       *
       * \param ig IntersectionGeometry of this boundary intersection.
       * \param x  Position in local coordinates where to evaluate.
       * \param y  The resulting value.
       */
      template<typename I>
      inline void
      evaluate(const IntersectionGeometry<I>& ig,
               const typename Traits::DomainType& x,
               typename Traits::RangeType& y) const
      {
        y = val;
      }

      //! get a reference to the GridView
      inline const GV& getGridView () const
      {
        return gv;
      }

    private:
      const GV gv;
      const typename Traits::RangeType val;
    };

  } // namspace PDELab
} // namspace Dune

#endif // DUNE_PDELAB_FUNCTION_CONST_HH
