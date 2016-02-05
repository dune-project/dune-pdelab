// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_FUNCTION_LAMBDAADAPTER_HH
#define DUNE_PDELAB_FUNCTION_LAMBDAADAPTER_HH

#include<utility>

namespace Dune {
  namespace PDELab {

    /******************************************************/
    /** \brief Adapter for globally defined functions     */
    /******************************************************/
    template<typename GV, typename RF, int n, typename F>
    class GlobalLambdaToGridFunctionAdapter
      : public Dune::PDELab::GridFunctionBase<Dune::PDELab::
                                              GridFunctionTraits<GV,RF,n,Dune::FieldVector<RF,n> >,
                                              GlobalLambdaToGridFunctionAdapter<GV,RF,n,F> >
    {
      GV gv;
      F f;
    public:
      typedef Dune::PDELab::
      GridFunctionTraits<GV,RF,n,Dune::FieldVector<RF,n> > Traits;

      //! construct from grid view
      GlobalLambdaToGridFunctionAdapter (const GV& gv_, F f_) : gv(gv_), f(f_) {}

      //! get a reference to the grid view
      inline const GV& getGridView () {return gv;}

      //! evaluate extended function on element
      inline void evaluate (const typename Traits::ElementType& e,
                            const typename Traits::DomainType& xl,
                            typename Traits::RangeType& y) const
      {
        const int dim = Traits::GridViewType::Grid::dimension;
        typename Traits::DomainType xg = e.geometry().global(xl);
        y = f(xg);
        return;
      }
    };

    template<typename T>
    struct LambdaAdapterGetDim {
      enum {dim=1};
    };

    template<typename T, int n>
    struct LambdaAdapterGetDim< FieldVector<T,n> > {
      enum {dim=n};
    };

    template<typename T>
    struct LambdaAdapterGetRangeFieldType {
      typedef T Type;
    };

    template<typename T, int n>
    struct LambdaAdapterGetRangeFieldType< FieldVector<T,n> > {
      typedef T Type;
    };

    /** \brief return a PDELab GridFunction defined by a lambda  */
    template<typename GV, typename LAMBDA>
    auto makeGridFunctionFromGlobalLambda (const GV& gv, LAMBDA lambda)
    {
      typedef typename GV::template Codim<0>::Entity::Geometry::GlobalCoordinate X;
      X x;
      typedef decltype(lambda(x)) ReturnType;
      typedef typename LambdaAdapterGetRangeFieldType<ReturnType>::Type RF;
      const int dim = LambdaAdapterGetDim<ReturnType>::dim;
      typedef GlobalLambdaToGridFunctionAdapter<GV,RF,dim,LAMBDA> TheType;
      return TheType(gv,lambda);
    }

    /** \brief return a PDELab GridFunction defined by a parameter class and a lambda  */
    template<typename GV, typename RF, int n, typename F>
    class LocalLambdaToGridFunctionAdapter
      : public Dune::PDELab::GridFunctionBase<Dune::PDELab::
                                              GridFunctionTraits<GV,RF,n,Dune::FieldVector<RF,n> >,
                                              LocalLambdaToGridFunctionAdapter<GV,RF,n,F> >
    {
      GV gv;
      F f;
    public:
      typedef Dune::PDELab::
      GridFunctionTraits<GV,RF,n,Dune::FieldVector<RF,n> > Traits;

      //! construct from grid view
      LocalLambdaToGridFunctionAdapter (const GV& gv_, F f_) : gv(gv_), f(f_) {}

      //! get a reference to the grid view
      inline const GV& getGridView () {return gv;}

      //! evaluate extended function on element
      inline void evaluate (const typename Traits::ElementType& e,
                            const typename Traits::DomainType& xl,
                            typename Traits::RangeType& y) const
      {
        y = f(e,xl);
        return;
      }
    };

    /** \brief return a PDELab GridFunction defined by a parameter class and a lambda  */
    template<typename GRIDVIEW, typename LAMBDA>
    auto makeGridFunctionFromLocalLambda (const GRIDVIEW& gridview, LAMBDA lambda)
    {
      typedef typename GRIDVIEW::template Codim<0>::Entity E;
      E e;
      typedef typename E::Geometry::LocalCoordinate X;
      X x;
      typedef decltype(lambda(e,x)) ReturnType;
      typedef typename LambdaAdapterGetRangeFieldType<ReturnType>::Type RF;
      const int dim = LambdaAdapterGetDim<ReturnType>::dim;
      typedef LocalLambdaToGridFunctionAdapter<GRIDVIEW,RF,dim,LAMBDA> TheType;
      return TheType(gridview,lambda);
    }

    /** \brief return a PDELab GridFunction defined by a parameter class and a lambda  */
    template<typename GV, typename RF, int n, typename F, typename P>
    class LocalLambdaToInstationaryGridFunctionAdapter
      : public Dune::PDELab::GridFunctionBase<Dune::PDELab::
                                              GridFunctionTraits<GV,RF,n,Dune::FieldVector<RF,n> >,
                                              LocalLambdaToInstationaryGridFunctionAdapter<GV,RF,n,F,P> >
    {
      GV gv;
      F f;
      P& p;
    public:
      typedef Dune::PDELab::
      GridFunctionTraits<GV,RF,n,Dune::FieldVector<RF,n> > Traits;

      //! construct from grid view
      LocalLambdaToInstationaryGridFunctionAdapter (const GV& gv_, F f_, P& p_) : gv(gv_), f(f_), p(p_) {}

      //! get a reference to the grid view
      inline const GV& getGridView () {return gv;}

      //! evaluate extended function on element
      inline void evaluate (const typename Traits::ElementType& e,
                            const typename Traits::DomainType& xl,
                            typename Traits::RangeType& y) const
      {
        y = f(e,xl);
        return;
      }

      // pass time to parameter object
      void setTime (RF t) {
        p.setTime(t);
      }
    };

    /** \brief return a PDELab GridFunction defined by a parameter class and a lambda  */
    template<typename GRIDVIEW, typename LAMBDA, typename PROBLEM>
    auto makeInstationaryGridFunctionFromLocalLambda (const GRIDVIEW& gridview, LAMBDA lambda, PROBLEM& problem)
    {
      typedef typename GRIDVIEW::template Codim<0>::Entity E;
      E e;
      typedef typename E::Geometry::LocalCoordinate X;
      X x;
      typedef decltype(lambda(e,x)) ReturnType;
      typedef typename LambdaAdapterGetRangeFieldType<ReturnType>::Type RF;
      const int dim = LambdaAdapterGetDim<ReturnType>::dim;
      typedef LocalLambdaToInstationaryGridFunctionAdapter<GRIDVIEW,RF,dim,LAMBDA,PROBLEM> TheType;
      return TheType(gridview,lambda,problem);
    }

    /******************************************************/
    /** \brief Adapter for globally defined boundary cond.*/
    /******************************************************/
    template<typename F>
    class GlobalLambdaToBoundaryConditionAdapter
      : public Dune::PDELab::DirichletConstraintsParameters
    {
      F f;
    public:
      //! construct from functor
      GlobalLambdaToBoundaryConditionAdapter (F f_) : f(f_) {}

      //! Test whether boundary is Dirichlet-constrained
      template<typename I>
      bool isDirichlet(const I & intersection,
                       const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord
                       ) const
      {
        Dune::FieldVector<typename I::ctype, I::dimension> xg = intersection.geometry().global(coord);
        return f(xg);
      }

      template<typename I>
      bool isNeumann(const I & ig,
                     const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord
                     ) const
      {
        return !isDirichlet( ig, coord );
      }

    };

    /** \brief get boundary condition from a lambda function */
    template<typename LAMBDA>
    auto makeBoundaryConditionFromGlobalLambda (LAMBDA lambda)
    {
      return GlobalLambdaToBoundaryConditionAdapter<LAMBDA>(lambda);
    }

    template<typename T>
    class LocalLambdaToBoundaryConditionAdapter :
      public Dune::PDELab::FluxConstraintsParameters,
      public Dune::PDELab::DirichletConstraintsParameters
    {
      const T t;

    public:

      LocalLambdaToBoundaryConditionAdapter(const T& t_ )
        : t( t_ )
      {}

      template<typename I>
      bool isDirichlet(const I & ig, const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord
                       ) const
      {
        return(t(ig.intersection(),coord));
      }

      template<typename I>
      bool isNeumann(const I & ig,
                     const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord
                     ) const
      {
        return !isDirichlet( ig, coord );
      }
    };

    /** \brief get boundary condition from a lambda function */
    template<typename LAMBDA>
    auto makeBoundaryConditionFromLocalLambda (LAMBDA lambda)
    {
      return LocalLambdaToBoundaryConditionAdapter<LAMBDA>(lambda);
    }

  }
}

#endif
