// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_FUNCTION_CALLABLEADAPTER_HH
#define DUNE_PDELAB_FUNCTION_CALLABLEADAPTER_HH

#include<utility>

namespace Dune {
  namespace PDELab {

    /** \brief Adapter for callables f(x) expecting a global coordinate x */
    template<typename GV, typename RF, int n, typename F>
    class GlobalCallableToGridFunctionAdapter
      : public Dune::PDELab::GridFunctionBase<Dune::PDELab::
                                              GridFunctionTraits<GV,RF,n,Dune::FieldVector<RF,n> >,
                                              GlobalCallableToGridFunctionAdapter<GV,RF,n,F> >
    {
      GV gv;
      F f;
    public:
      typedef Dune::PDELab::
      GridFunctionTraits<GV,RF,n,Dune::FieldVector<RF,n> > Traits;

      //! construct from grid view
      GlobalCallableToGridFunctionAdapter (const GV& gv_, F f_) : gv(gv_), f(f_) {}

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
    struct CallableAdapterGetDim {
      enum {dim=1};
    };

    template<typename T, int n>
    struct CallableAdapterGetDim< FieldVector<T,n> > {
      enum {dim=n};
    };

    template<typename T>
    struct CallableAdapterGetRangeFieldType {
      typedef T Type;
    };

    template<typename T, int n>
    struct CallableAdapterGetRangeFieldType< FieldVector<T,n> > {
      typedef T Type;
    };


    /** \brief Adapter for callables f(e,x) expecting an entity e and a global coordinate x */
    template<typename GV, typename RF, int n, typename F>
    class LocalCallableToGridFunctionAdapter
      : public Dune::PDELab::GridFunctionBase<Dune::PDELab::
                                              GridFunctionTraits<GV,RF,n,Dune::FieldVector<RF,n> >,
                                              LocalCallableToGridFunctionAdapter<GV,RF,n,F> >
    {
      GV gv;
      F f;
    public:
      typedef Dune::PDELab::
      GridFunctionTraits<GV,RF,n,Dune::FieldVector<RF,n> > Traits;

      //! construct from grid view
      LocalCallableToGridFunctionAdapter (const GV& gv_, F f_) : gv(gv_), f(f_) {}

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

#ifdef DOXYGEN
    //! \brief Create a GridFunction adapter from a callable
    /**
     * \param gv A GridView
     * \param f A callable of one of the two forms:
     *          1. f(x) taking a global coordinate x of type
     *          typename GV::template Codim<0>::Entity::Geometry::GlobalCoordinate.
     *          2. f(e,x) taking an Entity e
     *          coordinate x  or of the form f(e,x) taking an Entity e and a local
     *          coordinate x of type Entity::Geometry::LocalCoordinate.
     */
    WrapperConformingToGridFunctionInterface makeGridFunctionFromCallable (const GV& gv, F f)
    {}
#endif

#ifndef DOXYGEN
    /** \brief Create PDELab GridFunction from a callable f(x) that expects a global coordinate x */
    template <typename GV, typename F>
    auto makeGridFunctionFromCallable (const GV& gv, F f)
      -> typename std::enable_if<
        AlwaysTrue <
          decltype(f(std::declval<typename GV::template Codim<0>::Entity::Geometry::GlobalCoordinate>()))
          >::value,
        GlobalCallableToGridFunctionAdapter<
          GV,
          typename CallableAdapterGetRangeFieldType<
            decltype(f(std::declval<typename GV::template Codim<0>::Entity::Geometry::GlobalCoordinate>()))
            >::Type,
          CallableAdapterGetDim<
            decltype(f(std::declval<typename GV::template Codim<0>::Entity::Geometry::GlobalCoordinate>()))
            >::dim,
          F>
        >::type
    {
      typedef typename GV::template Codim<0>::Entity::Geometry::GlobalCoordinate X;
      X x;
      typedef decltype(f(x)) ReturnType;
      typedef typename CallableAdapterGetRangeFieldType<ReturnType>::Type RF;
      const int dim = CallableAdapterGetDim<ReturnType>::dim;
      typedef GlobalCallableToGridFunctionAdapter<GV,RF,dim,F> TheType;
      return TheType(gv,f);
    }

    /** \brief Create PDELab GridFunction from a callable f(e,x) that expects an entity e and a local coordinate x */
    template <typename GV, typename F>
    auto makeGridFunctionFromCallable (const GV& gv, F f)
      -> typename std::enable_if<
        AlwaysTrue <
          decltype(f(
                     std::declval<typename GV::template Codim<0>::Entity>(),
                     std::declval<typename GV::template Codim<0>::Entity::Geometry::LocalCoordinate>()
                     ))
          >::value,
        LocalCallableToGridFunctionAdapter<
          GV,
          typename CallableAdapterGetRangeFieldType<
            decltype(f(
                       std::declval<typename GV::template Codim<0>::Entity>(),
                       std::declval<typename GV::template Codim<0>::Entity::Geometry::LocalCoordinate>()
                       ))
            >::Type,
          CallableAdapterGetDim<
            decltype(f(
                       std::declval<typename GV::template Codim<0>::Entity>(),
                       std::declval<typename GV::template Codim<0>::Entity::Geometry::LocalCoordinate>()
                       ))
            >::dim,
          F>
        >::type
    {
      typedef typename GV::template Codim<0>::Entity E;
      E e;
      typedef typename E::Geometry::LocalCoordinate X;
      X x;
      typedef decltype(f(e,x)) ReturnType;
      typedef typename CallableAdapterGetRangeFieldType<ReturnType>::Type RF;
      const int dim = CallableAdapterGetDim<ReturnType>::dim;
      typedef LocalCallableToGridFunctionAdapter<GV,RF,dim,F> TheType;
      return TheType(gv,f);
    }
#endif // DOXYGEN

    /** \brief return a PDELab GridFunction defined by a parameter class and a lambda  */
    template<typename GV, typename RF, int n, typename F, typename P>
    class LocalCallableToInstationaryGridFunctionAdapter
      : public Dune::PDELab::GridFunctionBase<Dune::PDELab::
                                              GridFunctionTraits<GV,RF,n,Dune::FieldVector<RF,n> >,
                                              LocalCallableToInstationaryGridFunctionAdapter<GV,RF,n,F,P> >
    {
      GV gv;
      F f;
      P& p;
    public:
      typedef Dune::PDELab::
      GridFunctionTraits<GV,RF,n,Dune::FieldVector<RF,n> > Traits;

      //! construct from grid view
      LocalCallableToInstationaryGridFunctionAdapter (const GV& gv_, F f_, P& p_) : gv(gv_), f(f_), p(p_) {}

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
    auto makeInstationaryGridFunctionFromLocalCallable (const GRIDVIEW& gridview, LAMBDA lambda, PROBLEM& problem)
    {
      typedef typename GRIDVIEW::template Codim<0>::Entity E;
      E e;
      typedef typename E::Geometry::LocalCoordinate X;
      X x;
      typedef decltype(lambda(e,x)) ReturnType;
      typedef typename CallableAdapterGetRangeFieldType<ReturnType>::Type RF;
      const int dim = CallableAdapterGetDim<ReturnType>::dim;
      typedef LocalCallableToInstationaryGridFunctionAdapter<GRIDVIEW,RF,dim,LAMBDA,PROBLEM> TheType;
      return TheType(gridview,lambda,problem);
    }

    /******************************************************/
    /** \brief Adapter for globally defined boundary cond.*/
    /******************************************************/
    template<typename F>
    class GlobalCallableToBoundaryConditionAdapter
      : public Dune::PDELab::DirichletConstraintsParameters
    {
      F f;
    public:
      //! construct from functor
      GlobalCallableToBoundaryConditionAdapter (F f_) : f(f_) {}

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
    auto makeBoundaryConditionFromGlobalCallable (LAMBDA lambda)
    {
      return GlobalCallableToBoundaryConditionAdapter<LAMBDA>(lambda);
    }

    template<typename T>
    class LocalCallableToBoundaryConditionAdapter :
      public Dune::PDELab::FluxConstraintsParameters,
      public Dune::PDELab::DirichletConstraintsParameters
    {
      const T t;

    public:

      LocalCallableToBoundaryConditionAdapter(const T& t_ )
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
    auto makeBoundaryConditionFromLocalCallable (LAMBDA lambda)
    {
      return LocalCallableToBoundaryConditionAdapter<LAMBDA>(lambda);
    }

  }
}

#endif
