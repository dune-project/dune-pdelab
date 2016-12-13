// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_FUNCTION_CALLABLEADAPTER_HH
#define DUNE_PDELAB_FUNCTION_CALLABLEADAPTER_HH

#include <utility>
#include <type_traits>
#include <dune/common/fvector.hh>
#include <dune/common/typetraits.hh>

#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/constraints/common/constraintsparameters.hh>

namespace Dune {
  namespace PDELab {

    /************************
     * Grid function adapters
     ************************/

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
      GlobalCallableToGridFunctionAdapter (const GV& gv_, const F& f_) : gv(gv_), f(f_) {}

      //! get a reference to the grid view
      inline const GV& getGridView () const {return gv;}

      //! evaluate extended function on element
      inline void evaluate (const typename Traits::ElementType& e,
                            const typename Traits::DomainType& xl,
                            typename Traits::RangeType& y) const
      {
        typename Traits::DomainType xg = e.geometry().global(xl);
        y = f(xg);
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
      LocalCallableToGridFunctionAdapter (const GV& gv_, const F& f_) : gv(gv_), f(f_) {}

      //! get a reference to the grid view
      inline const GV& getGridView () const {return gv;}

      //! evaluate extended function on element
      inline void evaluate (const typename Traits::ElementType& e,
                            const typename Traits::DomainType& xl,
                            typename Traits::RangeType& y) const
      {
        y = f(e,xl);
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
     * \return A GlobalCallableToGridFunctionAdapter or a LocalCallableToGridFunctionAdapter.
     */
    template <typename GV, typename F>
    WrapperConformingToGridFunctionInterface makeGridFunctionFromCallable (const GV& gv, const F& f)
    {}
#endif

#ifndef DOXYGEN
    /** \brief Create PDELab GridFunction from a callable f(x) that expects a global coordinate x */
    template <typename GV, typename F>
    auto makeGridFunctionFromCallable (const GV& gv, const F& f)
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

    /** \brief Create PDELab GridFunction from a callable f(e,x) that expects
        an entity e and a local coordinate x */
    template <typename GV, typename F>
    auto makeGridFunctionFromCallable (const GV& gv, const F& f)
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


    /*************************************
     * Instationary grid function adapters
     *************************************/


    /** \brief return a PDELab GridFunction (with setTime method) defined by a parameter class and a callable f(x)
        global coordinates x */
    template<typename GV, typename RF, int n, typename F, typename P>
    class GlobalCallableToInstationaryGridFunctionAdapter
      : public Dune::PDELab::GridFunctionBase<Dune::PDELab::
                                              GridFunctionTraits<GV,RF,n,Dune::FieldVector<RF,n> >,
                                              GlobalCallableToInstationaryGridFunctionAdapter<GV,RF,n,F,P> >
    {
      GV gv;
      F f;
      P& p;
    public:
      typedef Dune::PDELab::
      GridFunctionTraits<GV,RF,n,Dune::FieldVector<RF,n> > Traits;

      //! construct from grid view
      GlobalCallableToInstationaryGridFunctionAdapter (const GV& gv_, const F& f_, P& p_)
        : gv(gv_), f(f_), p(p_)
      {}

      //! get a reference to the grid view
      inline const GV& getGridView () const {return gv;}

      //! evaluate extended function on element
      inline void evaluate (const typename Traits::ElementType& e,
                            const typename Traits::DomainType& xl,
                            typename Traits::RangeType& y) const
      {
        typename Traits::DomainType xg = e.geometry().global(xl);
        y = f(xg);
      }

      // pass time to parameter object
      void setTime (RF t) {
        p.setTime(t);
      }
    };

    /** \brief return a PDELab GridFunction (with setTime method) defined by a parameter class and a callable f(e,x)
        that expects an entity e and local coordinates x */
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
      LocalCallableToInstationaryGridFunctionAdapter (const GV& gv_, const F& f_, P& p_) : gv(gv_), f(f_), p(p_) {}

      //! get a reference to the grid view
      inline const GV& getGridView () const {return gv;}

      //! evaluate extended function on element
      inline void evaluate (const typename Traits::ElementType& e,
                            const typename Traits::DomainType& xl,
                            typename Traits::RangeType& y) const
      {
        y = f(e,xl);
      }

      // pass time to parameter object
      void setTime (RF t) {
        p.setTime(t);
      }
    };

#ifdef DOXYGEN
    //! \brief Create a GridFunction from callable and parameter class with setTime method
    /**
     * \param gv A GridView.
     * \param f A callable of one of the two forms:
     *          1. f(x) taking a global coordinate x of type
     *          typename GV::template Codim<0>::Entity::Geometry::GlobalCoordinate.
     *          2. f(e,x) taking an Entity e
     *          coordinate x  or of the form f(e,x) taking an Entity e and a local
     *          coordinate x of type Entity::Geometry::LocalCoordinate.
     * \param parameter Parameter class.
     * \return A GlobalCallableToInstationaryGridFunctionAdapter or a
     *         LocalCallableToInstationaryGridFunctionAdapter
     */
    template <typename GV, typename F>
    WrapperConformingToGridFunctionInterface makeInstationaryGridFunctionFromCallable (const GV& gv, const F& f)
    {}
#endif

#ifndef DOXYGEN
    /** \brief Create PDELab GridFunction with setTime method from a callable f(e,x)
        that expects an global coordinate x */
    template <typename GV, typename F, typename PARAM>
    auto makeInstationaryGridFunctionFromCallable (const GV& gv, const F& f, PARAM& param)
      -> typename std::enable_if<
        AlwaysTrue <
          decltype(f(std::declval<typename GV::template Codim<0>::Entity::Geometry::GlobalCoordinate>()))
          >::value,
        GlobalCallableToInstationaryGridFunctionAdapter<
          GV,
          typename CallableAdapterGetRangeFieldType<
            decltype(f(std::declval<typename GV::template Codim<0>::Entity::Geometry::GlobalCoordinate>()))
            >::Type,
          CallableAdapterGetDim<
            decltype(f(std::declval<typename GV::template Codim<0>::Entity::Geometry::GlobalCoordinate>()))
            >::dim,
          F,
          PARAM>
        >::type
    {
      typedef typename GV::template Codim<0>::Entity::Geometry::GlobalCoordinate X;
      X x;
      typedef decltype(f(x)) ReturnType;
      typedef typename CallableAdapterGetRangeFieldType<ReturnType>::Type RF;
      const int dim = CallableAdapterGetDim<ReturnType>::dim;
      typedef GlobalCallableToInstationaryGridFunctionAdapter<GV,RF,dim,F,PARAM> TheType;
      return TheType(gv,f,param);
    }

    /** \brief Create PDELab GridFunction with setTime method from a callable f(e,x)
        that expects an entity e and a local coordinate x */
    template <typename GV, typename F, typename PARAM>
    auto makeInstationaryGridFunctionFromCallable (const GV& gv, const F& f, PARAM& param)
      -> typename std::enable_if<
        AlwaysTrue <
          decltype(f(
                     std::declval<typename GV::template Codim<0>::Entity>(),
                     std::declval<typename GV::template Codim<0>::Entity::Geometry::LocalCoordinate>()
                     ))
          >::value,
        LocalCallableToInstationaryGridFunctionAdapter<
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
          F,
          PARAM>
        >::type
    {
      typedef typename GV::template Codim<0>::Entity E;
      E e;
      typedef typename E::Geometry::LocalCoordinate X;
      X x;
      typedef decltype(f(e,x)) ReturnType;
      typedef typename CallableAdapterGetRangeFieldType<ReturnType>::Type RF;
      const int dim = CallableAdapterGetDim<ReturnType>::dim;
      typedef LocalCallableToInstationaryGridFunctionAdapter<GV,RF,dim,F,PARAM> TheType;
      return TheType(gv,f,param);
    }
#endif // DOXYGEN


    /*****************************
     * Boundary condition adapters
     *****************************/

    /** \brief Adapter for boundary cond from a callable taking global coordinates*/
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

    /** \brief Adapter for boundary cond from a callable taking an entity and local coordinates*/
    template<typename F>
    class LocalCallableToBoundaryConditionAdapter :
      public Dune::PDELab::FluxConstraintsParameters,
      public Dune::PDELab::DirichletConstraintsParameters
    {
      const F f;

    public:

      LocalCallableToBoundaryConditionAdapter(const F& f_ )
        : f( f_ )
      {}

      template<typename I>
      bool isDirichlet(const I & ig, const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord
                       ) const
      {
        return(f(ig.intersection(),coord));
      }

      template<typename I>
      bool isNeumann(const I & ig,
                     const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord
                     ) const
      {
        return !isDirichlet( ig, coord );
      }
    };

#ifdef DOXYGEN
    //! \brief Create a BoundaryConditionAdapter from a callable
    /**
     * \param gv A GridView
     * \param f A callable of one of the two forms:
     *          1. f(x) taking a global coordinate x of type
     *          typename GV::template Codim<0>::Entity::Geometry::GlobalCoordinate.
     *          2. f(e,x) taking an Entity e
     *          coordinate x  or of the form f(e,x) taking an Entity e and a local
     *          coordinate x of type Entity::Geometry::LocalCoordinate.
     * \return A GlobalCallableToBoundaryConditionAdapter or a
     *         LocalCallableToBoundaryConditionAdapter.
     */
    template <typename GV, typename F>
    BoundaryConditionAdapter makebBoundaryConditionFromCallable (const GV& gv, const F& f)
#endif

#ifndef DOXYGEN
    /** \brief Create BoundaryConditionAdapter from a callable f(x) that expects a global coordinate x */
    template<typename GV, typename F>
    auto makeBoundaryConditionFromCallable (const GV& gv, const F& f)
      -> typename std::enable_if<
        AlwaysTrue <
          decltype(f(std::declval<typename GV::template Codim<0>::Entity::Geometry::GlobalCoordinate>()))
          >::value,
        GlobalCallableToBoundaryConditionAdapter<F>
        >::type
    {
      return GlobalCallableToBoundaryConditionAdapter<F>(f);
    }

    /** \brief Create BoundaryConditionAdapter from a callable f(e,x) that expects
        an entity e and a global coordinate x */
    template<typename GV, typename F>
    auto makeBoundaryConditionFromCallable (const GV& gv, const F& f)
      -> typename std::enable_if<
        AlwaysTrue <
          decltype(f(
                     std::declval<typename GV::template Codim<0>::Entity>(),
                     std::declval<typename GV::template Codim<0>::Entity::Geometry::LocalCoordinate>()
                     ))
      >::value,
        LocalCallableToBoundaryConditionAdapter<F>
        >::type
    {
      return LocalCallableToBoundaryConditionAdapter<F>(f);
    }
#endif


  }
}
#endif // DUNE_PDELAB_FUNCTION_CALLABLEADAPTER_HH
