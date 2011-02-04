//-*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_FUNCTION_HH
#define DUNE_PDELAB_FUNCTION_HH

#include <iostream>
#include <sstream>

#include <dune/common/deprecated.hh>
#include <dune/common/static_assert.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dune/grid/utility/hierarchicsearch.hh>

#include "typetree.hh"
#include "vtkexport.hh"
#include "geometrywrapper.hh"
#include "multitypetree.hh"

namespace Dune {
  namespace PDELab {

    //! \addtogroup PDELab_Function Function
    //! \ingroup PDELab
    //! \{

	//! traits class holding function signature, same as in local function
	template<class DF, int n, class D, class RF, int m, class R>
	struct FunctionTraits
	{
	  //! \brief Export type for domain field
	  typedef DF DomainFieldType;

	  //! \brief Enum for domain dimension
	  enum {
		//! \brief dimension of the domain
		dimDomain = n
	  };

	  //! \brief domain type
	  typedef D DomainType;

	  //! \brief Export type for range field
	  typedef RF RangeFieldType;

	  //! \brief Enum for range dimension
	  enum {
		//! \brief dimension of the range
		dimRange = m
	  };

	  //! \brief range type
	  typedef R RangeType;
	};


	//! a Function maps x in DomainType to y in RangeType
	template<class T, class Imp>
	class FunctionInterface
	{
	public:
	  //! \brief Export type traits
	  typedef T Traits;

	  /** \brief Evaluate all basis function at given position

		  Evaluates all shape functions at the given position and returns
		  these values in a vector.
	  */
	  inline void evaluate (const typename Traits::DomainType& x,
							typename Traits::RangeType& y) const
	  {
		asImp().evaluate(x,y);
	  }

	private:
	  Imp& asImp () {return static_cast<Imp &> (*this);}
	  const Imp& asImp () const {return static_cast<const Imp &>(*this);}
	};

    //! Default class for additional methods in instationary functions
    class InstationaryFunctionDefaults
    {
    public:
      //! set time for subsequent evaluation
      /**
       * This method set the time for subsequent calls to any of the
       * evaluation methods.
       *
       * \note This default method does nothing, it just ensures setTime() can
       *       be called without ill effects.
       * \note Function implementation are free to restrict the types of
       *       acceptable parameters.  This should be noted in the function
       *       classes documentation.
       */
      template<typename Time>
      inline void setTime(Time t)
      { }
    };

    template<typename GV>
    struct PowerCompositeGridFunctionTraits
    {
      typedef GV GridViewType;

	  //! \brief codim 0 entity
	  typedef typename GV::Traits::template Codim<0>::Entity ElementType;

    };

	//! traits class holding function signature, same as in local function
	template<class GV, class RF, int m, class R>
	struct GridFunctionTraits
	  : public FunctionTraits<typename GV::Grid::ctype, GV::dimension,
				    		  Dune::FieldVector<typename GV::Grid::ctype,
                                                GV::dimension>,
							  RF, m, R>
      , public PowerCompositeGridFunctionTraits<GV>
	{
	};

	//! a GridFunction maps x in DomainType to y in RangeType
	template<class T, class Imp>
	class GridFunctionInterface
	{
	public:
	  //! \brief Export type traits
	  typedef T Traits;

	  /** \brief Evaluate the GridFunction at given position

		  Evaluates components of the grid function at the given position and
		  returns these values in a vector.

          \param[in]  e The entity to evaluate on
          \param[in]  x The position in entity-local coordinates
          \param[out] y The result of the evaluation
	  */
	  inline void evaluate (const typename Traits::ElementType& e,
							const typename Traits::DomainType& x,
							typename Traits::RangeType& y) const
	  {
		asImp().evaluate(e,x,y);
	  }

      //! get a reference to the GridView
      /* \note This is deprecated in favor of "getGridView() const" */
      DUNE_DEPRECATED
	  inline const typename Traits::GridViewType& getGridView ()
	  {
		return asImp().getGridView();
	  }

      //! get a reference to the GridView
	  inline const typename Traits::GridViewType& getGridView () const
	  {
		return asImp().getGridView();
	  }

	private:
	  Imp& asImp () {return static_cast<Imp &> (*this);}
	  const Imp& asImp () const {return static_cast<const Imp &>(*this);}
	};

	//! traits class holding function signature, same as in local function
	template<class GV, class RF, int m, class R>
	struct BoundaryGridFunctionTraits
	  : public FunctionTraits<typename GV::Grid::ctype, GV::dimension-1,
				    		  Dune::FieldVector<typename GV::Grid::ctype,
                                                GV::dimension-1>,
							  RF, m, R>
	{
	  //! \brief Export grid view type in addition
	  typedef GV GridViewType;
	};


	//! a BoundaryGridFunction allows evaluation on boundary intersections
    // T are BoundaryGridFunctionTraits
	template<class T, class Imp>
	class BoundaryGridFunctionInterface
	{
	public:
	  //! \brief Export type traits
	  typedef T Traits;

	  /** \brief Evaluate the GridFunction at given position

		  Evaluates components of the grid function at the given position and
		  returns these values in a vector.

          \param[in]  ig geometry of intersection with boundary
          \param[in]  x The position in entity-local coordinates
          \param[out] y The result of the evaluation
	  */
      template<typename I>
	  inline void evaluate (const IntersectionGeometry<I>& ig,
							const typename Traits::DomainType& x,
							typename Traits::RangeType& y) const
	  {
		asImp().evaluate(ig,x,y);
	  }

      //! get a reference to the GridView
	  inline const typename Traits::GridViewType& getGridView () const
	  {
		return asImp().getGridView();
	  }

	private:
	  Imp& asImp () {return static_cast<Imp &> (*this);}
	  const Imp& asImp () const {return static_cast<const Imp &>(*this);}
	};

    //! \addtogroup PDELab_FunctionAdapters Function Adapters
    //! \{

	/** \brief make a GridFunction from a Function
     *
     *  \tparam G The GridView type
     *  \tparam T The function type
     */
	template<typename G, typename T>
	class FunctionToGridFunctionAdapter :
	  public GridFunctionInterface<GridFunctionTraits<
									 G,
									 typename T::Traits::RangeFieldType,
									 T::Traits::dimRange,
									 typename T::Traits::RangeType>,
								   FunctionToGridFunctionAdapter<G,T> >
	{
	public:
	  typedef GridFunctionTraits<G,
								 typename T::Traits::RangeFieldType,
								 T::Traits::dimRange,
								 typename T::Traits::RangeType> Traits;
      dune_static_assert(
       (is_same<typename T::Traits::DomainFieldType,
                typename Traits::DomainFieldType>::value),
       "GridView's and wrapped Functions DomainFieldType don't match");
      dune_static_assert(
       T::Traits::dimDomain==Traits::dimDomain,
       "GridView's and wrapped Functions dimDomain don't match");
      dune_static_assert(
       (is_same<typename T::Traits::DomainType,
                typename Traits::DomainType>::value),
       "GridView's and wrapped Functions DomainType don't match");

      /** \brief Create a FunctionToGridFunctionAdapter
       *
       *  \param g_ The GridView
       *  \param t_ The function
       */
	  FunctionToGridFunctionAdapter (const G& g_, const T& t_) : g(g_), t(t_) {}

	  inline void evaluate (const typename Traits::ElementType& e,
							const typename Traits::DomainType& x,
							typename Traits::RangeType& y) const
	  {
		t.evaluate(e.geometry().global(x),y);
	  }

	  inline const typename Traits::GridViewType& getGridView () const
	  {
		return g;
	  }

	private:
	  const G& g;
	  const T& t;
	};

	/** \brief make a Function from a GridFunction
     *
     *  \tparam GF The GridFunction type
     */
	template<typename GF>
	class GridFunctionToFunctionAdapter
      : public FunctionInterface<typename GF::Traits, GridFunctionToFunctionAdapter<GF> >
	{
	public:
	  //! \brief Export type traits
	  typedef typename GF::Traits Traits;

      //! make a GridFunctionToFunctionAdapter
      GridFunctionToFunctionAdapter(const GF &gf_)
        : gf(gf_)
        , hsearch(gf.getGridView().grid(), gf.getGridView().indexSet())
      { }

	  /** \brief Evaluate all basis function at given position

		  Evaluates all shape functions at the given position and returns
		  these values in a vector.
	  */
	  inline void evaluate (const typename Traits::DomainType& x,
							typename Traits::RangeType& y) const
	  {
        typename Traits::GridViewType::Grid::Traits::template Codim<0>::EntityPointer
          ep = hsearch.findEntity(x);
        gf.evaluate(*ep, ep->geometry().local(x), y);
	  }

	private:
      const GF &gf;
      const Dune::HierarchicSearch<typename Traits::GridViewType::Grid,
                                   typename Traits::GridViewType::IndexSet> hsearch;
	};


	/** \brief make a Function in local coordinates from a Function in global coordinates
     *
     *  \tparam T Type of the global function
     *  \tparam E Type of the grid's element
     */
	template<typename T, typename E>
	class GlobalFunctionToLocalFunctionAdapter :
	  public FunctionInterface<typename T::Traits,
							   GlobalFunctionToLocalFunctionAdapter<T,E> >
	{
	public:
	  typedef typename T::Traits Traits;

      /** \brief Create a GlobalFunctionToLocalFunctionAdapter
       *
       *  \param t_ Global function
       *  \param e_ Grid's element where the local function is defined
       */
	  GlobalFunctionToLocalFunctionAdapter (const T& t_, const E& e_) : t(t_), e(e_) {}

	  /** \brief Evaluate the local function at the given position

          \param[in]  x The position in local coordinates
          \param[out] y The result of the evaluation
	  */
	  inline void evaluate (const typename Traits::DomainType& x,
							typename Traits::RangeType& y) const
	  {
		t.evaluate(e.geometry().global(x),y);
	  }

	private:
	  const T& t;
	  const E& e;
	};


    /** \brief make a LocalFunction from a GridFunction using local coordinates
     *
     *  \tparam T type of the GridFunction
     */
	template<typename T> // T: GridFunction, E: Entity
	class GridFunctionToLocalFunctionAdapter :
	  public FunctionInterface<typename T::Traits,
							   GridFunctionToLocalFunctionAdapter<T> >
	{
	public:
	  typedef typename T::Traits Traits;

      /** \brief Create a GridFunctionToLocalFunctionAdapter
       *
       *  \param t_ GridFunction
       *  \param e_ Grid's element where the local function is defined
       */
	  GridFunctionToLocalFunctionAdapter (const T& t_,
										  const typename Traits::ElementType& e_)
		: t(t_), e(e_) {}

	  /** \brief Evaluate the local function at the given position

          \param[in]  x The position in local coordinates
          \param[out] y The result of the evaluation
	  */
	  inline void evaluate (const typename Traits::DomainType& x,
							typename Traits::RangeType& y) const
	  {
		t.evaluate(e,x,y);
	  }

	private:
	  const T& t;
	  const typename Traits::ElementType& e;
	};


	//! a Function maps x in DomainType to y in RangeType
	template<class T>
	class SelectComponentAdapter : public FunctionInterface<FunctionTraits<typename T::Traits::DomainFieldType,T::Traits::dimDomain,typename T::Traits::DomainType,typename T::Traits::RangeFieldType,1,Dune::FieldVector<typename T::Traits::RangeFieldType,1> > , SelectComponentAdapter<T> >
	{
      typedef FunctionInterface<FunctionTraits<typename T::Traits::DomainFieldType,T::Traits::dimDomain,typename T::Traits::DomainType,typename T::Traits::RangeFieldType,1,Dune::FieldVector<typename T::Traits::RangeFieldType,1> > , SelectComponentAdapter<T> > BaseT;
	public:
	  //! \brief Export type traits
	  typedef typename BaseT::Traits Traits;

      SelectComponentAdapter (const T& t_, int k_) : t(t_), k(k_) {}

	  /** \brief Evaluate all basis function at given position

		  Evaluates all shape functions at the given position and returns
		  these values in a vector.
	  */
	  inline void evaluate (const typename Traits::DomainType& x,
							typename Traits::RangeType& y) const
	  {
        typename T::Traits::RangeType Y;
        t.evaluate(x,Y);
        y = Y[k];
	  }

      //! set component to be selected
      void select (int k_)
      {
        k = k_;
      }

	private:
      const T& t;
      int k;
	};

	//! Takes a BoundaryGridFunction and acts as a single component
	template<class T>
	class BoundaryGridFunctionSelectComponentAdapter
      : public BoundaryGridFunctionInterface<BoundaryGridFunctionTraits<typename T::Traits::GridViewType,
                                                                        typename T::Traits::RangeFieldType,1,
                                                                        Dune::FieldVector<typename T::Traits::RangeFieldType,1> > ,
                                             BoundaryGridFunctionSelectComponentAdapter<T> >
	{
      typedef BoundaryGridFunctionInterface<BoundaryGridFunctionTraits<typename T::Traits::GridViewType,
                                                                       typename T::Traits::RangeFieldType,1,
                                                                       Dune::FieldVector<typename T::Traits::RangeFieldType,1> > ,
                                            BoundaryGridFunctionSelectComponentAdapter<T> > BaseT;
    public:
	  //! \brief Export type traits
	  typedef typename BaseT::Traits Traits;

      BoundaryGridFunctionSelectComponentAdapter (const T& t_, int k_) : t(t_), k(k_) {}

	  /** \brief Evaluate all basis function at given position

		  Evaluates all shape functions at the given position and returns
		  these values in a vector.
	  */
      template<typename I>
	  inline void evaluate (const IntersectionGeometry<I>& ig,
							const typename Traits::DomainType& x,
							typename Traits::RangeType& y) const
	  {
        typename T::Traits::RangeType Y;
        t.evaluate(ig,x,Y);
        y = Y[k];
	  }

      //! get a reference to the GridView
	  inline const typename Traits::GridViewType& getGridView () const
	  {
		return t.getGridView();
	  }


      //! set component to be selected
      void select (int k_)
      {
        k = k_;
      }

	private:
      const T& t;
      int k;
	};

    //! \}

	//============================
	// Function tree
	//============================

    //! \addtogroup GridFunctionTree
    //! \{

	/** \brief leaf of a function tree
     *
     *  Classes derived from this class implement a \ref GridFunctionTree.
     *
     *  \tparam T   Traits class holding the functions signature
     *  \tparam Imp Class implementing the function.  Imp must be derived from
     *              GridFunctionBase in some way (Barton-Nackman-Trick).
     */
	template<class T, class Imp>
	class GridFunctionBase
      : public GridFunctionInterface<T,Imp>
      , public TypeTree::LeafNode
	{
	public:
      //! Type of the GridView
	  typedef typename T::GridViewType GridViewType;
	};


	/** \brief leaf of a function tree
     *
     *  Classes derived from this class implement a \ref GridFunctionTree.
     *
     *  \tparam T   Traits class holding the functions signature
     *  \tparam Imp Class implementing the function.  Imp must be derived from
     *              GridFunctionBase in some way (Barton-Nackman-Trick).
     */
	template<class T, class Imp>
	class BoundaryGridFunctionBase
      : public BoundaryGridFunctionInterface<T,Imp>
      , public TypeTree::LeafNode
	{
	public:
      //! Type of the GridView
	  typedef typename T::GridViewType GridViewType;
	};

    struct PowerGridFunctionTag {};

	/** \brief product of identical functions
     *
     *  This collects k instances of T in a \ref GridFunctionTree.
     *
     *  \tparam T The type of the children of this node in the tree.
     *  \tparam k The number of children this node has.
     */
	template<class T, std::size_t k>
	class PowerGridFunction
      : public TypeTree::PowerNode<T,k>
	{

      typedef TypeTree::PowerNode<T,k> BaseT;

	public:

      typedef PowerCompositeGridFunctionTraits<typename T::GridViewType> Traits;

      typedef PowerGridFunctionTag ImplementationTag;

      //! record the GridView
	  typedef typename T::GridViewType GridViewType;

      //! Construct a PowerGridFunction with k clones of the function t
	  PowerGridFunction (T& t)
        : BaseT(t) {}

      /** \brief Initialize all children with different function objects
       *
       *  This constructor is only available in the non-specialized version
       *
       *  \param t Points to an array of pointers to function objects of type
       *           T.  The function pointed to by the first pointer will be
       *           used to initialize the first child, the second pointer for
       *           the second child and so on.
       */
	  // TODO: PowerGridFunction (T** t) : ...

#ifdef DOXYGEN
      /** \brief Initialize all children with different function objects
       *
       *  Currently there exist specializations for 2 <= k <= 9.  Each
       *  specialization has a constructor which takes the initializers for
       *  its children as arguments.
       *
       *  @param t0 The initializer for the first child.
       *  @param t1 The initializer for the second child.
       *  @param ... more initializers
       */
      PowerGridFunction (T& t0, T& t1, ...)
      {
      }

#else

      template<std::size_t K = 2>
      PowerGridFunction (typename enable_if<K == BaseT::CHILDREN,T>::type& c0,
                         T& c1)
        : BaseT(c0,c1)
      {
      }

      template<std::size_t K = 3>
      PowerGridFunction (typename enable_if<K == BaseT::CHILDREN,T>::type& c0,
                         T& c1,
                         T& c2)
        : BaseT(c0,c1,c2)
      {
      }

      template<std::size_t K = 4>
      PowerGridFunction (typename enable_if<K == BaseT::CHILDREN,T>::type& c0,
                         T& c1,
                         T& c2,
                         T& c3)
        : BaseT(c0,c1,c2,c3)
      {
      }

      template<std::size_t K = 5>
      PowerGridFunction (typename enable_if<K == BaseT::CHILDREN,T>::type& c0,
                         T& c1,
                         T& c2,
                         T& c3,
                         T& c4)
        : BaseT(c0,c1,c2,c3,c4)
      {
      }

      template<std::size_t K = 6>
      PowerGridFunction (typename enable_if<K == BaseT::CHILDREN,T>::type& c0,
                         T& c1,
                         T& c2,
                         T& c3,
                         T& c4,
                         T& c5)
        : BaseT(c0,c1,c2,c3,c4,c5)
      {
      }

      template<std::size_t K = 7>
      PowerGridFunction (typename enable_if<K == BaseT::CHILDREN,T>::type& c0,
                         T& c1,
                         T& c2,
                         T& c3,
                         T& c4,
                         T& c5,
                         T& c6)
        : BaseT(c0,c1,c2,c3,c4,c5,c6)
      {
      }

      template<std::size_t K = 8>
      PowerGridFunction (typename enable_if<K == BaseT::CHILDREN,T>::type& c0,
                         T& c1,
                         T& c2,
                         T& c3,
                         T& c4,
                         T& c5,
                         T& c6,
                         T& c7)
        : BaseT(c0,c1,c2,c3,c4,c5,c6,c7)
      {
      }

      template<std::size_t K = 9>
      PowerGridFunction (typename enable_if<K == BaseT::CHILDREN,T>::type& c0,
                         T& c1,
                         T& c2,
                         T& c3,
                         T& c4,
                         T& c5,
                         T& c6,
                         T& c7,
                         T& c8)
        : BaseT(c0,c1,c2,c3,c4,c5,c6,c7,c8)
      {
      }

      template<std::size_t K = 10>
      PowerGridFunction (typename enable_if<K == BaseT::CHILDREN,T>::type& c0,
                         T& c1,
                         T& c2,
                         T& c3,
                         T& c4,
                         T& c5,
                         T& c6,
                         T& c7,
                         T& c8,
                         T& c9)
        : BaseT(c0,c1,c2,c3,c4,c5,c6,c7,c8,c9)
      {
      }

#endif // DOXYGEN
	};

    struct CompositeGridFunctionTag {};

    /** \brief composite functions
     *
     *  Collect instances of possibly different function types Tn within a
     *  \ref GridFunctionTree.  This impolements a \ref GridFunctionTree
     *
     *  \tparam Tn The base types.  Tn==EmptyChild means that slot n is
     *             unused.  Currently, up to 9 slots are supported, making 8
     *             the maximum n.
     */
	template<DUNE_TYPETREE_COMPOSITENODE_TEMPLATE_CHILDREN>
	class CompositeGridFunction
	  : public DUNE_TYPETREE_COMPOSITENODE_BASETYPE
	{

      typedef DUNE_TYPETREE_COMPOSITENODE_BASETYPE BaseT;

	public:

      typedef CompositeGridFunctionTag ImplementationTag;

      typedef PowerCompositeGridFunctionTraits<typename BaseT::template Child<0>::Type::GridViewType> Traits;

      //! record the GridView
	  typedef typename BaseT::template Child<0>::Type::GridViewType GridViewType;

	  CompositeGridFunction (DUNE_TYPETREE_COMPOSITENODE_CONSTRUCTOR_SIGNATURE)
		: BaseT(DUNE_TYPETREE_COMPOSITENODE_CHILDVARIABLES_THROUGH_FUNCTION(TypeTree::assertGridViewType<typename BaseT::template Child<0>::Type>))
	  {
	  }

#ifdef DOXYGEN
      /** \brief Initialize all children
       *
       *  @param t0 The initializer for the first child.
       *  @param t1 The initializer for the second child.
       *  @param ... more initializers
       *
       *  The actual number of arguments for this constructor corresponds to
       *  the number of slots used in the template parameter list of the class.
       */
	  CompositeGridFunction (T0& t0, T1& t1, ...) {}
#endif //DOXYGEN
	};

    //========================================================
    // helper template to turn an ordinary GridFunction into a
    // GridFunctionTree leaf
    //========================================================
    //! Turn an ordinary GridFunction into a GridFunctionTree leaf
    /**
     *  \tparam Imp Class implementing the function.
     */
    template<class Imp>
    class GridFunctionBaseAdapter
      : public GridFunctionBase<typename Imp::Traits,
                                GridFunctionBaseAdapter<Imp> >
    {
      const Imp &imp;

    public:
      //! construct a GridFunctionBaseAdapter
      /**
       * \param imp_ The underlying ordinary GridFunction.  A reference to
       *             this Object is stored, so the object must be valid for as
       *             long as this GridFunctionBaseAdapter is used.
       */
      GridFunctionBaseAdapter(const Imp& imp_)
        : imp(imp_)
      { }

      //! Evaluate the GridFunction at given position
      /**
       * Evaluates components of the grid function at the given position and
       * returns these values in a vector.
       *
       * \param[in]  e The entity to evaluate on
       * \param[in]  x The position in entity-local coordinates
       * \param[out] y The result of the evaluation
       */
      inline void evaluate (const typename Imp::Traits::ElementType& e,
                            const typename Imp::Traits::DomainType& x,
                            typename Imp::Traits::RangeType& y) const
      {
        imp.evaluate(e,x,y);
      }

      //! get a reference to the GridView
      inline const typename Imp::Traits::GridViewType& getGridView () const
      {
        return imp.getGridView();
      }
    };

	//=======================================
	// helper template for analytic functions
	//=======================================

	//! function signature for analytic functions on a grid
	template<typename GV, typename RF, int m>
	struct AnalyticGridFunctionTraits
	  : public GridFunctionTraits<GV, RF, m, Dune::FieldVector<RF,m> >
	{
	};

	/** \brief an analytic grid function
     *
     *  This is a convenience class which eases the creation of analytic
     *  GridFunctions.  Classes derived from it need only implement a method
     *  evaluateGlobal(const DomainType &x_global, RangeType &y) to have a
     *  full-fledged GridFunction.
     *
     *  \tparam T   The Traits class
     *  \tparam Imp Class implementing the function.  Imp must be derived from
     *              AnalyticGridFunctionBase in some way
     *              (Barton-Nackman-Trick).
     */
	template<typename T, typename Imp>
	class AnalyticGridFunctionBase
	  : public GridFunctionBase<T,AnalyticGridFunctionBase<T,Imp> >
	{
	public:
	  typedef T Traits;

      //! Construct an Analytic GridFunctionBase given a GridView g_
	  AnalyticGridFunctionBase (const typename Traits::GridViewType& g_) : g(g_) {}

      //! \copydoc GridFunctionBase::evaluate()
	  inline void evaluate (const typename Traits::ElementType& e,
							const typename Traits::DomainType& x,
							typename Traits::RangeType& y) const
	  {
		asImp().evaluateGlobal(e.geometry().global(x),y);
	  }

	  inline const typename Traits::GridViewType& getGridView () const
	  {
		return g;
	  }

	private:
	  const typename Traits::GridViewType& g;
	  Imp& asImp () {return static_cast<Imp &> (*this);}
	  const Imp& asImp () const {return static_cast<const Imp &>(*this);}
	};


    // Adapter takes a vector-valued grid function and provides evaluation
    // of normal flux on the interior of faces.
    template<typename T>
    class NormalFluxGridFunctionAdapter
      : public Dune::PDELab::GridFunctionInterface<Dune::PDELab::GridFunctionTraits<typename T::Traits::GridViewType,
                                                                                    typename T::Traits::RangeFieldType,
                                                                                    1,
                                                                                    Dune::FieldVector<typename T::Traits::RangeFieldType,1>
                                                                                    >,
                                                   NormalFluxGridFunctionAdapter<T> >
      , public TypeTree::LeafNode
    {
    public:
      typedef Dune::PDELab::GridFunctionTraits<typename T::Traits::GridViewType,typename T::Traits::RangeFieldType,1,Dune::FieldVector<typename T::Traits::RangeFieldType,1> > Traits;
      typedef Dune::PDELab::GridFunctionInterface<Traits,NormalFluxGridFunctionAdapter<T> > BaseT;

      NormalFluxGridFunctionAdapter (const T& t_) : t(stackobject_to_shared_ptr(t_)) {}


      inline void evaluate (const typename Traits::ElementType& e,
                            const typename Traits::DomainType& x,
                            typename Traits::RangeType& y) const
      {
        // ensure correct size
        dune_static_assert((static_cast<int>(T::Traits::GridViewType::dimension)==static_cast<int>(T::Traits::dimRange)),"number of components must equal dimension");

        // evaluate velocity
        typename T::Traits::RangeType v;
        t->evaluate(e,x,v);

        // implementation only handles triangles so far
        if (!e.geometry().type().isTriangle())
          DUNE_THROW(Dune::NotImplemented, "only implemented for triangles");

        // start and end corner in local numbering
        int n0, n1;

        typename Traits::DomainType nu;

        // determine outer unit normal
        if (std::abs(x[0])<1E-10)
          {
            // edge 1
            n0 = 2;
            n1 = 0;

            nu = e.geometry().corner(n1);
            nu -= e.geometry().corner(n0);
            typename Traits::DomainFieldType temp = nu[0];
            nu[0] = nu[1];
            nu[1] = -temp;
            nu /= nu.two_norm();
            y = v[0]*nu[0]+v[1]*nu[1];
            return;
          }

        if (std::abs(x[1])<1E-10)
          {
            // edge 2
            n0 = 0;
            n1 = 1;

            nu = e.geometry().corner(n1);
            nu -= e.geometry().corner(n0);
            typename Traits::DomainFieldType temp = nu[0];
            nu[0] = nu[1];
            nu[1] = -temp;
            nu /= nu.two_norm();
            y = v[0]*nu[0]+v[1]*nu[1];
            return;
          }

        if (std::abs(x[0]+x[1]-1.0)<1E-10)
          {
            // edge 0
            n0 = 1;
            n1 = 2;

            nu = e.geometry().corner(n1);
            nu -= e.geometry().corner(n0);
            typename Traits::DomainFieldType temp = nu[0];
            nu[0] = nu[1];
            nu[1] = -temp;
            nu /= nu.two_norm();
            y = v[0]*nu[0]+v[1]*nu[1];
            return;
          }

        DUNE_THROW(Dune::Exception, "x needs to be on an edge");
      }

      //! get a reference to the GridView
      inline const typename Traits::GridViewType& getGridView () const
      {
        return t->getGridView();
      }

    private:
      shared_ptr<T const> t;
    };

    // Adapter takes a vector-valued grid function and applies
    // backward Piola transformation on each element
    template<typename T>
    class PiolaBackwardAdapter
      : public Dune::PDELab::GridFunctionInterface<typename T::Traits,PiolaBackwardAdapter<T> >
      , public TypeTree::LeafNode
    {
    public:
      typedef typename T::Traits::GridViewType GridViewType;
      typedef typename T::Traits Traits;
      typedef Dune::PDELab::GridFunctionInterface<Traits,PiolaBackwardAdapter<T> > BaseT;

      PiolaBackwardAdapter (const T& t_) : t(stackobject_to_shared_ptr(t_)) {}


      inline void evaluate (const typename Traits::ElementType& e,
                            const typename Traits::DomainType& x,
                            typename Traits::RangeType& y) const
      {
        // evaluate velocity
        typename T::Traits::RangeType v;
        t->evaluate(e,x,v);

        // apply Piola transformation
        Dune::FieldMatrix<typename Traits::DomainFieldType,Traits::dimRange,Traits::dimRange>
          J = e.geometry().jacobianInverseTransposed(x);
        y = 0;
        J.umtv(v,y);
        y *= e.geometry().integrationElement(x);
      }

      //! get a reference to the GridView
      inline const typename Traits::GridViewType& getGridView () const
      {
        return t->getGridView();
      }

    private:
      shared_ptr<T const> t;
    };


	//==========================
	// template metaprograms
	//==========================

    namespace {

      //! implement VisitingFunctor for vtkwriter_tree_addvertexdata
      template<typename VTKWriter>
      struct AddGridFunctionsToVTKWriter
        : public TypeTree::TreeVisitor
        , public TypeTree::DynamicTraversal
      {

        VTKWriter& w;
        const std::string s;

        AddGridFunctionsToVTKWriter(VTKWriter& w_, const std::string & s_) :
          w(w_), s(s_) {}

        template<typename T, typename TreePath>
        void leaf(const T& t, TreePath treePath) {
          std::stringstream name;
          name << s;
          for (std::size_t i=0; i < treePath.size(); ++i)
            name << "_" << treePath.element(i);
          w.addVertexData(new VTKGridFunctionAdapter<T>(t,name.str()));
        }
      };

    } // anonymous namespace

    /** \brief add vertex data from a \ref GridFunctionTree to a VTKWriter
     *
     *  \tparam GV The GridView for the VTKWriter
     *  \tparam T  The \ref GridFunctionTree
     */
	template<typename GV, typename T>
	void vtkwriter_tree_addvertexdata (Dune::VTKWriter<GV>& w, const T& t, std::string s = "data")
	{
      AddGridFunctionsToVTKWriter<Dune::VTKWriter<GV> > visitor(w,s);
      TypeTree::applyToTree(t,visitor);
	}

    //! \} GridFunctionTree

    //! \} Function

  } // namespace PDELab
} // namespace Dune

#endif
