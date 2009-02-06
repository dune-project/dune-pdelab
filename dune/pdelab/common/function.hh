// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_FUNCTION_HH
#define DUNE_PDELAB_FUNCTION_HH

#include <iostream>

#include <dune/common/static_assert.hh>
#include <dune/common/typetraits.hh>

#include"countingptr.hh"
#include"multitypetree.hh"
#include"cpstoragepolicy.hh"
#include"vtkexport.hh"

namespace Dune {
  namespace PDELab {

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


	//! traits class holding function signature, same as in local function
	template<class GV, class RF, int m, class R>
	struct GridFunctionTraits
	  : public FunctionTraits<typename GV::Grid::ctype, GV::dimension,
				    		  Dune::FieldVector<typename GV::Grid::ctype,
                                                GV::dimension>,
							  RF, m, R>
	{
	  //! \brief Export grid view type in addition
	  typedef GV GridViewType;
	  
	  //! \brief codim 0 entity
	  typedef typename GV::Traits::template Codim<0>::Entity ElementType;
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
	  inline const typename Traits::GridViewType& getGridView ()
	  {
		return asImp().getGridView();
	  }

	private:
	  Imp& asImp () {return static_cast<Imp &> (*this);}
	  const Imp& asImp () const {return static_cast<const Imp &>(*this);}
	};


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

	  inline const typename Traits::GridViewType& getGridView ()
	  {
		return g;
	  }

	private:
	  const G& g;
	  const T& t;
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

	//============================
	// Function tree
	//============================

	/** \brief leaf of a function tree
     *
     *  Classes derived from this class implement a \ref GridFunctionTree.
     *
     *  \tparam T   Traits class holding the functions signature
     *  \tparam Imp Class implementing the function.  Imp must be derived from
     *              GridFunctionBase in some way (Barton-Nackman-Trick).
     */
	template<class T, class Imp>
	class GridFunctionBase : public GridFunctionInterface<T,Imp>, 
							 public Countable, 
							 public LeafNode
	{
	public:
      //! Type of the GridView
	  typedef typename T::GridViewType GridViewType;
	};


	/** \brief product of identical functions
     *
     *  This collects k instances of T in a \ref GridFunctionTree.
     *
     *  \tparam T The type of the children of this node in the tree.
     *  \tparam k The number of children this node has.
     */
	template<class T, int k>
	class PowerGridFunction : public PowerNode<T,k,CountingPointerStoragePolicy>,
							  public Countable
	{
	public:
      //! record the GridView
	  typedef typename T::GridViewType GridViewType;

      //! Construct a PowerGridFunction with k clones of the function t
	  PowerGridFunction (T& t) : PowerNode<T,k,CountingPointerStoragePolicy>(t) {}

      /** \brief Initialize all children with different function objects
       *
       *  This constructor is only available in the non-specialized version
       *
       *  \param t Points to an array of pointers to function objects of type
       *           T.  The function pointed to by the first pointer will be
       *           used to initialize the first child, the second pointer for
       *           the second child and so on.
       */
	  PowerGridFunction (T** t) : PowerNode<T,k,CountingPointerStoragePolicy>(t) {}

#ifdef DOXYGEN
      /** \brief Initialize all children with different function objects
       *
       *  Currently there exist specializations for 2 <= k <= 9.  Each
       *  specialization has a constructor which takes the initializers for
       *  its children as arguments.
       *
       *  @param tn The initializer for the nth child.
       */
      PowerGridFunction (T& t0, T& t1, ...)
      {
      }
#endif // DOXYGEN
	};

    // for k=2
	template<class T>
	class PowerGridFunction<T,2> : public PowerNode<T,2,CountingPointerStoragePolicy>,
								   public Countable
	{
	public:
	  typedef typename T::GridViewType GridViewType;

	  PowerGridFunction (T& t) 
		: PowerNode<T,2,CountingPointerStoragePolicy>(t) 
	  {}

	  PowerGridFunction (T& t0, T& t1) 
		: PowerNode<T,2,CountingPointerStoragePolicy>(t0,t1) 
	  {}
	};

    // for k=3
	template<class T>
	class PowerGridFunction<T,3> : public PowerNode<T,3,CountingPointerStoragePolicy>,
								   public Countable
	{
	public:
	  typedef typename T::GridViewType GridViewType;

	  PowerGridFunction (T& t) 
		: PowerNode<T,3,CountingPointerStoragePolicy>(t) 
	  {}

	  PowerGridFunction (T& t0, T& t1, T& t2) 
		: PowerNode<T,3,CountingPointerStoragePolicy>(t0,t1,t2) 
	  {}
	};

    // for k=4
	template<class T>
	class PowerGridFunction<T,4> : public PowerNode<T,4,CountingPointerStoragePolicy>,
								   public Countable
	{
	public:
	  typedef typename T::GridViewType GridViewType;

	  PowerGridFunction (T& t) 
		: PowerNode<T,4,CountingPointerStoragePolicy>(t) 
	  {}

	  PowerGridFunction (T& t0, T& t1, T& t2, T& t3) 
		: PowerNode<T,4,CountingPointerStoragePolicy>(t0,t1,t2,t3) 
	  {}
	};

    // for k=5
	template<class T>
	class PowerGridFunction<T,5> : public PowerNode<T,5,CountingPointerStoragePolicy>,
								   public Countable
	{
	public:
	  typedef typename T::GridViewType GridViewType;

	  PowerGridFunction (T& t) 
		: PowerNode<T,5,CountingPointerStoragePolicy>(t) 
	  {}

	  PowerGridFunction (T& t0, T& t1, T& t2, T& t3, T& t4) 
		: PowerNode<T,5,CountingPointerStoragePolicy>(t0,t1,t2,t3,t4) 
	  {}
	};

    // for k=6
	template<class T>
	class PowerGridFunction<T,6> : public PowerNode<T,6,CountingPointerStoragePolicy>,
								   public Countable
	{
	public:
	  typedef typename T::GridViewType GridViewType;

	  PowerGridFunction (T& t) 
		: PowerNode<T,6,CountingPointerStoragePolicy>(t) 
	  {}

	  PowerGridFunction (T& t0, T& t1, T& t2, T& t3, T& t4, T& t5) 
		: PowerNode<T,6,CountingPointerStoragePolicy>(t0,t1,t2,t3,t4,t5) 
	  {}
	};

    // for k=7
	template<class T>
	class PowerGridFunction<T,7> : public PowerNode<T,7,CountingPointerStoragePolicy>,
								   public Countable
	{
	public:
	  typedef typename T::GridViewType GridViewType;

	  PowerGridFunction (T& t) 
		: PowerNode<T,7,CountingPointerStoragePolicy>(t) 
	  {}

	  PowerGridFunction (T& t0, T& t1, T& t2, T& t3, T& t4, T& t5, T& t6) 
		: PowerNode<T,7,CountingPointerStoragePolicy>(t0,t1,t2,t3,t4,t5,t6) 
	  {}
	};

    // for k=8
	template<class T>
	class PowerGridFunction<T,8> : public PowerNode<T,8,CountingPointerStoragePolicy>,
								   public Countable
	{
	public:
	  typedef typename T::GridViewType GridViewType;

	  PowerGridFunction (T& t) 
		: PowerNode<T,8,CountingPointerStoragePolicy>(t) 
	  {}

	  PowerGridFunction (T& t0, T& t1, T& t2, T& t3, T& t4, T& t5, T& t6, T& t7) 
		: PowerNode<T,8,CountingPointerStoragePolicy>(t0,t1,t2,t3,t4,t5,t6,t7) 
	  {}
	};

    // for k=9
	template<class T>
	class PowerGridFunction<T,9> : public PowerNode<T,9,CountingPointerStoragePolicy>,
								   public Countable
	{
	  typedef typename T::GridViewType GridViewType;

	public:
	  PowerGridFunction (T& t) 
		: PowerNode<T,9,CountingPointerStoragePolicy>(t) 
	  {}

	  PowerGridFunction (T& t0, T& t1, T& t2, T& t3, T& t4, T& t5, T& t6, T& t7, T& t8) 
		: PowerNode<T,9,CountingPointerStoragePolicy>(t0,t1,t2,t3,t4,t5,t6,t7,t8) 
	  {}
	};

    /** \brief composite functions
     *
     *  Collect instances of possibly different function types Tn within a
     *  \ref GridFunctionTree.  This impolements a \ref GridFunctionTree
     *
     *  \tparam Tn The base types.  Tn==EmptyChild means that slot n is
     *             unused.  Currently, up to 9 slots are supported, making 8
     *             the maximum n.
     */
	template<typename T0, typename T1, typename T2=EmptyChild, typename T3=EmptyChild,
			 typename T4=EmptyChild, typename T5=EmptyChild, typename T6=EmptyChild,
			 typename T7=EmptyChild, typename T8=EmptyChild>
	class CompositeGridFunction
	  : public CompositeNode<CountingPointerStoragePolicy,T0,T1,T2,T3,T4,T5,T6,T7,T8>,
		public Countable
	{
	public:
      //! record the GridView
	  typedef typename T0::GridViewType GridViewType;

	  CompositeGridFunction (T0& t0, T1& t1, T2& t2, T3& t3, T4& t4, 
							 T5& t5, T6& t6, T7& t7, T8& t8)
		: CompositeNode<CountingPointerStoragePolicy,
						T0,T1,T2,T3,T4,T5,T6,T7,T8>(t0,t1,t2,t3,t4,t5,t6,t7,t8)
	  {
		dune_static_assert((is_same<typename T0::GridViewType,typename T1::GridViewType>::value),  
						   "GridViewType must be equal in all components of composite grid function");
		dune_static_assert((is_same<typename T0::GridViewType,typename T2::GridViewType>::value),  
						   "GridViewType must be equal in all components of composite grid function");
		dune_static_assert((is_same<typename T0::GridViewType,typename T3::GridViewType>::value),  
						   "GridViewType must be equal in all components of composite grid function");
		dune_static_assert((is_same<typename T0::GridViewType,typename T4::GridViewType>::value),  
						   "GridViewType must be equal in all components of composite grid function");
		dune_static_assert((is_same<typename T0::GridViewType,typename T5::GridViewType>::value),  
						   "GridViewType must be equal in all components of composite grid function");
		dune_static_assert((is_same<typename T0::GridViewType,typename T6::GridViewType>::value),  
						   "GridViewType must be equal in all components of composite grid function");
		dune_static_assert((is_same<typename T0::GridViewType,typename T7::GridViewType>::value),  
						   "GridViewType must be equal in all components of composite grid function");
		dune_static_assert((is_same<typename T0::GridViewType,typename T8::GridViewType>::value),  
						   "GridViewType must be equal in all components of composite grid function");
	  } 

#ifdef DOXYGEN
      /** \brief Initialize all children
       *
       *  @param tn The initializer for the nth child.
       *
       *  The actual number of arguments for this constructor corresponds to
       *  the number of slots used in the template parameter list of the class.
       */
	  CompositeGridFunction (T0& t0, T1& t1, ...) {}
#endif //DOXYGEN
	};


	template<typename T0, typename T1> // 2 children
	class CompositeGridFunction<T0,T1,EmptyChild,EmptyChild,EmptyChild,
								EmptyChild,EmptyChild,EmptyChild,EmptyChild>
	  : public CompositeNode<CountingPointerStoragePolicy,
							 T0,T1,EmptyChild,EmptyChild,EmptyChild,
							 EmptyChild,EmptyChild,EmptyChild,EmptyChild>,
		public Countable

	{
	public:
	  typedef typename T0::GridViewType GridViewType;

	  CompositeGridFunction (T0& t0, T1& t1)
		: CompositeNode<CountingPointerStoragePolicy,
						T0,T1,EmptyChild,EmptyChild,EmptyChild,
						EmptyChild,EmptyChild,EmptyChild,EmptyChild>(t0,t1)
	  {
		dune_static_assert((is_same<typename T0::GridViewType,typename T1::GridViewType>::value),  
						   "GridViewType must be equal in all components of composite grid function");
	  } 
	};
	
	template<typename T0, typename T1, typename T2> // 3 children
	class CompositeGridFunction<T0,T1,T2,EmptyChild,EmptyChild,
								EmptyChild,EmptyChild,EmptyChild,EmptyChild>
	  : public CompositeNode<CountingPointerStoragePolicy,
							 T0,T1,T2,EmptyChild,EmptyChild,
							 EmptyChild,EmptyChild,EmptyChild,EmptyChild>,
		public Countable

	{
	public:
	  typedef typename T0::GridViewType GridViewType;

	  CompositeGridFunction (T0& t0, T1& t1, T2& t2)
		: CompositeNode<CountingPointerStoragePolicy,
						T0,T1,T2,EmptyChild,EmptyChild,
						EmptyChild,EmptyChild,EmptyChild,EmptyChild>(t0,t1,t2)
	  {
		dune_static_assert((is_same<typename T0::GridViewType,typename T1::GridViewType>::value),  
						   "GridViewType must be equal in all components of composite grid function");
		dune_static_assert((is_same<typename T0::GridViewType,typename T2::GridViewType>::value),  
						   "GridViewType must be equal in all components of composite grid function");
	  } 
	};
	
	template<typename T0, typename T1, typename T2, typename T3> // 4 children
	class CompositeGridFunction<T0,T1,T2,T3,EmptyChild,
								EmptyChild,EmptyChild,EmptyChild,EmptyChild>
	  : public CompositeNode<CountingPointerStoragePolicy,
							 T0,T1,T2,T3,EmptyChild,
							 EmptyChild,EmptyChild,EmptyChild,EmptyChild>,
		public Countable

	{
	public:
	  typedef typename T0::GridViewType GridViewType;

	  CompositeGridFunction (T0& t0, T1& t1, T2& t2, T3& t3)
		: CompositeNode<CountingPointerStoragePolicy,
						T0,T1,T2,T3,EmptyChild,
						EmptyChild,EmptyChild,EmptyChild,EmptyChild>(t0,t1,t2,t3)
	  {
		dune_static_assert((is_same<typename T0::GridViewType,typename T1::GridViewType>::value),  
						   "GridViewType must be equal in all components of composite grid function");
		dune_static_assert((is_same<typename T0::GridViewType,typename T2::GridViewType>::value),  
						   "GridViewType must be equal in all components of composite grid function");
		dune_static_assert((is_same<typename T0::GridViewType,typename T3::GridViewType>::value),  
						   "GridViewType must be equal in all components of composite grid function");
	  } 
	};
	
	template<typename T0, typename T1, typename T2, typename T3, typename T4> // 5 children
	class CompositeGridFunction<T0,T1,T2,T3,T4,
								EmptyChild,EmptyChild,EmptyChild,EmptyChild>
	  : public CompositeNode<CountingPointerStoragePolicy,
							 T0,T1,T2,T3,T4,
							 EmptyChild,EmptyChild,EmptyChild,EmptyChild>,
		public Countable

	{
	public:
	  typedef typename T0::GridViewType GridViewType;

	  CompositeGridFunction (T0& t0, T1& t1, T2& t2, T3& t3, T4& t4)
		: CompositeNode<CountingPointerStoragePolicy,
						T0,T1,T2,T3,T4,
						EmptyChild,EmptyChild,EmptyChild,EmptyChild>(t0,t1,t2,t3,t4)
	  {
		dune_static_assert((is_same<typename T0::GridViewType,typename T1::GridViewType>::value),  
						   "GridViewType must be equal in all components of composite grid function");
		dune_static_assert((is_same<typename T0::GridViewType,typename T2::GridViewType>::value),  
						   "GridViewType must be equal in all components of composite grid function");
		dune_static_assert((is_same<typename T0::GridViewType,typename T3::GridViewType>::value),  
						   "GridViewType must be equal in all components of composite grid function");
		dune_static_assert((is_same<typename T0::GridViewType,typename T4::GridViewType>::value),  
						   "GridViewType must be equal in all components of composite grid function");
	  } 
	};
	
	template<typename T0, typename T1, typename T2, typename T3, typename T4,
			 typename T5> // 6 children
	class CompositeGridFunction<T0,T1,T2,T3,T4,
								T5,EmptyChild,EmptyChild,EmptyChild>
	  : public CompositeNode<CountingPointerStoragePolicy,
							 T0,T1,T2,T3,T4,
							 T5,EmptyChild,EmptyChild,EmptyChild>,
		public Countable

	{
	public:
	  typedef typename T0::GridViewType GridViewType;

	  CompositeGridFunction (T0& t0, T1& t1, T2& t2, T3& t3, T4& t4, T5& t5)
		: CompositeNode<CountingPointerStoragePolicy,
						T0,T1,T2,T3,T4,
						T5,EmptyChild,EmptyChild,EmptyChild>(t0,t1,t2,t3,t4,t5)
	  {
		dune_static_assert((is_same<typename T0::GridViewType,typename T1::GridViewType>::value),  
						   "GridViewType must be equal in all components of composite grid function");
		dune_static_assert((is_same<typename T0::GridViewType,typename T2::GridViewType>::value),  
						   "GridViewType must be equal in all components of composite grid function");
		dune_static_assert((is_same<typename T0::GridViewType,typename T3::GridViewType>::value),  
						   "GridViewType must be equal in all components of composite grid function");
		dune_static_assert((is_same<typename T0::GridViewType,typename T4::GridViewType>::value),  
						   "GridViewType must be equal in all components of composite grid function");
		dune_static_assert((is_same<typename T0::GridViewType,typename T5::GridViewType>::value),  
						   "GridViewType must be equal in all components of composite grid function");
	  } 
	};
	
	template<typename T0, typename T1, typename T2, typename T3, typename T4,
			 typename T5, typename T6> // 7 children
	class CompositeGridFunction<T0,T1,T2,T3,T4,
								T5,T6,EmptyChild,EmptyChild>
	  : public CompositeNode<CountingPointerStoragePolicy,
							 T0,T1,T2,T3,T4,
							 T5,T6,EmptyChild,EmptyChild>,
		public Countable

	{
	public:
	  typedef typename T0::GridViewType GridViewType;

	  CompositeGridFunction (T0& t0, T1& t1, T2& t2, T3& t3, T4& t4, T5& t5, T6& t6)
		: CompositeNode<CountingPointerStoragePolicy,
						T0,T1,T2,T3,T4,
						T5,T6,EmptyChild,EmptyChild>(t0,t1,t2,t3,t4,t5,t6)
	  {
		dune_static_assert((is_same<typename T0::GridViewType,typename T1::GridViewType>::value),  
						   "GridViewType must be equal in all components of composite grid function");
		dune_static_assert((is_same<typename T0::GridViewType,typename T2::GridViewType>::value),  
						   "GridViewType must be equal in all components of composite grid function");
		dune_static_assert((is_same<typename T0::GridViewType,typename T3::GridViewType>::value),  
						   "GridViewType must be equal in all components of composite grid function");
		dune_static_assert((is_same<typename T0::GridViewType,typename T4::GridViewType>::value),  
						   "GridViewType must be equal in all components of composite grid function");
		dune_static_assert((is_same<typename T0::GridViewType,typename T5::GridViewType>::value),  
						   "GridViewType must be equal in all components of composite grid function");
		dune_static_assert((is_same<typename T0::GridViewType,typename T6::GridViewType>::value),  
						   "GridViewType must be equal in all components of composite grid function");
	  } 
	};
	
	template<typename T0, typename T1, typename T2, typename T3, typename T4,
			 typename T5, typename T6, typename T7> // 8 children
	class CompositeGridFunction<T0,T1,T2,T3,T4,
								T5,T6,T7,EmptyChild>
	  : public CompositeNode<CountingPointerStoragePolicy,
							 T0,T1,T2,T3,T4,
							 T5,T6,T7,EmptyChild>,
		public Countable

	{
	public:
	  typedef typename T0::GridViewType GridViewType;

	  CompositeGridFunction (T0& t0, T1& t1, T2& t2, T3& t3, T4& t4, T5& t5, T6& t6, T7& t7)
		: CompositeNode<CountingPointerStoragePolicy,
						T0,T1,T2,T3,T4,
						T5,T6,T7,EmptyChild>(t0,t1,t2,t3,t4,t5,t6,t7)
	  {
		dune_static_assert((is_same<typename T0::GridViewType,typename T1::GridViewType>::value),  
						   "GridViewType must be equal in all components of composite grid function");
		dune_static_assert((is_same<typename T0::GridViewType,typename T2::GridViewType>::value),  
						   "GridViewType must be equal in all components of composite grid function");
		dune_static_assert((is_same<typename T0::GridViewType,typename T3::GridViewType>::value),  
						   "GridViewType must be equal in all components of composite grid function");
		dune_static_assert((is_same<typename T0::GridViewType,typename T4::GridViewType>::value),  
						   "GridViewType must be equal in all components of composite grid function");
		dune_static_assert((is_same<typename T0::GridViewType,typename T5::GridViewType>::value),  
						   "GridViewType must be equal in all components of composite grid function");
		dune_static_assert((is_same<typename T0::GridViewType,typename T6::GridViewType>::value),  
						   "GridViewType must be equal in all components of composite grid function");
		dune_static_assert((is_same<typename T0::GridViewType,typename T7::GridViewType>::value),  
						   "GridViewType must be equal in all components of composite grid function");
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

	  inline const typename Traits::GridViewType& getGridView ()
	  {
		return g;
	  }
  
	private:
	  const typename Traits::GridViewType& g;
	  Imp& asImp () {return static_cast<Imp &> (*this);}
	  const Imp& asImp () const {return static_cast<const Imp &>(*this);}
	};


	//==========================
	// template metaprograms
	//==========================

	template<typename T, bool isleaf>
	struct GridFunctionTreeVisitNodeMetaProgram;

	template<typename T, int n, int i>
	struct GridFunctionTreeVisitChildMetaProgram // visit child of inner node
	{
	  template<typename GV> 
	  static void vtkwriter_tree_addvertexdata (Dune::VTKWriter<GV>& w, const T& t, std::string s)
	  {
		typedef typename T::template Child<i>::Type C;
		std::string cs(s);
		cs += "_";
		std::stringstream out;
		out << i;
		cs += out.str();
		GridFunctionTreeVisitNodeMetaProgram<C,C::isLeaf>::vtkwriter_tree_addvertexdata(w,t.template getChild<i>(),cs);
		GridFunctionTreeVisitChildMetaProgram<T,n,i+1>::vtkwriter_tree_addvertexdata(w,t,s);
	  }
	};

	template<typename T, int n>
	struct GridFunctionTreeVisitChildMetaProgram<T,n,n> // end of child recursion
	{
	  template<typename GV> 
	  static void vtkwriter_tree_addvertexdata (Dune::VTKWriter<GV>& w, const T& t, std::string s)
	  {
		return;
	  }
	};

	template<typename T, bool isleaf> 
	struct GridFunctionTreeVisitNodeMetaProgram // visit inner node
	{
	  template<typename GV> 
	  static void vtkwriter_tree_addvertexdata (Dune::VTKWriter<GV>& w, const T& t, std::string s)
	  {
		GridFunctionTreeVisitChildMetaProgram<T,T::CHILDREN,0>::vtkwriter_tree_addvertexdata(w,t,s);
	  }
	};

	template<typename T> 
	struct GridFunctionTreeVisitNodeMetaProgram<T,true> // visit leaf node 
	{
	  template<typename GV> 
	  static void vtkwriter_tree_addvertexdata (Dune::VTKWriter<GV>& w, const T& t, std::string s)
	  {
		w.addVertexData(new VTKGridFunctionAdapter<T>(t,s));
	  }
	};

	template<typename GV, typename T> 
	void vtkwriter_tree_addvertexdata (Dune::VTKWriter<GV>& w, const T& t)
	{
	  std::string s="data";
	  GridFunctionTreeVisitNodeMetaProgram<T,T::isLeaf>::vtkwriter_tree_addvertexdata(w,t,s);
	}
  }
}

#endif
