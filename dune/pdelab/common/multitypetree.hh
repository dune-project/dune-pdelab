// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_MULTITYPETREE_HH
#define DUNE_PDELAB_MULTITYPETREE_HH

#include <iostream>
#include <string>
#include <sstream>

#include <dune/common/static_assert.hh>
#include <dune/common/tuples.hh>
#include <dune/common/typetraits.hh>

#include "multitypetreeutility.hh"

namespace Dune {
  namespace PDELab {

    /** \addtogroup MultiTypeTree MultiTypeTree
     *  \ingroup PDELab
     *  \{
     *
     *  A \ref MultiTypeTree is used to collect different types or multiple
     *  copies of the same type in a tree structure.  This is an invasive
     *  procedure: Types in the tree must be classes or structs derived from
     *  certain base classes, depending on the type of node they represent.
     *
     *  There a three node types:
     *  \li LeafNode is the base class for nodes at the tree's leaf, which do
     *      not have any children.
     *  \li PowerNode is the base class for nodes where all children are of
     *      the same type.
     *  \li CompositeNode is the base class for nodes where the children have
     *      possibly different types.
     *
     *  Because the way the tree is build up from subtrees, child nodes have
     *  no information about their parents.  In fact, as far as the type
     *  system is concerned, a particular node can be a child of multiple
     *  other nodes, or can even be a child of the same node multiple times
     *  (this is the case with the PowerNode).  If needed, this kind of
     *  functionality can be provided on top of the MultiTypeTree, for
     *  instance the GridFunctionSubSpace does something like that.
     *
     *  The MultiTypeTree is more than a simple tree of C++ types -- it can
     *  also store data.  Each non-leaf node can store objects of the types of
     *  its child classes.  How exactly this storing is done is determined by
     *  a StoragePolicy.  For instance, the CopyStoragePolicy stores a copy of
     *  each child object.  The CountingPointerStoragePolicy stores a pointer
     *  to each child object -- multiple child pointers can point to the same
     *  object, saving space.
     *
     *  For storing data the stored child types must be copy-constructible.
     *  The StoragePolicy may add other requirements: CopyStoragePolicy
     *  requires that the children are assignable,
     *  CountingPointerStoragePolicy requires that the children are derived
     *  from Countable.
     */

	//==========================
	// LeafNode
	//==========================

    /** \brief Base class for leaf nodes in the \ref MultiTypeTree
     *
     *  Every leaf type in the \ref MultiTypeTree should be derived from this
     *  class.
     *
     *  A LeafNode models a \ref MultiTypeTree
     */
	class LeafNode
	{
	public:
      //! Mark this class as a leaf in the \ref MultiTypeTree
	  enum { isLeaf = true /**< */ };
      //! Mark this class as a non power in the \ref MultiTypeTree
	  enum { isPower = false /**< */ };
      //! Mark this class as a non composite in the \ref MultiTypeTree
	  enum { isComposite = false /**< */ };
      //! Leafs have no children in the \ref MultiTypeTree
	  enum { CHILDREN = 0 /**< */ };
	};

	//==========================
	// StoragePolicy
	//==========================

    /** \addtogroup StoragePolicy StoragePolicy
     *  \{
     *
     *  In general, there is an object t of type T which should be stored
     *  somehow.  We use a storage object s ("store" for short) of type S to
     *  store t.  In the CopyStoragePolicy T==S, so each store for t will have
     *  its own copy of t.  Other possibilities are S==T& where each store
     *  holds a reference to the original t, or S==T* where each store holds a
     *  pointer to the original t.  It gets interesting when use
     *  S==shared_ptr\<T\> or, in the case of the \ref MultiTypeTree, S=CountingPointer\<T\>
     *  (which is implemented in CountingPointerStoragePolicy).
     */

    /** \brief Default storage policy for the \ref MultiTypeTree
     *
     *  This class determines that elements of the \ref MultiTypeTree are
     *  stored as copies.  For an alternative look at
     *  CountingPointerStoragePolicy, it will all make much more sense then.
     */
	class CopyStoragePolicy
	{
	public:
      /** \brief Determine the storage type S for an object of type T
       *
       * \tparam T The type of the object you want to store
       */
	  template<typename T>
	  struct Storage
	  {
        //! The storage type S for an object of type T is the type T itself
		typedef T Type;
	  };

      //! convert an object of type T to something assignable to its storage type
	  template<typename T>
	  static T& convert (T& t)
	  {
		return t;
	  }
	  
      /** \brief set a store from an object
       *
       *  \param[out] s The store to assign to
       *  \param[in]  t The object to assign
       */
	  template<typename T>
	  static void set (T& s, const T& t)
	  {
		s = t;
	  }
	  
      //! get the object from a store
	  template<typename T>
	  static T& get (T& s)
	  {
		return s;
	  }

      //! get the const object from a const store
	  template<typename T>
	  static const T& get (const T& s)
	  {
		return s;
	  }
	};

    //! \} group StoragePolicy

	//==========================
	// PowerNode
	//==========================

	template<typename T, int k, typename P=CopyStoragePolicy> 
	class PowerNode;

	template<typename T, int k, typename P=CopyStoragePolicy>
	class PowerNodeBase
	{
	public:
	  friend class PowerNode<T,k,P>;

	  enum { isLeaf = false };
	  enum { isPower = true /**< */ };
	  enum { isComposite = false /**< */ };
	  enum { CHILDREN = k };

	  template<int i>
	  struct Child
	  {
		typedef T Type;
	  };

	  template<int i>
	  T& getChild ()
	  {
		return P::get(c[i]);
	  }

	  template<int i>
	  const T& getChild () const
	  {
		return P::get(c[i]);
	  }

	  template<int i>
	  void setChild (T& t)
	  {
		P::set(c[i],t);
	  }

      template<int i>
      typename enable_if<AlwaysTrue<integral_constant<int, i> >::Value &&
                         !IsConst<T>::value>::type
      setChild (const T& t)
	  {
		P::set(c[i],t);
	  }

	  T& getChild (int i)
	  {
		return P::get(c[i]);
	  }

	  const T& getChild (int i) const
	  {
		return P::get(c[i]);
	  }

	  void setChild (int i, T& t)
	  {
		P::set(c[i],t);
	  }

      template<typename U>
      typename enable_if<!IsConst<U>::value>::type
      setChild (int i, const U& t)
	  {
		P::set(c[i],t);
	  }

	private:
	  typename P::template Storage<T>::Type c[k];
	};

    /** \brief Collect k instances of type T within a \ref MultiTypeTree
     *
     *  A PowerNode models a \ref MultiTypeTree
     *
     *  \tparam T The base type
     *  \tparam k The number of instances this node should collect
     *  \tparam P The StoragePolicy to use
     */
	template<typename T, int k, typename P>
	class PowerNode : public PowerNodeBase<T,k,P>
	{
	public:
      //! constructor without arguments
	  PowerNode ()
	  {}

      //! initialize all children with the same object t
	  PowerNode (T& t)
	  {
		for (int i=0; i<k; i++)
		  P::set(this->c[i],t);
	  }

      template<typename U>
      PowerNode (typename enable_if<!IsConst<U>::value, const U>::type& t)
	  {
		for (int i=0; i<k; i++)
		  P::set(this->c[i],t);
	  }

      /** \brief Initialize all children with different objects
       *
       *  This constructor is only available in the non-specialized version
       *
       *  \param t Points to an array of pointers to objects of type T.  The
       *           object pointed to by the first pointer will be used to
       *           initialize the first child, the second pointer for the
       *           second child and so on.
       */
	  PowerNode (T** t)
	  {
		for (int i=0; i<k; i++)
		  P::set(this->c[i],*(t[i]));
	  }

#ifdef DOXYGEN
      /** \brief Initialize all children with different objects
       *
       *  Currently there exist specializations for 2 <= k <= 10.  Each
       *  specialization has a constructor which takes the initializers for
       *  its children as arguments.
       *
       *  @param tn The initializer for the nth child.
       */
      PowerNode (T& t0, T& t1, ...)
      {
      }
#endif // DOXYGEN
	};

	template<typename T, typename P>
	class PowerNode<T,2,P> : public PowerNodeBase<T,2,P>
	{
	public:
	  PowerNode ()
	  {}

	  PowerNode (T& t)
	  {
		for (int i=0; i<2; i++)
		  P::set(this->c[i],t);
	  }

      template<typename U>
      PowerNode (typename enable_if<!IsConst<U>::value, const U>::type& t)
	  {
		for (int i=0; i<2; i++)
		  P::set(this->c[i],t);
	  }

	  PowerNode (T& t0, T& t1)
	  {
		P::set(this->c[0],t0);
		P::set(this->c[1],t1);
	  }

      template<typename U>
      PowerNode (typename enable_if<!IsConst<U>::value, const U>::type& t0,
                 const U& t1)
	  {
		P::set(this->c[0],t0);
		P::set(this->c[1],t1);
	  }
	};

	template<typename T, typename P>
	class PowerNode<T,3,P> : public PowerNodeBase<T,3,P>
	{
	public:
	  PowerNode ()
	  {}

	  PowerNode (T& t)
	  {
		for (int i=0; i<3; i++)
		  P::set(this->c[i],t);
	  }

      template<typename U>
      PowerNode (typename enable_if<!IsConst<U>::value, const U>::type& t)
	  {
		for (int i=0; i<3; i++)
		  P::set(this->c[i],t);
	  }

	  PowerNode (T& t0, T& t1, T& t2)
	  {
		P::set(this->c[0],t0);
		P::set(this->c[1],t1);
		P::set(this->c[2],t2);
	  }

      template<typename U>
      PowerNode (typename enable_if<!IsConst<U>::value, const U>::type& t0,
                 const U& t1, const U& t2)
	  {
		P::set(this->c[0],t0);
		P::set(this->c[1],t1);
		P::set(this->c[2],t2);
	  }
	};

	template<typename T, typename P>
	class PowerNode<T,4,P> : public PowerNodeBase<T,4,P>
	{
	public:
	  PowerNode ()
	  {}

	  PowerNode (T& t)
	  {
		for (int i=0; i<4; i++)
		  P::set(this->c[i],t);
	  }

      template<typename U>
      PowerNode (typename enable_if<!IsConst<U>::value, const U>::type& t)
	  {
		for (int i=0; i<4; i++)
		  P::set(this->c[i],t);
	  }

	  PowerNode (T& t0, T& t1, T& t2, T& t3)
	  {
		P::set(this->c[0],t0);
		P::set(this->c[1],t1);
		P::set(this->c[2],t2);
		P::set(this->c[3],t3);
	  }

      template<typename U>
      PowerNode (typename enable_if<!IsConst<U>::value, const U>::type& t0,
                 const U& t1, const U& t2, const U& t3)
	  {
		P::set(this->c[0],t0);
		P::set(this->c[1],t1);
		P::set(this->c[2],t2);
		P::set(this->c[3],t3);
	  }
	};

	template<typename T, typename P>
	class PowerNode<T,5,P> : public PowerNodeBase<T,5,P>
	{
	public:
	  PowerNode ()
	  {}

	  PowerNode (T& t)
	  {
		for (int i=0; i<5; i++)
		  P::set(this->c[i],t);
	  }

      template<typename U>
      PowerNode (typename enable_if<!IsConst<U>::value, const U>::type& t)
	  {
		for (int i=0; i<5; i++)
		  P::set(this->c[i],t);
	  }

	  PowerNode (T& t0, T& t1, T& t2, T& t3, T& t4)
	  {
		P::set(this->c[0],t0);
		P::set(this->c[1],t1);
		P::set(this->c[2],t2);
		P::set(this->c[3],t3);
		P::set(this->c[4],t4);
	  }

      template<typename U>
      PowerNode (typename enable_if<!IsConst<U>::value, const U>::type& t0,
                 const U& t1, const U& t2, const U& t3, const U& t4)
	  {
		P::set(this->c[0],t0);
		P::set(this->c[1],t1);
		P::set(this->c[2],t2);
		P::set(this->c[3],t3);
		P::set(this->c[4],t4);
	  }
	};

	template<typename T, typename P>
	class PowerNode<T,6,P> : public PowerNodeBase<T,6,P>
	{
	public:
	  PowerNode ()
	  {}

	  PowerNode (T& t)
	  {
		for (int i=0; i<6; i++)
		  P::set(this->c[i],t);
	  }

      template<typename U>
      PowerNode (typename enable_if<!IsConst<U>::value, const U>::type& t)
	  {
		for (int i=0; i<6; i++)
		  P::set(this->c[i],t);
	  }

	  PowerNode (T& t0, T& t1, T& t2, T& t3, T& t4, T& t5)
	  {
		P::set(this->c[0],t0);
		P::set(this->c[1],t1);
		P::set(this->c[2],t2);
		P::set(this->c[3],t3);
		P::set(this->c[4],t4);
		P::set(this->c[5],t5);
	  }

      template<typename U>
      PowerNode (typename enable_if<!IsConst<U>::value, const U>::type& t0,
                 const U& t1, const U& t2, const U& t3, const U& t4,
                 const U& t5)
	  {
		P::set(this->c[0],t0);
		P::set(this->c[1],t1);
		P::set(this->c[2],t2);
		P::set(this->c[3],t3);
		P::set(this->c[4],t4);
		P::set(this->c[5],t5);
	  }
	};

	template<typename T, typename P>
	class PowerNode<T,7,P> : public PowerNodeBase<T,7,P>
	{
	public:
	  PowerNode ()
	  {}

	  PowerNode (T& t)
	  {
		for (int i=0; i<7; i++)
		  P::set(this->c[i],t);
	  }

      template<typename U>
      PowerNode (typename enable_if<!IsConst<U>::value, const U>::type& t)
	  {
		for (int i=0; i<7; i++)
		  P::set(this->c[i],t);
	  }

	  PowerNode (T& t0, T& t1, T& t2, T& t3, T& t4, T& t5, T& t6)
	  {
		P::set(this->c[0],t0);
		P::set(this->c[1],t1);
		P::set(this->c[2],t2);
		P::set(this->c[3],t3);
		P::set(this->c[4],t4);
		P::set(this->c[5],t5);
		P::set(this->c[6],t6);
	  }

      template<typename U>
      PowerNode (typename enable_if<!IsConst<U>::value, const U>::type& t0,
                 const U& t1, const U& t2, const U& t3, const U& t4,
                 const U& t5, const U& t6)
	  {
		P::set(this->c[0],t0);
		P::set(this->c[1],t1);
		P::set(this->c[2],t2);
		P::set(this->c[3],t3);
		P::set(this->c[4],t4);
		P::set(this->c[5],t5);
		P::set(this->c[6],t6);
	  }
	};

	template<typename T, typename P>
	class PowerNode<T,8,P> : public PowerNodeBase<T,8,P>
	{
	public:
	  PowerNode ()
	  {}

	  PowerNode (T& t)
	  {
		for (int i=0; i<8; i++)
		  P::set(this->c[i],t);
	  }

      template<typename U>
      PowerNode (typename enable_if<!IsConst<U>::value, const U>::type& t)
	  {
		for (int i=0; i<8; i++)
		  P::set(this->c[i],t);
	  }

	  PowerNode (T& t0, T& t1, T& t2, T& t3, T& t4, T& t5, T& t6, T& t7)
	  {
		P::set(this->c[0],t0);
		P::set(this->c[1],t1);
		P::set(this->c[2],t2);
		P::set(this->c[3],t3);
		P::set(this->c[4],t4);
		P::set(this->c[5],t5);
		P::set(this->c[6],t6);
		P::set(this->c[7],t7);
	  }

      template<typename U>
      PowerNode (typename enable_if<!IsConst<U>::value, const U>::type& t0,
                 const U& t1, const U& t2, const U& t3, const U& t4,
                 const U& t5, const U& t6, const U& t7)
	  {
		P::set(this->c[0],t0);
		P::set(this->c[1],t1);
		P::set(this->c[2],t2);
		P::set(this->c[3],t3);
		P::set(this->c[4],t4);
		P::set(this->c[5],t5);
		P::set(this->c[6],t6);
		P::set(this->c[7],t7);
	  }
	};

	template<typename T, typename P>
	class PowerNode<T,9,P> : public PowerNodeBase<T,9,P>
	{
	public:
	  PowerNode ()
	  {}

	  PowerNode (T& t)
	  {
		for (int i=0; i<9; i++)
		  P::set(this->c[i],t);
	  }

      template<typename U>
      PowerNode (typename enable_if<!IsConst<U>::value, const U>::type& t)
	  {
		for (int i=0; i<9; i++)
		  P::set(this->c[i],t);
	  }

	  PowerNode (T& t0, T& t1, T& t2, T& t3, T& t4, T& t5, T& t6, T& t7, T& t8)
	  {
		P::set(this->c[0],t0);
		P::set(this->c[1],t1);
		P::set(this->c[2],t2);
		P::set(this->c[3],t3);
		P::set(this->c[4],t4);
		P::set(this->c[5],t5);
		P::set(this->c[6],t6);
		P::set(this->c[7],t7);
		P::set(this->c[8],t8);
	  }

      template<typename U>
      PowerNode (typename enable_if<!IsConst<U>::value, const U>::type& t0,
                 const U& t1, const U& t2, const U& t3, const U& t4,
                 const U& t5, const U& t6, const U& t7, const U& t8)
	  {
		P::set(this->c[0],t0);
		P::set(this->c[1],t1);
		P::set(this->c[2],t2);
		P::set(this->c[3],t3);
		P::set(this->c[4],t4);
		P::set(this->c[5],t5);
		P::set(this->c[6],t6);
		P::set(this->c[7],t7);
		P::set(this->c[8],t8);
	  }
	};

	template<typename T, typename P>
	class PowerNode<T,10,P> : public PowerNodeBase<T,10,P>
	{
	public:
	  PowerNode ()
	  {}

	  PowerNode (T& t)
	  {
		for (int i=0; i<10; i++)
		  P::set(this->c[i],t);
	  }

      template<typename U>
      PowerNode (typename enable_if<!IsConst<U>::value, const U>::type& t)
	  {
		for (int i=0; i<10; i++)
		  P::set(this->c[i],t);
	  }

	  PowerNode (T& t0, T& t1, T& t2, T& t3, T& t4, T& t5, T& t6, T& t7, T& t8, T& t9)
	  {
		P::set(this->c[0],t0);
		P::set(this->c[1],t1);
		P::set(this->c[2],t2);
		P::set(this->c[3],t3);
		P::set(this->c[4],t4);
		P::set(this->c[5],t5);
		P::set(this->c[6],t6);
		P::set(this->c[7],t7);
		P::set(this->c[8],t8);
		P::set(this->c[9],t9);
	  }

      template<typename U>
      PowerNode (typename enable_if<!IsConst<U>::value, const U>::type& t0,
                 const U& t1, const U& t2, const U& t3, const U& t4,
                 const U& t5, const U& t6, const U& t7, const U& t8,
                 const U& t9)
	  {
		P::set(this->c[0],t0);
		P::set(this->c[1],t1);
		P::set(this->c[2],t2);
		P::set(this->c[3],t3);
		P::set(this->c[4],t4);
		P::set(this->c[5],t5);
		P::set(this->c[6],t6);
		P::set(this->c[7],t7);
		P::set(this->c[8],t8);
		P::set(this->c[9],t9);
	  }
	};

	//==========================
	// CompositeNode
	//==========================

	struct EmptyChild {};

	template<typename P, typename OT, typename ST>
	class CompositeNodeBase
	{
	public:

	  enum { isLeaf = false };
	  enum { isPower = false /**< */ };
	  enum { isComposite = true /**< */ };
	  enum { CHILDREN = Dune::tuple_size<OT>::value };

	  CompositeNodeBase (ST& c_) : c(c_) {}

	  template<int i>
	  struct Child
	  {
		typedef typename Dune::tuple_element<i,OT>::type Type;
	  };

	  template<int i>
	  typename Dune::tuple_element<i,OT>::type& getChild ()
	  {
		return P::get(Dune::get<i>(c));
	  }

	  template<int i>
	  const typename Dune::tuple_element<i,OT>::type& getChild () const
	  {
		return P::get(Dune::get<i>(c));
	  }

	  template<int i>
	  void setChild (typename Dune::tuple_element<i,OT>::type& t)
	  {
		P::set(Dune::get<i>(c),t);
	  }	

	private:
	  ST& c;
	};

    /** \brief Collect instances of possibly different types Tn within a \ref
     *         MultiTypeTree
     *
     *  A CompositeNode models a \ref MultiTypeTree
     *
     *  \tparam P  The StoragePolicy to use
     *  \tparam Tn The base types.  Tn==EmptyChild means that slot n is
     *             unused.  Currently, up to 9 slots are supported, making 8
     *             the maximum n.
     */
	template<typename P,
			 typename T0, typename T1, typename T2=EmptyChild, typename T3=EmptyChild,
			 typename T4=EmptyChild, typename T5=EmptyChild, typename T6=EmptyChild,
			 typename T7=EmptyChild, typename T8=EmptyChild>
	class CompositeNode 
#ifndef DOXYGEN
	  : public CompositeNodeBase<P,
								 Dune::tuple<T0,T1,T2,T3,T4,T5,T6,T7,T8>,
								 Dune::tuple<typename P::template Storage<T0>::Type,
											 typename P::template Storage<T1>::Type,
											 typename P::template Storage<T2>::Type,
											 typename P::template Storage<T3>::Type,
											 typename P::template Storage<T4>::Type,
											 typename P::template Storage<T5>::Type,
											 typename P::template Storage<T6>::Type,
											 typename P::template Storage<T7>::Type,
											 typename P::template Storage<T8>::Type> >
#endif //!DOXYGEN
	{
	  typedef Dune::tuple<T0,T1,T2,T3,T4,T5,T6,T7,T8> OT;
#ifndef DOXYGEN
	  typedef Dune::tuple<typename P::template Storage<T0>::Type,
						  typename P::template Storage<T1>::Type,
						  typename P::template Storage<T2>::Type,
						  typename P::template Storage<T3>::Type,
						  typename P::template Storage<T4>::Type,
						  typename P::template Storage<T5>::Type,
						  typename P::template Storage<T6>::Type,
						  typename P::template Storage<T7>::Type,
						  typename P::template Storage<T8>::Type> ST;
#endif //!DOXYGEN
	typedef CompositeNodeBase<P,OT,ST> BaseT;

	ST c;

	public:
	  CompositeNode () : BaseT(c) {}

	  CompositeNode (T0& t0, T1& t1, T2& t2, T3& t3, T4& t4, 
					 T5& t5, T6& t6, T7& t7, T8& t8)
		: BaseT(c), 
		  c(P::convert(t0),P::convert(t1),P::convert(t2),P::convert(t3),P::convert(t4),
			P::convert(t5),P::convert(t6),P::convert(t7),P::convert(t8))
	  {} 

#ifdef DOXYGEN
      /** \brief Initialize all children
       *
       *  @param tn The initializer for the nth child.
       *
       *  The actual number of arguments for this constructor corresponds to
       *  the number of slots used in the template parameter list of the class.
       */
	  CompositeNode (T0& t0, T1& t1, ...) {}
#endif //DOXYGEN
	};

	// 2 children
	template<typename P,
			 typename T0, typename T1>
	class CompositeNode<P,T0,T1,EmptyChild,EmptyChild,EmptyChild,
						EmptyChild,EmptyChild,EmptyChild,EmptyChild>
	  : public CompositeNodeBase<P,
								 Dune::tuple<T0,T1>,
								 Dune::tuple<typename P::template Storage<T0>::Type,
											 typename P::template Storage<T1>::Type> >
	{
	  typedef Dune::tuple<T0,T1> OT;
	  typedef Dune::tuple<typename P::template Storage<T0>::Type,
						  typename P::template Storage<T1>::Type> ST;
      typedef CompositeNodeBase<P,OT,ST> BaseT;
	
	  ST c;

    public:
	  CompositeNode () : BaseT(c) {}

	  CompositeNode (T0& t0, T1& t1)
		: BaseT(c), 
		  c(P::convert(t0),P::convert(t1))
	  {} 
	};

	// 3 children
	template<typename P,
			 typename T0, typename T1, typename T2>
	class CompositeNode<P,T0,T1,T2,EmptyChild,EmptyChild,
						EmptyChild,EmptyChild,EmptyChild,EmptyChild>
	  : public CompositeNodeBase<P,
								 Dune::tuple<T0,T1,T2>,
								 Dune::tuple<typename P::template Storage<T0>::Type,
											 typename P::template Storage<T1>::Type,
											 typename P::template Storage<T2>::Type> >
	{
	  typedef Dune::tuple<T0,T1,T2> OT;
	  typedef Dune::tuple<typename P::template Storage<T0>::Type,
						  typename P::template Storage<T1>::Type,
						  typename P::template Storage<T2>::Type> ST;
      typedef CompositeNodeBase<P,OT,ST> BaseT;
	
	  ST c;

    public:
	  CompositeNode () : BaseT(c) {}
      
      CompositeNode (T0& t0, T1& t1, T2& t2)
		: BaseT(c), 
		  c(P::convert(t0),P::convert(t1),P::convert(t2))
	  {} 
	};

	// 4 children
	template<typename P,
			 typename T0, typename T1, typename T2, typename T3>
	class CompositeNode<P,T0,T1,T2,T3,EmptyChild,
						EmptyChild,EmptyChild,EmptyChild,EmptyChild>
	  : public CompositeNodeBase<P,
								 Dune::tuple<T0,T1,T2,T3>,
								 Dune::tuple<typename P::template Storage<T0>::Type,
											 typename P::template Storage<T1>::Type,
											 typename P::template Storage<T2>::Type,
											 typename P::template Storage<T3>::Type> >
	{
	  typedef Dune::tuple<T0,T1,T2,T3> OT;
	  typedef Dune::tuple<typename P::template Storage<T0>::Type,
						  typename P::template Storage<T1>::Type,
						  typename P::template Storage<T2>::Type,
						  typename P::template Storage<T3>::Type> ST;
      typedef CompositeNodeBase<P,OT,ST> BaseT;
	
	  ST c;

    public:
	  CompositeNode () : BaseT(c) {}
      
      CompositeNode (T0& t0, T1& t1, T2& t2, T3& t3)
		: BaseT(c), 
		  c(P::convert(t0),P::convert(t1),P::convert(t2),P::convert(t3))
	  {} 
	};

	// 5 children
	template<typename P,
			 typename T0, typename T1, typename T2, typename T3, typename T4>
	class CompositeNode<P,T0,T1,T2,T3,T4,
						EmptyChild,EmptyChild,EmptyChild,EmptyChild>
	  : public CompositeNodeBase<P,
								 Dune::tuple<T0,T1,T2,T3,T4>,
								 Dune::tuple<typename P::template Storage<T0>::Type,
											 typename P::template Storage<T1>::Type,
											 typename P::template Storage<T2>::Type,
											 typename P::template Storage<T3>::Type,
											 typename P::template Storage<T4>::Type> >
	{
	  typedef Dune::tuple<T0,T1,T2,T3,T4> OT;
	  typedef Dune::tuple<typename P::template Storage<T0>::Type,
						  typename P::template Storage<T1>::Type,
						  typename P::template Storage<T2>::Type,
						  typename P::template Storage<T3>::Type,
						  typename P::template Storage<T4>::Type> ST;
      typedef CompositeNodeBase<P,OT,ST> BaseT;
	
	  ST c;

    public:
	  CompositeNode () : BaseT(c) {}
      
      CompositeNode (T0& t0, T1& t1, T2& t2, T3& t3, T4& t4)
		: BaseT(c), 
		  c(P::convert(t0),P::convert(t1),P::convert(t2),P::convert(t3),P::convert(t4))
	  {} 
	};

	// 6 children
	template<typename P,
			 typename T0, typename T1, typename T2, typename T3, typename T4,
			 typename T5>
	class CompositeNode<P,T0,T1,T2,T3,T4,
						T5,EmptyChild,EmptyChild,EmptyChild>
	  : public CompositeNodeBase<P,
								 Dune::tuple<T0,T1,T2,T3,T4,T5>,
								 Dune::tuple<typename P::template Storage<T0>::Type,
											 typename P::template Storage<T1>::Type,
											 typename P::template Storage<T2>::Type,
											 typename P::template Storage<T3>::Type,
											 typename P::template Storage<T4>::Type,
											 typename P::template Storage<T5>::Type> >
	{
	  typedef Dune::tuple<T0,T1,T2,T3,T4,T5> OT;
	  typedef Dune::tuple<typename P::template Storage<T0>::Type,
						  typename P::template Storage<T1>::Type,
						  typename P::template Storage<T2>::Type,
						  typename P::template Storage<T3>::Type,
						  typename P::template Storage<T4>::Type,
						  typename P::template Storage<T5>::Type> ST;
      typedef CompositeNodeBase<P,OT,ST> BaseT;
	
	  ST c;

    public:
	  CompositeNode () : BaseT(c) {}
      
      CompositeNode (T0& t0, T1& t1, T2& t2, T3& t3, T4& t4, T5& t5)
		: BaseT(c), 
		  c(P::convert(t0),P::convert(t1),P::convert(t2),P::convert(t3),P::convert(t4),
			P::convert(t5))
	  {} 
	};

	// 7 children
	template<typename P,
			 typename T0, typename T1, typename T2, typename T3, typename T4,
			 typename T5, typename T6>
	class CompositeNode<P,T0,T1,T2,T3,T4,
						T5,T6,EmptyChild,EmptyChild>
	  : public CompositeNodeBase<P,
								 Dune::tuple<T0,T1,T2,T3,T4,T5,T6>,
								 Dune::tuple<typename P::template Storage<T0>::Type,
											 typename P::template Storage<T1>::Type,
											 typename P::template Storage<T2>::Type,
											 typename P::template Storage<T3>::Type,
											 typename P::template Storage<T4>::Type,
											 typename P::template Storage<T5>::Type,
											 typename P::template Storage<T6>::Type> >
	{
	  typedef Dune::tuple<T0,T1,T2,T3,T4,T5,T6> OT;
	  typedef Dune::tuple<typename P::template Storage<T0>::Type,
						  typename P::template Storage<T1>::Type,
						  typename P::template Storage<T2>::Type,
						  typename P::template Storage<T3>::Type,
						  typename P::template Storage<T4>::Type,
						  typename P::template Storage<T5>::Type,
						  typename P::template Storage<T6>::Type> ST;
      typedef CompositeNodeBase<P,OT,ST> BaseT;
	
	  ST c;

    public:
	  CompositeNode () : BaseT(c) {}
      
      CompositeNode (T0& t0, T1& t1, T2& t2, T3& t3, T4& t4, 
				   T5& t5, T6& t6)
		: BaseT(c), 
		  c(P::convert(t0),P::convert(t1),P::convert(t2),P::convert(t3),P::convert(t4),
			P::convert(t5),P::convert(t6))
	  {} 
	};

	// 8 children
	template<typename P,
			 typename T0, typename T1, typename T2, typename T3, typename T4,
			 typename T5, typename T6, typename T7>
	class CompositeNode<P,T0,T1,T2,T3,T4,
						T5,T6,T7,EmptyChild>
	  : public CompositeNodeBase<P,
								 Dune::tuple<T0,T1,T2,T3,T4,T5,T6,T7>,
								 Dune::tuple<typename P::template Storage<T0>::Type,
											 typename P::template Storage<T1>::Type,
											 typename P::template Storage<T2>::Type,
											 typename P::template Storage<T3>::Type,
											 typename P::template Storage<T4>::Type,
											 typename P::template Storage<T5>::Type,
											 typename P::template Storage<T6>::Type,
											 typename P::template Storage<T7>::Type> >
	{
	  typedef Dune::tuple<T0,T1,T2,T3,T4,T5,T6,T7> OT;
	  typedef Dune::tuple<typename P::template Storage<T0>::Type,
						  typename P::template Storage<T1>::Type,
						  typename P::template Storage<T2>::Type,
						  typename P::template Storage<T3>::Type,
						  typename P::template Storage<T4>::Type,
						  typename P::template Storage<T5>::Type,
						  typename P::template Storage<T6>::Type,
						  typename P::template Storage<T7>::Type> ST;
      typedef CompositeNodeBase<P,OT,ST> BaseT;
	
	  ST c;

    public:
	  CompositeNode () : BaseT(c) {}
      
      CompositeNode (T0& t0, T1& t1, T2& t2, T3& t3, T4& t4, 
                     T5& t5, T6& t6, T7& t7)
		: BaseT(c), 
		  c(P::convert(t0),P::convert(t1),P::convert(t2),P::convert(t3),P::convert(t4),
			P::convert(t5),P::convert(t6),P::convert(t7))
	  {} 
	};

    //! \} group MultiTypeTree

  }
}

#endif
