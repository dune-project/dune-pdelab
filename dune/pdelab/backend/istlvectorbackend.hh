#ifndef DUNE_ISTLVECTORBACKEND_HH
#define DUNE_ISTLVECTORBACKEND_HH

#include<dune/istl/bvector.hh>

namespace Dune {
  namespace PDELab {

	//! ISTL backend for FunctionSpace
	template<int BLOCKSIZE=1>
	class ISTLVectorBackend
	{
	public:
	  //! container construction
	  template<typename T, typename E>
	  class VectorContainer : public Dune::BlockVector< Dune::FieldVector<E,BLOCKSIZE> >
	  {
	  public:
		typedef E ElementType;
		typedef Dune::BlockVector< Dune::FieldVector<E,BLOCKSIZE> > BaseT;

		VectorContainer (const T& t) : BaseT(t.globalSize()/BLOCKSIZE) 
		{}
		VectorContainer (const T& t, const E& e) : BaseT(t.globalSize()/BLOCKSIZE) 
		{
		  BaseT::operator=(e);
		}
		VectorContainer& operator= (const E& e)
		{
		  BaseT::operator=(e);
		  return *this;
		}

		// for debugging and AMG access
		BaseT& base ()
		{
		  return *this;
		}

		const BaseT& base () const
		{
		  return *this;
		}
	  };

	  // extract type of container element 
	  template<class C>
	  struct Value
	  {
		typedef typename C::block_type::block_type Type;
	  };

	  //! The size type
	  typedef typename Dune::BlockVector< Dune::FieldVector<float,BLOCKSIZE> >::size_type size_type;

	  // get const_reference to container element
	  // note: this method does not depend on T!
	  template<typename C>
	  static const typename C::field_type& const_access (const C& c, size_type i)
	  {
		return c[i/BLOCKSIZE][i%BLOCKSIZE];
	  }

	  // get non const_reference to container element 
	  // note: this method does not depend on T!
	  template<typename C>
	  static typename C::field_type& access (C& c, size_type i)
	  {
		return c[i/BLOCKSIZE][i%BLOCKSIZE];
	  }
	};

  } // namespace PDELab
} // namespace Dune

#endif
