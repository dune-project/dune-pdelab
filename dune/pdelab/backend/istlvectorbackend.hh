// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ISTLVECTORBACKEND_HH
#define DUNE_ISTLVECTORBACKEND_HH

#include<vector>

#include<dune/common/fvector.hh>
#include<dune/istl/bvector.hh>

namespace Dune {
  namespace PDELab {

	//! ISTL backend for FunctionSpace
	template<int BLOCKSIZE=1>
	class ISTLVectorBackend
	{
	public:
      enum{
        //! \brief export the block size
        BlockSize = BLOCKSIZE
      };
      
	  //! container construction
	  template<typename T, typename E>
	  class VectorContainer : public Dune::BlockVector< Dune::FieldVector<E,BLOCKSIZE> >
	  {
	  public:
		typedef E ElementType;
		typedef Dune::BlockVector< Dune::FieldVector<E,BLOCKSIZE> > BaseT;
		typedef ISTLVectorBackend<BLOCKSIZE> Backend;
        typedef typename FieldTraits<E>::real_type real_type;
        
		VectorContainer (const T& t_) : BaseT(t_.globalSize()/BLOCKSIZE)
		{}
		VectorContainer (const T& t_, const E& e) : BaseT(t_.globalSize()/BLOCKSIZE)
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

		template<typename X>
		void std_copy_to (std::vector<X>& x) const
		{
		  size_t n = this->size()*BLOCKSIZE;
		  x.resize(n);
		  for (size_t i=0; i<n; i++)
			x[i] = (*this)[i/BLOCKSIZE][i%BLOCKSIZE];
		}

		template<typename X>
		void std_copy_from (const std::vector<X>& x)
		{
		  size_t n = this->size()*BLOCKSIZE;
		  x.resize(n);
		  for (size_t i=0; i<n; i++)
			(*this)[i/BLOCKSIZE][i%BLOCKSIZE] = x[i];
		}
	  };

	  // extract type of container element 
	  template<class C>
	  struct Value
	  {
		typedef typename C::field_type Type;
	  };

	  //! The size type
	  typedef typename Dune::BlockVector< Dune::FieldVector<float,BLOCKSIZE> >::size_type size_type;

	  // get const_reference to container element
	  // note: this method does not depend on T!
	  template<typename C>
	  static const typename C::field_type& access (const C& c, size_type i)
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
