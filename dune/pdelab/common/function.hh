#ifndef DUNE_PDELAB_FUNCTION_HH
#define DUNE_PDELAB_FUNCTION_HH

#include <iostream>

namespace Dune {
  namespace PDELab {

	// traits class holding function signature, same as in local function
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


	// a Function maps x in DomainType to y in RangeType
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


	// a GridFunction maps x in DomainType to y in RangeType
	template<class T, class Imp>
	class GridFunctionInterface
	{
	public:
	  //! \brief Export type traits
	  typedef T Traits;  

	  /** \brief Evaluate all basis function at given position

		  Evaluates all shape functions at the given position and returns 
		  these values in a vector.
	  */
	  template<typename E>
	  inline void evaluate (const E& e, const typename Traits::DomainType& x,
							typename Traits::RangeType& y) const
	  {  
		asImp().evaluate(e,x,y);
	  }

	private:
	  Imp& asImp () {return static_cast<Imp &> (*this);}
	  const Imp& asImp () const {return static_cast<const Imp &>(*this);}
	};

	// make a GridFunction from a Function
	template<typename T>
	class FunctionToGridFunctionAdapter : 
	  public GridFunctionInterface<typename T::Traits,FunctionToGridFunctionAdapter<T> >
	{
	public:
	  typedef typename T::Traits Traits;

	  FunctionToGridFunctionAdapter (const T& t_) : t(t_) {}

	  template<typename E>
	  inline void evaluate (const E& e, const typename Traits::DomainType& x,
							typename Traits::RangeType& y) const
	  {  
		t.evaluate(e.geometry().global(x),y);
	  }

	private:
	  const T& t; // store reference
	};

	// make a Function in local coordinates from a Function in global coordinates
	template<typename T, typename E>
	class GlobalFunctionToLocalFunctionAdapter : 
	  public FunctionInterface<typename T::Traits,GlobalFunctionToLocalFunctionAdapter<T,E> >
	{
	public:
	  typedef typename T::Traits Traits;

	  GlobalFunctionToLocalFunctionAdapter (const T& t_, const E& e_) : t(t_), e(e_) {}
	  
	  inline void evaluate (const typename Traits::DomainType& x,
							typename Traits::RangeType& y) const
	  {  
		t.evaluate(e.geometry().global(x),y);
	  }

	private:
	  const T& t;
	  const E& e;
	};


	// make a Function from GridFunction using local coordinates
	template<typename T, typename E> // T: GridFunction, E: Entity
	class GridFunctionToLocalFunctionAdapter : 
	  public FunctionInterface<typename T::Traits,GridFunctionToLocalFunctionAdapter<T,E> >
	{
	public:
	  typedef typename T::Traits Traits;

	  GridFunctionToLocalFunctionAdapter (const T& t_, const E& e_) : t(t_), e(e_) {}

	  inline void evaluate (const typename Traits::DomainType& x,
							typename Traits::RangeType& y) const
	  {  
		t.evaluate(e,x,y);
	  }

	private:
	  const T& t;
	  const E& e;
	};

  }
}

#endif
