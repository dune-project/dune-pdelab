#ifndef DUNE_PDELAB_COUNTINGPTR_HH
#define DUNE_PDELAB_COUNTINGPTR_HH

#include <iostream>

namespace Dune {
  namespace PDELab {


	class NondeletingMemoryManagementPolicy
	{
	public:
	  template<typename T>
	  static void delete_action (T* p)
	  {}
	};

	class DeletingMemoryManagementPolicy
	{
	public:
	  template<typename T>
	  static void delete_action (T* p)
	  {
		if (p->reference_counter_zero())
		  delete p;
	  }
	};

	template<typename T, typename P=NondeletingMemoryManagementPolicy>
	class CP
	{
	  T* p;

	public:
	  CP ()
	  {
		p = 0;
	  }

	  explicit CP (T* p_)
	  {
		p = p_;
		if (p!=0)
		  p->reference_counter_increment();
	  }

	  CP (const CP<T>& cp)
	  {
		p = cp.p;
		if (p!=0)
		  p->reference_counter_increment();
	  }

	  ~CP ()
	  {
		if (p!=0)
		  {
			p->reference_counter_decrement();
			P::delete_action(p);
		  }
	  }

	  CP<T>& operator= (T* p_)
	  {
		if (p!=p_)
		  {
			if (p!=0)
			  p->reference_counter_decrement();
			p = p_;
			if (p!=0)
			  p->reference_counter_increment();
		  }
		return *this;
	  }

	  CP<T>& operator= (const CP<T>& cp)
	  {
		if (p!=cp.p)
		  {
			if (p!=0)
			  p->reference_counter_decrement();
			p = cp.p;
			if (p!=0)
			  p->reference_counter_increment();
		  }
		return *this;		
	  }

	  T* operator-> () const
	  {
		return p;
	  }

	  T& operator* () const
	  {
		return *p;
	  }

	  bool operator== (const CP<T>& cp) const
	  {
		return p==cp.p;
	  }

	  bool operator!= (const CP<T>& cp) const
	  {
		return p!=cp.p;
	  }

	};

	class CountableException
	{
	  int counter;
	public:
	  CountableException (int i) : counter(i) {}
	  int get_counter () const
	  {
		return counter;
	  }
	};

	class Countable 
	{
	  mutable int counter;

	public:

	  Countable () : counter(0) 
	  {
		std::cout << "creating object at " << this << std::endl;
	  }

	  // copy constructor: new object, no pointer exists
	  Countable (const Countable& x)
	  {
		std::cout << "creating object at " << this << std::endl;
		counter = 0;
	  } 

	  // number of pointers does not change
	  Countable& operator= (const Countable& x)
	  {
	  } 

	  void reference_counter_increment () const
	  {
		counter++;
	  }

	  void reference_counter_decrement () const
	  {
		counter--;
	  }

	  bool reference_counter_zero () const
	  {
		return counter==0;
	  }

	  int get_reference_counter () const
	  {
		return counter;
	  }

	  ~Countable ()
	  {
		if (counter!=0)
		  {
			std::cout << counter << " counting pointer(s) point to object at "
					  << this << " while it is deleted" << std::endl;
			//			throw CountableException(counter);

		  }
	  }
	};

  }
}

#endif
