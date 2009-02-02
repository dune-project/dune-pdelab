#ifndef DUNE_PDELAB_CPSTORAGEPOLICY_HH
#define DUNE_PDELAB_CPSTORAGEPOLICY_HH

#include <iostream>

#include "countingptr.hh"

namespace Dune {
  namespace PDELab {

	class CountingPointerStoragePolicy
	{
	public:
	  template<typename T>
	  struct Storage
	  {
		typedef CP<T> Type;
	  };

	  template<typename T>
	  static CP<T> convert (T& t) // convert to something assignable to storage type
	  {
		return CP<T>(&t);
	  }
	  
	  template<typename T>
	  static void set (CP<T>& s, T& t)
	  {
		s = &t;
	  }
	  
	  template<typename T>
	  static T& get (CP<T>& s)
	  {
		return *s;
	  }

	  template<typename T>
	  static const T& get (const CP<T>& s)
	  {
		return *s;
	  }
	};

  }
}

#endif
