// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_CPSTORAGEPOLICY_HH
#define DUNE_PDELAB_CPSTORAGEPOLICY_HH

#include <iostream>

#include "countingptr.hh"

namespace Dune {
  namespace PDELab {

    //! \addtogroup MultiTypeTree MultiTypeTree
    //! \{
    //!   \addtogroup StoragePolicy
    //!   \{

    /** \brief Storage policy for the \ref MultiTypeTree using CountingPointer objects
     *
     *  This class determines that elements of the \ref MultiTypeTree are
     *  stored as pointers to the original object which are managed be CountingPointer
     *  objects.  See CopyStoragePolicy for a motivation.
     */
	class CountingPointerStoragePolicy
	{
	public:
      /** \brief Determine the storage type S for an object of type T
       *
       * \tparam T The type of the object you want to store
       */
	  template<typename T>
	  struct Storage
	  {
        //! The storage type S for an object of type T is a CountingPointer object
		typedef CountingPointer<T> Type;
	  };

      //! convert an object of type T to something assignable to its storage type
	  template<typename T>
	  static CountingPointer<T> convert (T& t)
	  {
		return CountingPointer<T>(&t);
	  }
	  
      /** \brief set a store from an object
       *
       *  \param[out] s The store to assign to
       *  \param[in]  t The object to assign
       */
	  template<typename T>
	  static void set (CountingPointer<T>& s, T& t)
	  {
		s = &t;
	  }
	  
      //! get the object from a store
	  template<typename T>
	  static T& get (CountingPointer<T>& s)
	  {
		return *s;
	  }

      //! get the const object from a const store
	  template<typename T>
	  static const T& get (const CountingPointer<T>& s)
	  {
		return *s;
	  }
	};

    //!   \} group StoragePolicy
    //! \} group MultiTypeTree

  }
}

#endif
