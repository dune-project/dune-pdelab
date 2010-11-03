// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_LOCALVECTOR_HH
#define DUNE_PDELAB_LOCALVECTOR_HH

#include <vector>

#include "localindex.hh"

/** \file
    \author Christian Engwer
    A local Vector class, which can be tagged by a tag from localfunctionspacetags.hh
 */

namespace Dune {
  namespace PDELab {
    
    /**
       \addtogroup PDELAB_StrictTrialAndTest Strict Trial and Test space handling
       \ingroup GridFunctionSpace
       \{
    */

    /**
       \brief a simple container to store local vector entry

       This makes it possible to use strict type checking, even for integers.
     */
    template<typename T, typename TAG>
    class LocalVector : public std::vector<T>
    {
    private:
      typedef std::vector<T> BaseT; 
    public:
      typedef typename BaseT::value_type  value_type;
      typedef typename BaseT::size_type  v_size_type;
      typedef typename BaseT::reference    reference;
      typedef typename BaseT::const_reference const_reference;
      typedef LocalIndex<v_size_type, TAG> size_type;

      /**
	 \{
	 pass contructors to the base class
      */
      LocalVector() {}
      LocalVector(v_size_type i) : BaseT(i) {}
      LocalVector(v_size_type i, const value_type & v) : BaseT(i,v) {}
      /** \} */

      /**
	 \{
	 Access Operators operator[](LocalIndex), automatically hides
	 operator[](unsigned int)
       */
      reference operator[] (size_type n) { return BaseT::operator[](n.i); }
      const_reference operator[] (size_type n) const { return BaseT::operator[](n.i); }
      /** \} */
    private:
    };

#ifndef DOXYGEN
    // specialization for AnySpaceTag
    template<typename T>
    class LocalVector<T, AnySpaceTag> : public std::vector<T>
    {
    private:
      typedef std::vector<T> BaseT; 
    public:
      typedef typename BaseT::value_type  value_type;
      typedef typename BaseT::size_type    size_type;
      typedef typename BaseT::reference    reference;
      typedef typename BaseT::const_reference const_reference;

      /**
	 \{
	 pass contructors to the base class
      */
      LocalVector() {}
      LocalVector(size_type i) : BaseT(i) {}
      LocalVector(size_type i, const value_type & v) : BaseT(i,v) {}
      /** \} */
    };
#endif

    /**
       \}
     */

  } // end namespace PDELab
} // end namespace Dune

#endif // DUNE_PDELAB_LOCALVECTOR_HH
