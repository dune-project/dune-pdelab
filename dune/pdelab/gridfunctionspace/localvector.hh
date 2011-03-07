// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_LOCALVECTOR_HH
#define DUNE_PDELAB_LOCALVECTOR_HH

#include <vector>
#include <dune/common/deprecated.hh>

#include "localindex.hh"

/** \file
    \author Christian Engwer
    A local Vector class, which can be tagged by a tag from localfunctionspacetags.hh
 */

namespace Dune {
  namespace PDELab {

#ifndef DOXYGEN

    //! A proxy class to make the new weighted residual and jacobian containers backwards compatible.
    /**
     * This proxy for a single entry of a weighted local vector / local matrix only supports
     * the two operations operator+=() and operator-=(). In both cases, the argument will
     * be multiplied by the weight associated with the corresponding weighted container.
     */
    template<typename T, typename W>
    class WeightedContainerEntryProxy
    {

    public:

      template<typename U>
      void operator+=(const U& r)
      {
        _value += _weight * r;
      }

      template<typename U>
      void operator-=(const U& r)
      {
        _value -= _weight * r;
      }

      WeightedContainerEntryProxy(T& value, W weight)
        : _value(value)
        , _weight(weight)
      {}

    private:

      void operator=(const WeightedContainerEntryProxy&);

      T& _value;
      const W _weight;
    };


#endif // DOXYGEN

    /**
       \addtogroup PDELAB_StrictTrialAndTest Strict Trial and Test space handling
       \ingroup GridFunctionSpace
       \{
    */

    //! An accumulate-only view on a local vector that automatically takes into account an accumulation weight.
    template<typename C>
    class WeightedVectorAccumulationView
    {
    public:

      //! The type of the underlying container.
      typedef C Container;

      //! The value type of the entries.
      typedef typename Container::value_type value_type;

      //! The type of the weight applied when accumulating contributions.
      typedef typename Container::weight_type weight_type;

      //! A special wrapper type to enable backwards compatibility with current containers when directly accessing entries.
      typedef WeightedContainerEntryProxy<value_type,weight_type> reference;

      //! The size_type of the underlying container.
      typedef typename Container::size_type size_type;

      //! Returns a reference proxy to an entry of the underlying container.
      /**
       * \returns A proxy that wraps the entry in the underlying container and
       *          only allows application of operator+=() and operator-=(). In
       *          both cases, the argument will automatically be multiplied by
       *          the weight associated with this view.
       *
       * \deprecated Direct access to the individual container entries is deprecated
       *             and will be removed in a future release. Please use accumulate()
       *             or a combination of weight() and rawAccumulate() to update the
       *             entries of this container.
       */
      reference operator[] (size_type n) DUNE_DEPRECATED
      {
        _modified = true;
        return reference(_container[n],_weight);
      }

      //! Returns the weight associated with this view.
      /**
       * \note This can be used together with rawAccumulate() to avoid applying the weight at
       * each loop iteration.
       */
      weight_type weight()
      {
        _modified = true;
        return _weight;
      }

      //! Resets the weighting coefficient of the view.
      /**
       * \warning Only call this method when you know what you are doing! It is especially not meant
       *          to be called from inside local operators.
       */
      void setWeight(weight_type weight)
      {
        _weight = weight;
      }

      //! Applies the current weight to v and adds the result to the n-th entry of the container.
      void accumulate(size_type n, value_type v)
      {
        _modified = true;
        _container[n] += _weight * v;
      }

      //! Adds v to the n-th entry of the underlying container without applying the current weight.
      /**
       * \warning When using this method, you must take care of applying the weight yourself, otherwise
       *          your program may exhibit strange behavior or even calculate wrong results!
       */
      void rawAccumulate(size_type n, value_type v)
      {
        _container[n] += v;
      }

      //! Constructor
      WeightedVectorAccumulationView(C& container, weight_type weight)
        : _container(container)
        , _weight(weight)
        , _modified(false)
      {}

      //! Returns the size of the underlying container.
      size_type size() const
      {
        return _container.size();
      }

      //! Returns whether this view has been written to.
      bool modified() const
      {
        return _modified;
      }

      //! Resets the modification state of the view to not modified.
      /**
       * \warning Never call this method from within a local operator, or
       *          your local residual / matrix contributions will be lost!
       */
      void resetModified()
      {
        _modified = false;
      }

    private:
      C& _container;
      weight_type _weight;
      bool _modified;
    };

    /**
       \brief a simple container to store local vector entry

       This makes it possible to use strict type checking, even for integers.

       \tparam T   The type of values to store in the vector.
       \tparam TAG Tag type for differentiating between trial and test space vectors.
       \tparam W   The type of weight applied in a WeightedAccumulationView.
     */
    template<typename T, typename TAG, typename W = T>
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
      typedef W weight_type;

      //! An accumulate-only view of the vector that automatically takes into account an accumulation weight.
      typedef WeightedVectorAccumulationView<LocalVector> WeightedAccumulationView;

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

      //! Returns a weighted accumulate-only view of this vector with the given weight.
      WeightedAccumulationView weightedAccumulationView(weight_type weight)
      {
        return WeightedAccumulationView(*this,weight);
      }

    private:
    };

#ifndef DOXYGEN
    // specialization for AnySpaceTag
    template<typename T, typename W>
    class LocalVector<T, AnySpaceTag, W> : public std::vector<T>
    {
    private:
      typedef std::vector<T> BaseT;
    public:
      typedef typename BaseT::value_type  value_type;
      typedef typename BaseT::size_type    size_type;
      typedef typename BaseT::reference    reference;
      typedef typename BaseT::const_reference const_reference;
      typedef W weight_type;
      typedef WeightedVectorAccumulationView<LocalVector> WeightedAccumulationView;

      WeightedAccumulationView weightedAccumulationView(weight_type weight)
      {
        return WeightedAccumulationView(*this,weight);
      }

      /**
	 \{
	 pass contructors to the base class
      */
      LocalVector() {}
      explicit LocalVector(size_type i) : BaseT(i) {}
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
