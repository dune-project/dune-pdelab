// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_LOCALVECTOR_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_LOCALVECTOR_HH

#include <vector>
#include <algorithm>
#include <functional>
#include <dune/pdelab/gridfunctionspace/localfunctionspacetags.hh>

/** \file
    \author Christian Engwer
    A local Vector class, which can be tagged by a tag from localfunctionspacetags.hh
 */

namespace Dune {
  namespace PDELab {

    /**
     * \addtogroup PDELab
     * \{
     */

    //! An accumulate-only view on a local vector that automatically takes into account an accumulation weight.
    template<typename C>
    class WeightedVectorAccumulationView
    {
    public:

      //! The type of the underlying LocalVector.
      typedef C Container;

      //! The type of the storage container underlying the LocalVector.
      typedef typename Container::BaseContainer BaseContainer;

      //! The value type of the entries.
      typedef typename Container::value_type value_type;

      //! The type of the weight applied when accumulating contributions.
      typedef typename Container::weight_type weight_type;

      //! \brief Export this type for uniform handling of the containers
      //!        themselves and their views.
      typedef WeightedVectorAccumulationView WeightedAccumulationView;

      //! \brief Returns a WeighedAccumulationView with some weight in
      //!        addition to this view's weight
      WeightedAccumulationView weightedAccumulationView(weight_type weight)
      {
        return WeightedAccumulationView(container(),weight*this->weight());
      }

      //! The size_type of the underlying container.
      typedef typename Container::size_type size_type;

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

      //! Applies the current weight to v and adds the result to the n-th degree of freedom of the lfs.
      template<typename LFS>
      void accumulate(const LFS& lfs, size_type n, value_type v)
      {
        _modified = true;
        _container(lfs,n) += _weight * v;
      }

      //! Adds v to the n-th degree of freedom of the lfs without applying the current weight.
      /**
       * \warning When using this method, you must take care of applying the weight yourself, otherwise
       *          your program may exhibit strange behavior or even calculate wrong results!
       */
      template<typename LFS>
      void rawAccumulate(const LFS& lfs, size_type n, value_type v)
      {
        _modified = true;
        _container(lfs,n) += v;
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

      //! Returns the container (of type LocalVector) that this view is based on.
      Container& container()
      {
        _modified = true;
        return _container;
      }

      //! Returns the container (of type LocalVector) that this view is based on (const version).
      const Container& container() const
      {
        return _container;
      }

      //! Returns the storage container of the underlying LocalVector.
      BaseContainer& base()
      {
        _modified = true;
        return _container.base();
      }

      //! Returns the storage container of the underlying LocalVector (const version).
      const BaseContainer& base() const
      {
        return _container.base();
      }

    private:
      C& _container;
      weight_type _weight;
      bool _modified;
    };


    //! A container for storing data associated with the degrees of freedom of a LocalFunctionSpace.
    /**
     * This container acts as a wrapper around a std::vector-like container and supports accessing
     * its entries indexed by pairs of (LocalFunctionSpace,DOF of LocalFunctionSpace). If requested
     * by specifying a non-default LFSFlavorTag, the container will also assert that a LocalFunctionSpace
     * of the matching kind (trial or test space) is used to access its content.
     *
     * \tparam T            The type of values to store in the vector.
     * \tparam LFSFlavorTag Tag type for differentiating between trial and test space vectors.
     * \tparam W            The type of weight applied in a WeightedAccumulationView.
     */
    template<typename T, typename LFSFlavorTag = AnySpaceTag, typename W = T>
    class LocalVector
    {
    public:

      //! The type of the underlying storage container.
      typedef std::vector<T> BaseContainer;

      //! The value type of this container.
      typedef typename BaseContainer::value_type  value_type;

      //! The size type of this container.
      typedef typename BaseContainer::size_type    size_type;

      //! The reference type of this container.
      typedef typename BaseContainer::reference    reference;

      //! The const reference type of this container.
      typedef typename BaseContainer::const_reference const_reference;

      //! The weight type of this container.
      /**
       * A value of this type will be used to assign a weight to contributions in
       * a WeightedAccumulationView.
       */
      typedef W weight_type;

      //! An accumulate-only view of this container that automatically applies a weight to all contributions.
      typedef WeightedVectorAccumulationView<LocalVector> WeightedAccumulationView;

      //! Returns a WeighedAccumulationView of this container with the given weight.
      WeightedAccumulationView weightedAccumulationView(weight_type weight)
      {
        return WeightedAccumulationView(*this,weight);
      }

      //! Access the value in this container associated with the i-th degree of freedom of the LocalFunctionSpace lfs.
      /**
       * \param lfs The LocalFunctionSpace for which to retrieve a value. This must be the LFS that has been used to
       *            load the values into this vector or one of its children (right now, this is not checked).
       * \param i   The index of the degree of freedom of the LocalFunctionSpace that will be returned.
       */
      template<typename LFS>
      reference operator()(const LFS& lfs, size_type i)
      {
        return _container[lfs.localIndex(i)];
      }

      //! Access the value in this container associated with the i-th degree of freedom of the LocalFunctionSpace lfs (const version).
      /**
       * \param lfs The LocalFunctionSpace for which to retrieve a value. This must be the LFS that has been used to
       *            load the values into this vector or one of its children (right now, this is not checked).
       * \param i   The index of the degree of freedom of the LocalFunctionSpace that will be returned.
       */
      template<typename LFS>
      const_reference operator()(const LFS& lfs, size_type i) const
      {
        return _container[lfs.localIndex(i)];
      }

      //! Assigns v to all entries.
      LocalVector& operator=(const value_type& v)
      {
        std::fill(_container.begin(),_container.end(),v);
        return *this;
      }

      //! Multiplies all entries by v.
      LocalVector& operator*=(const value_type& v)
      {
        using namespace std::placeholders;
        std::transform(
          _container.begin(),
          _container.end(),
          _container.begin(),
          std::bind(std::multiplies<value_type>(),v,_1)
          );
        return *this;
      }

      //! The size of the container.
      size_type size() const
      {
        return _container.size();
      }

      //! Resize the container.
      void resize(size_type size)
      {
        _container.resize(size);
      }

      //! Resize the container to size and assign the passed value to all entries.
      void assign(size_type size, const T& value)
      {
        _container.assign(size,value);
      }

      //! Returns the underlying, std::vector-like storage container.
      BaseContainer& base()
      {
        return _container;
      }

      //! Returns the underlying, std::vector-like storage container (const version).
      const BaseContainer& base() const
      {
        return _container;
      }

      //! Default constructor.
      LocalVector()
      {}

      //! Construct a LocalVector with size n.
      explicit LocalVector(size_type n)
        : _container(n)
      {}

      //! Construct a LocalVector with size n and initialize all entries with v.
      LocalVector(size_type n, const value_type & v)
        : _container(n,v)
      {}

    private:

      BaseContainer _container;

    };


    template<typename C>
    C& accessBaseContainer(C& c)
    {
      return c;
    }

    template<typename T, typename Tag, typename W>
    typename LocalVector<T,Tag,W>::BaseContainer& accessBaseContainer(LocalVector<T,Tag,W>& c)
    {
      return c.base();
    }

    template<typename C>
    typename WeightedVectorAccumulationView<C>::BaseContainer& accessBaseContainer
    (WeightedVectorAccumulationView<C>& c)
    {
      return c.base();
    }

    template<typename C>
    const C& accessBaseContainer(const C& c)
    {
      return c;
    }

    template<typename T, typename Tag, typename W>
    const typename LocalVector<T,Tag,W>::BaseContainer& accessBaseContainer(const LocalVector<T,Tag,W>& c)
    {
      return c.base();
    }

    template<typename C>
    const typename WeightedVectorAccumulationView<C>::BaseContainer& accessBaseContainer
    (const WeightedVectorAccumulationView<C>& c)
    {
      return c.base();
    }

    /**
     * \} group PDELab
     */

  } // end namespace PDELab
} // end namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_LOCALVECTOR_HH
