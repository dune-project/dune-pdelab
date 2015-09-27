// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_COMMON_MULTIINDEX_HH
#define DUNE_PDELAB_COMMON_MULTIINDEX_HH

#include <dune/common/reservedvector.hh>
#include <dune/geometry/typeindex.hh>

#include <dune/common/hash.hh>

#include <algorithm>
#include <iomanip>

namespace Dune {

  namespace PDELab {


    //! A class for representing multi-indices.
    /**
     * A MultiIndex represents an ordered tuple of indices.
     *
     * \tparam T  the type of the index entries.
     * \tparam n  the maximum number of indices in the MultiIndex.
     */
    template<typename T, std::size_t n>
    class MultiIndex
      : public ReservedVector<T,n>
    {

      typedef ReservedVector<T,n> base_type;

    public:

      //! The maximum possible depth of the MultiIndex.
      static const std::size_t max_depth = n;

      class View
      {

        friend class MultiIndex;

      public:

        //! The maximum possible depth of the MultiIndex.
        static const std::size_t max_depth = n;

        typedef typename base_type::value_type value_type;
        typedef typename base_type::pointer pointer;
        typedef typename base_type::const_reference reference;
        typedef typename base_type::const_reference const_reference;
        typedef typename base_type::size_type size_type;
        typedef typename base_type::difference_type difference_type;
        typedef typename base_type::const_iterator iterator;
        typedef typename base_type::const_iterator const_iterator;

      private:

        View(const MultiIndex& mi, size_type size)
          : _mi(mi)
          , _size(size)
        {}

      public:

        void clear()
        {
          _size = 0;
        }

        reference front()
        {
          return _mi.front();
        }

        const_reference front() const
        {
          return _mi.front();
        }

        reference back()
        {
          return _mi[_size-1];
        }

        const_reference back() const
        {
          return _mi[_size-1];
        }

        reference operator[](size_type i)
        {
          assert(i < _size);
          return _mi[i];
        }

        const_reference operator[](size_type i) const
        {
          assert(i < _size);
          return _mi[i];
        }

        void resize(size_type s)
        {
          assert(s <= _mi.size());
          _size = s;
        }

        View back_popped() const
        {
          assert(_size > 0);
          return View(_mi,_size-1);
        }

        size_type size() const
        {
          return _size;
        }

        bool empty() const
        {
          return _size == 0;
        }

        friend std::ostream& operator<< (std::ostream& s, const View& mi)
        {
          s << "(";
          // fill up to maximum depth for consistent formatting
          for (std::size_t i = mi.size(); i < max_depth; ++i)
            s << "  -";
          for (typename ReservedVector<T,n>::const_iterator it = mi._mi.begin(); it != mi._mi.begin() + mi.size(); ++it)
            s << std::setw(3) << *it;
          s << ")";
          return s;
        }

      private:
        const MultiIndex& _mi;
        size_type _size;

      };

      MultiIndex()
      {}

      MultiIndex(const View& view)
        : base_type(static_cast<const base_type&>(view._mi))
      {
        this->resize(view.size());
      }

      void set(typename ReservedVector<T,n>::value_type index)
      {
        this->clear();
        this->push_back(index);
      }

      //! Writes a pretty representation of the MultiIndex to the given std::ostream.
      friend std::ostream& operator<< (std::ostream& s, const MultiIndex& mi)
      {
        s << "(";
        // fill up to maximum depth for consistent formatting
        for (std::size_t i = mi.size(); i < max_depth; ++i)
          s << "  -";
        for (typename ReservedVector<T,n>::const_iterator it = mi.begin(); it != mi.end(); ++it)
          s << std::setw(3) << *it;
        s << ")";
        return s;
      }

      View view() const
      {
        return View(*this,this->size());
      }

      View view(std::size_t size) const
      {
        return View(*this,size);
      }

      //! Tests whether two MultiIndices are equal.
      /**
       * \note Only MultiIndices of identical max_depth are comparable.
       */
      bool operator== (const MultiIndex& r) const
      {
        return
          this->size() == r.size() &&
          std::equal(this->begin(),this->end(),r.begin());
      }

      //! Tests whether two MultiIndices are not equal.
      bool operator!= (const MultiIndex& r) const
      {
        return !(*this == r);
      }

#if 0
      bool operator< (const MultiIndex& r) const
      {
        // FIXME: think about natural ordering
        return _c.size() < _r.size();
        return std::lexicographical_compare(_c.begin(),_c.end(),r._c.begin(),r._c.end());
      }
#endif

    };


    template<typename T, std::size_t n>
    inline std::size_t hash_value(const MultiIndex<T,n>& mi)
    {
      return hash_range(mi.begin(),mi.end());
    }


  } // namespace PDELab
} // namespace Dune

DUNE_DEFINE_HASH(DUNE_HASH_TEMPLATE_ARGS(typename T, std::size_t n),DUNE_HASH_TYPE(Dune::PDELab::MultiIndex<T,n>))

#endif // DUNE_PDELAB_COMMON_MULTIINDEX_HH
