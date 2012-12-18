// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_ISTL_VECTORITERATOR_HH
#define DUNE_PDELAB_BACKEND_ISTL_VECTORITERATOR_HH

#include <iterator>

#include <dune/common/static_assert.hh>
#include <dune/common/nullptr.hh>
#include <dune/common/tuples.hh>
#include <dune/pdelab/backend/istl/tags.hh>

namespace Dune {

  namespace PDELab {

    namespace istl {

      namespace impl {

#if HAVE_VARIADIC_TEMPLATES

        template<typename T, bool is_const, typename Tag, typename... Iterators>
        struct _extract_iterators;

        template<typename T, typename... Iterators>
        struct _extract_iterators<T,true,tags::block_vector,Iterators...>
          : public _extract_iterators<typename T::block_type,
                                      true,
                                      typename tags::container<typename T::block_type>::type::base_tag,
                                      Iterators..., typename T::const_iterator
                                      >
        {};

        template<typename T, typename... Iterators>
        struct _extract_iterators<T,false,tags::block_vector,Iterators...>
          : public _extract_iterators<typename T::block_type,
                                      false,
                                      typename tags::container<typename T::block_type>::type::base_tag,
                                      Iterators..., typename T::iterator
                                      >
        {};

        template<typename T, typename... Iterators>
        struct _extract_iterators<T,true,tags::field_vector,Iterators...>
        {
          typedef tuple<Iterators...,typename T::const_iterator> type;
        };

        template<typename T, typename... Iterators>
        struct _extract_iterators<T,false,tags::field_vector,Iterators...>
        {
          typedef tuple<Iterators...,typename T::iterator> type;
        };

        template<typename V>
        struct extract_iterators
          : public _extract_iterators<V,false,typename tags::container<V>::type::base_tag>
        {};

        template<typename V>
        struct extract_iterators<const V>
          : public _extract_iterators<V,true,typename tags::container<V>::type::base_tag>
        {};

#endif // HAVE_VARIADIC_TEMPLATES

        template<typename V>
        struct vector_iterator_base
          : public std::iterator<std::forward_iterator_tag,
                                 typename V::field_type,
                                 typename std::ptrdiff_t,
                                 typename V::field_type*,
                                 typename V::field_type&
                                 >
        {
          typedef V vector;
          typedef V& vector_reference;
        };

        template<typename V>
        struct vector_iterator_base<const V>
          : public std::iterator<std::forward_iterator_tag,
                                 typename V::field_type,
                                 typename std::ptrdiff_t,
                                 const typename V::field_type*,
                                 const typename V::field_type&
                                 >
        {
          typedef V vector;
          typedef const V& vector_reference;
        };

      }

      template<typename V>
      struct vector_iterator
        : public impl::vector_iterator_base<V>
      {

        typedef impl::vector_iterator_base<V> BaseT;
        typedef typename BaseT::vector vector;
        typedef typename BaseT::vector_reference vector_reference;
        typedef typename impl::extract_iterators<V>::type Iterators;

        typedef typename tags::container<V>::type::base_tag vector_tag;


      public:

        vector_iterator(vector_reference vector, bool at_end)
          : _vector(vector)
          , _at_end(at_end)
          , _current(nullptr)
        {
          if (!_at_end)
            {
              get<0>(_iterators) = vector.begin();
              get<0>(_end) = vector.end();
              if (!start(vector_tag(),level<0>()))
                _at_end = true;
            }
        }


        typename BaseT::pointer operator->() const
        {
          assert(!_at_end);
          return _current;
        }

        typename BaseT::reference operator*() const
        {
          assert(!_at_end);
          return *_current;
        }

        vector_iterator& operator++()
        {
          increment();
          return *this;
        }

        vector_iterator operator++(int)
        {
          vector_iterator tmp(*this);
          increment();
          return tmp;
        }

        bool operator==(const vector_iterator& r) const
        {
          if (!_at_end)
            {
              if (r._at_end)
                return false;
              return _current == r._current;
            }
          else
            return r._at_end;
        }

        bool operator!=(const vector_iterator& r) const
        {
          return !operator==(r);
        }

      private:

        template<std::size_t l>
        struct level
          : public integral_constant<std::size_t,l>
        {};

        void increment()
        {
          assert(!_at_end);
          if (!advance(vector_tag(),level<0>()))
            _at_end = true;
        }

        template<std::size_t l>
        bool start(tags::field_vector, level<l>)
        {
          typedef typename tuple_element<l,Iterators>::type iterator;
          iterator& it = get<l>(_iterators);
          const iterator& end = get<l>(_end);

          if (it == end)
            return false;

          _current = &(*it);

          return true;
        }

        template<std::size_t l>
        bool start(tags::block_vector, level<l>)
        {
          typedef typename tuple_element<l,Iterators>::type iterator;
          iterator& it = get<l>(_iterators);
          const iterator& end = get<l>(_end);

          while (it != end)
            {
              get<l+1>(_iterators) = it->begin();
              get<l+1>(_end) = it->end();

              if (start(container_tag(*it),level<l+1>()))
                return true;

              ++it;
            }

          return false;
        }


        template<std::size_t l>
        bool advance(tags::field_vector, level<l>)
        {
          typedef typename tuple_element<l,Iterators>::type iterator;
          iterator& it = get<l>(_iterators);
          const iterator& end = get<l>(_end);

          ++it;

          if (it == end)
            return false;

          _current = &(*it);

          return true;
        }


        template<std::size_t l>
        bool advance(tags::block_vector, level<l>)
        {
          typedef typename tuple_element<l,Iterators>::type iterator;
          iterator& it = get<l>(_iterators);

          if (advance(container_tag(*it),level<l+1>()))
            return true;

          ++it;

          return start(tags::block_vector(), level<l>());
        }


        vector_reference _vector;
        bool _at_end;
        typename BaseT::pointer _current;
        Iterators _iterators;
        Iterators _end;

      };

    } // namespace istl
  } // namespace PDELab
} // namespace Dune



#endif // DUNE_PDELAB_BACKEND_ISTL_VECTORITERATOR_HH
