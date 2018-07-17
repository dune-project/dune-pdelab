// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_ASSEMBLER_LOCALVIEWPROXY_HH
#define DUNE_PDELAB_ASSEMBLER_LOCALVIEWPROXY_HH

#include <tuple>
#include <utility>

namespace Dune {
  namespace PDELab {

    template<typename Field_, typename Weight_>
    class DOFAccumulationProxy
    {

    public:

      using Field  = Field_;
      using Weight = Weight_;

      Weight weight() const
      {
        return _weight;
      }

      template<typename V>
      void operator+=(const V& value)
      {
        *_data += _weight * value;
      }

      template<typename V>
      void operator-=(const V& value)
      {
        *_data -= _weight * value;
      }

      DOFAccumulationProxy(Field& field, Weight weight) noexcept
        : _data(&field)
        , _weight(weight)
      {}

    private:

      Field_* _data;
      Weight _weight;

    };


    template<typename LV, typename LFS>
    class LocalVectorProxy
    {

    public:

      using size_type = typename LV::size_type;
      using value_type = typename LV::value_type;
      using weight_type = typename LV::weight_type;

      using Proxy = DOFAccumulationProxy<value_type,weight_type>;

      class iterator
        : public RandomAccessIteratorFacade<iterator,std::pair<Proxy,size_type>,std::pair<Proxy,size_type>>
      {

        friend class RandomAccessIteratorFacade<iterator,std::pair<Proxy,size_type>,std::pair<Proxy,size_type>>;

      public:

        iterator()
          : _view(nullptr)
          , _weight(0)
          , _index(0)
        {}

        iterator(LocalVectorProxy& view, weight_type weight, std::size_t index)
          : _view(&view)
          , _weight(weight)
          , _index(index)
        {}

        bool equals(const iterator& other) const
        {
          assert(_view == other._view and _weight == other._weight);
          return _index == other._index;
        }

        void increment()
        {
          ++_index;
        }

        void decrement()
        {
          --_index;
        }

        void advance(int n)
        {
          _index += n;
        }

        std::ptrdiff_t distanceTo(iterator& other) const
        {
          return other._index - _index;
        }

        std::pair<Proxy,size_type> dereference() const
        {
          assert(_view);
          return {{(*_view)[_index],_weight},_index};
        }

      private:

        LocalVectorProxy* _view;
        weight_type _weight;
        std::size_t _index;

      };

      iterator begin()
      {
        return {*this,weight(),0};
      }

      iterator end()
      {
        return {*this,weight(),size()};
      }

      weight_type weight() const
      {
        return _weight;
      }

      size_type size() const
      {
        return space().size();
      }

      value_type& operator[](size_type i)
      {
        return view()(space(),i);
      }

      value_type& operator()(const LFS& lfs, size_type i)
      {
        assert(&lfs == _lfs);
        return view()(lfs,i);
      }

      void accumulate(size_type i, const value_type& v)
      {
        view()(space(),i) += v * weight();
      }

      void accumulate(const LFS& lfs, size_type i, const value_type& v)
      {
        assert(&lfs == _lfs);
        view()(lfs,i) += v * weight();
      }

      LV& view()
      {
        return *_lv;
      }

      const LV& view() const
      {
        return *_lv;
      }

      const LFS& space() const
      {
        return *_lfs;
      }

      LocalVectorProxy(LV& lv, const LFS& lfs, weight_type weight)
        : _lv(&lv)
        , _lfs(&lfs)
        , _weight(weight)
      {}

    private:

      LV* _lv;
      const LFS* _lfs;
      weight_type _weight;

    };



    template<typename LM, typename TestSpace, typename TrialSpace>
    class LocalMatrixProxy
    {

    public:

      using size_type = typename LM::size_type;
      using value_type = typename LM::value_type;
      using weight_type = typename LM::weight_type;

      using Proxy = DOFAccumulationProxy<value_type,weight_type>;

      class iterator
        : public RandomAccessIteratorFacade<iterator,std::tuple<Proxy,size_type,size_type>,std::tuple<Proxy,size_type,size_type>>
      {

        friend class RandomAccessIteratorFacade<iterator,std::tuple<Proxy,size_type,size_type>,std::tuple<Proxy,size_type,size_type>>;

      public:

        iterator()
          : _view(nullptr)
          , _weight(0)
          , _columns(0)
          , _index(0)
        {}

        iterator(LocalMatrixProxy& view, weight_type weight, std::size_t columns, std::size_t index)
          : _view(&view)
          , _weight(weight)
          , _columns(columns)
          , _index(index)
        {}

        bool equals(const iterator& other) const
        {
          assert(_view == other._view and _weight == other._weight and _columns == other._columns);
          return _index == other._index;
        }

        void increment()
        {
          ++_index;
        }

        void decrement()
        {
          --_index;
        }

        void advance(int n)
        {
          _index += n;
        }

        std::ptrdiff_t distanceTo(iterator& other) const
        {
          assert(_view == other._view and _weight == other._weight and _columns == other._columns);
          return other._index - _index;
        }

        std::tuple<Proxy,size_type,size_type> dereference() const
        {
          assert(_view);
          size_type i = _index / _columns;
          size_type j = _index % _columns;
          return {{(*_view)(i,j),_weight},i,j};
        }

      private:

        LocalMatrixProxy* _view;
        weight_type _weight;
        std::size_t _columns;
        std::size_t _index;

      };

      iterator begin()
      {
        return {*this,weight(),trialSpace().size(),0};
      }

      iterator end()
      {
        return {*this,weight(),trialSpace().size(),size()};
      }

      weight_type weight() const
      {
        return _weight;
      }

      size_type size() const
      {
        return trialSpace().size() * testSpace().size();
      }

      size_type trialSize() const
      {
        return trialSpace().size();
      }

      size_type testSize() const
      {
        return testSpace().size();
      }

      value_type& operator()(size_type i, size_type j)
      {
        assert(i < testSpace().size());
        assert(j < trialSpace().size());
        return view()(testSpace(),i,trialSpace(),j);
      }

      value_type& operator()(
        const TestSpace& test_space, size_type i,
        const TrialSpace& trial_space, size_type j
        )
      {
        assert(&test_space == _test_space);
        assert(&trial_space == _trial_space);
        return view()(test_space,i,trial_space,j);
      }

      LM& view()
      {
        return *_lm;
      }

      const LM& view() const
      {
        return *_lm;
      }

      const TestSpace& testSpace() const
      {
        return *_test_space;
      }

      const TrialSpace& trialSpace() const
      {
        return *_trial_space;
      }

      void accumulate(size_type i, size_type j, const value_type& v)
      {
        view()(testSpace(),i,trialSpace(),j) += v * weight();
      }

      void accumulate(const TestSpace& test_space, size_type i, const TrialSpace& trial_space, size_type j, const value_type& v)
      {
        assert(&test_space == _test_space);
        assert(&trial_space == _trial_space);
        view()(test_space,i,trial_space,j) += v * weight();
      }

      LocalMatrixProxy(LM& lm, const TestSpace& test_space, const TrialSpace& trial_space, weight_type weight)
        : _lm(&lm)
        , _test_space(&test_space)
        , _trial_space(&trial_space)
        , _weight(weight)
      {}

    private:

      LM* _lm;
      const TestSpace* _test_space;
      const TrialSpace* _trial_space;
      weight_type _weight;

    };


  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ASSEMBLER_LOCALVIEWPROXY_HH
