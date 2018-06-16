// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_ASSEMBLER_QUADRATURERULE_HH
#define DUNE_PDELAB_ASSEMBLER_QUADRATURERULE_HH

#include <dune/common/iteratorfacades.hh>

#include <dune/geometry/quadraturerules.hh>

namespace Dune {
  namespace PDELab {


    template<typename Rule_>
    class QuadraturePoint
    {

    public:

      using Rule             = Rule_;
      using Native           = typename Rule::Native::value_type;
      using LocalCoordinate  = typename Rule::LocalCoordinate;
      using CellCoordinate   = typename Rule::CellCoordinate;
      using GlobalCoordinate = typename Rule::GlobalCoordinate;
      using Index            = typename Rule::Index;
      using Field            = typename Rule::Field;

      const Native& native() const
      {
        return _qp;
      }

      const LocalCoordinate& local() const
      {
        return _qp.position();
      }

      Index index() const
      {
        return _index;
      }

      Field weight() const
      {
        return _qp.weight() * _rule.integrationElement(*this);
      }

      const CellCoordinate& cell() const
      {
        return inside();
      }

      const CellCoordinate& inside() const
      {
        if (!_inside)
          return _inside.emplace(_rule.inside(local()));
        else
          return *_inside;
      }

      const CellCoordinate& outside() const
      {
        if (!_outside)
          return _outside.emplace(_rule.outside(local()));
        else
          return *_outside;
      }

      const GlobalCoordinate& global() const
      {
        if (!_global)
          return _global.emplace(_rule.global(local()));
        else
          return *_global;
      }

      const Rule& rule() const
      {
        return _rule;
      }

      QuadraturePoint(const Native& qp, Index index, const Rule& rule)
        : _qp(qp)
        , _index(index)
        , _rule(rule)
      {}

    private:

      const Native& _qp;
      Index _index;
      mutable std::optional<CellCoordinate> _inside;
      mutable std::optional<CellCoordinate> _outside;
      mutable std::optional<GlobalCoordinate> _global;
      const Rule& _rule;

    };


    template<typename QR, typename Embedding>
    class QuadratureRule
    {

      friend class Dune::PDELab::QuadraturePoint<QuadratureRule>;

    public:

      using Native           = QR;
      using Field            = typename Native::value_type::Field;
      using LocalCoordinate  = typename Native::value_type::Vector;
      using CellCoordinate   = typename Embedding::Local::GlobalCoordinate;
      using GlobalCoordinate = typename Embedding::Global::GlobalCoordinate;

      //! The size type used by the container.
      using size_type        = typename Native::size_type;
      using Index            = size_type;

      using QuadraturePoint  = Dune::PDELab::QuadraturePoint<QuadratureRule>;
      using value_type       = QuadraturePoint;

      class iterator
        : public RandomAccessIteratorFacade<iterator,QuadraturePoint,QuadraturePoint>
      {

        friend class RandomAccessIteratorFacade<iterator,QuadraturePoint,QuadraturePoint>;

        using NativeIterator = typename Native::const_iterator;

      public:

        // Add support for returning non-references from iterator.
        // We need a little bit of magic to make operator->() work for this iterator
        // because we return a temporary object from dereference(), and the standard
        // implementation of operator->() in the facade tries to take the address of
        // that temporary, which the compiler will vehemently object to... ;-)
        //
        // So I borrowed the following neat little trick from Boost's iterator library:
        // The proxy object stores a copy of the temporary View object, and operator()->
        // returns the proxy object to the caller. As mandated by the standard, the compiler
        // will then attempt to repeat the operator->() on the returned object and get the
        // address of the copy stored in the (temporary) proxy object. That proxy object
        // is guaranteed to live until the next sequence point, and that is precisely as
        // long as we have to guarantee the validity of the pointer to our View object.
        // Problem solved - and another example of how difficult it is to get this low-level
        // stuff implemented on the same level as Boost...
        struct proxy
        {

          explicit proxy(QuadratureRule&& v)
            : _tmp(v)
          {}

          QuadratureRule* operator->()
          {
            return &_tmp;
          }

          QuadratureRule _tmp;
        };

        // The proxy object will stand in as a pointer
        using pointer = proxy;

        iterator()
          : _iterator()
          , _index(0)
          , _rule(nullptr)
        {}

        iterator(NativeIterator it, std::size_t index, const QuadratureRule& rule)
          : _iterator(it)
          , _index(index)
          , _rule(&rule)
        {}

        bool equals(const iterator& other) const
        {
          return _iterator == other._iterator;
        }

        void increment()
        {
          ++_iterator;
          ++_index;
        }

        void decrement()
        {
          --_iterator;
          --_index;
        }

        void advance(int n)
        {
          _iterator += n;
          _index += n;
        }

        std::ptrdiff_t distanceTo(iterator& other) const
        {
          return other._iterator - _iterator;
        }

        QuadraturePoint dereference() const
        {
          return QuadraturePoint(*_iterator,_index,*_rule);
        }

        pointer operator->() const
        {
          return pointer(dereference());
        }

      private:

        NativeIterator _iterator;
        std::size_t _index;
        const QuadratureRule* _rule;

      };

      //! A const iterator over the quadrature points.
      using const_iterator   = iterator;

      //! Returns the maximum polynomial order up to which this rule is exact.
      int order() const
      {
        return _quadrature_rule->order();
      }

      //! Returns the geometry type that this rule is valid for.
      GeometryType type() const
      {
        return _quadrature_rule->type();
      }

      //! Returns the number of quadrature points.
      size_type size() const
      {
        return _quadrature_rule->size();
      }

      //! Returns an iterator pointing to the first quadrature point.
      const_iterator begin() const
      {
        auto it = _quadrature_rule->begin();
        return {it,0,*this};
      }

      //! Returns an iterator pointing after the last quadrature point.
      const_iterator end() const
      {
        return {_quadrature_rule->end(),size(),*this};
      }

#ifndef DOXYGEN

      QuadratureRule(const QR& quadrature_rule, Embedding&& embedding)
        : _quadrature_rule(&quadrature_rule)
        , _embedding(std::move(embedding))
      {}

#endif

    private:

      CellCoordinate inside(const LocalCoordinate& coord) const
      {
        return _embedding.inside().global(coord);
      }

      CellCoordinate outside(const LocalCoordinate& coord) const
      {
        return _embedding.outside().global(coord);
      }

      GlobalCoordinate global(const LocalCoordinate& coord) const
      {
        return _embedding.global().global(coord);
      }

      Field integrationElement(const QuadraturePoint& qp) const
      {
        return _embedding.global().integrationElement(qp.local());
      }

      const QR* _quadrature_rule;
      Embedding _embedding;

    };

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ASSEMBLER_QUADRATURERULE_HH
