// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_ASSEMBLER_CELLDATA_HH
#define DUNE_PDELAB_ASSEMBLER_CELLDATA_HH

#include <cassert>
#include <type_traits>
#include <optional>

#include <dune/geometry/identitygeometry.hh>

#include <dune/typetree/childextraction.hh>

#include <dune/pdelab/common/intersectiontype.hh>
#include <dune/pdelab/common/quadraturerules.hh>
#include <dune/pdelab/assembler/utility.hh>
#include <dune/pdelab/assembler/finiteelementwrapper.hh>
#include <dune/pdelab/assembler/quadraturerule.hh>

namespace Dune {
  namespace PDELab {

    namespace CellFlavor {

      struct Inside
      {
        template<typename QP>
        static constexpr auto& quadratureCoordinate(const QP& qp)
        {
          return qp.inside();
        }
      };

      struct Outside
      {
        template<typename QP>
        static constexpr auto& quadratureCoordinate(const QP& qp)
        {
          return qp.outside();
        }
      };

    }

    template<typename Context>
    class CellGridData
      : public Context
    {

    public:

      using EntitySet = typename Context::EntitySet;
      using Entity    = typename EntitySet::template Codim<0>::Entity;
      using Index     = typename EntitySet::IndexSet::Index;
      using Geometry  = typename Entity::Geometry;
      using Embedding = CellEmbedding<Geometry>;
      static constexpr int dimWorld = Geometry::coorddimension;

      const Entity& entity() const
      {
        assert(bound());
        return *_entity;
      }

      Index entityIndex() const {
        assert(bound());
        return _entity_index;
      }

      Index uniqueIndex() const {
        assert(bound());
        return _unique_index;
      }

      const Geometry& geometry() const
      {
        assert(bound());
        if (not _geometry)
          return _geometry.emplace(entity().geometry());
        else
          return *_geometry;
      }

      typename Embedding::GlobalCoordinate global(const typename Embedding::LocalCoordinate& local) const
      {
        return geometry().global(local);
      }

      template<typename P>
      auto global(const P& p) const
      {
        return p.global();
      }

      bool bound() const
      {
        return _entity;
      }

      Embedding embedding() const
      {
        return {geometry()};
      }

      typename EntitySet::Field volume() const
      {
        return geometry().volume();
      }

      EmbeddedPoint<Embedding> centroid() const
      {
        auto ref_el = referenceElement(geometry());
        return {ref_el.position(0,0),embedding()};
      }

      using Context::bind;
      using Context::unbind;

      Context* bind(const Entity& entity, Index entity_index, Index unique_index)
      {
        _entity = &entity;
        _entity_index = entity_index;
        _unique_index = unique_index;
        return this;
      }

      Context* unbind(const Entity& entity, Index entity_index, Index unique_index)
      {
        _entity = nullptr;
        _entity_index = EntitySet::IndexSet::invalidIndex();
        _unique_index = EntitySet::IndexSet::invalidIndex();
        _geometry.reset();
        return this;
      }

      CellGridData(Context&& ctx)
        : Context(std::move(ctx))
        , _entity(nullptr)
        , _entity_index(EntitySet::IndexSet::invalidIndex())
        , _unique_index(EntitySet::IndexSet::invalidIndex())
      {}

    private:

      const Entity* _entity;
      Index _entity_index;
      Index _unique_index;
      mutable std::optional<Geometry> _geometry;

    };

    template<typename Context>
    auto cellGridData(Context&& ctx)
    {
      return CellGridData<Context>{std::move(ctx)};
    }


    template<typename Context>
    class CellDomainData
      : public Context
    {

    public:

      class CellDomain
      {

        friend class CellDomainData;

      public:

        using EntitySet = typename Context::EntitySet;
        using Entity    = typename EntitySet::template Codim<0>::Entity;
        using Field     = typename EntitySet::ctype;

      private:

        using Index = typename EntitySet::IndexSet::Index;

      public:

        // re-use the cell embedding
        using Embedding = typename Context::Embedding;

        Embedding embedding() const
        {
          return {_data.embedding()};
        }

        auto centroid() const
        {
          return _data.centroid();
        }

        const Entity& entity() const
        {
          return _data.inside().entity();
        }

        Field volume() const
        {
          return _data.inside().volume();
        }

        bool bound() const
        {
          return _data.inside().bound();
        }

        auto quadratureRule(std::size_t order, QuadratureType::Enum quadrature_type = QuadratureType::GaussLegendre) const
        {
          auto& rule = QuadratureRules<typename Embedding::Global::ctype,Embedding::Global::mydimension>::rule(embedding().global().type(),order,quadrature_type);
          return QuadratureRule(rule,embedding());
        }

        CellDomain(const CellDomainData& data)
          : _data(data)
        {}

      private:

        const CellDomainData& _data;

      };

      CellDomain cellDomain() const
      {
        return {*this};
      }

      CellDomainData(Context&& ctx)
        : Context(std::move(ctx))
      {}

    };

    template<typename Context>
    auto cellDomainData(Context&& ctx)
    {
      return CellDomainData<Context>{std::move(ctx)};
    }


    template<typename Context>
    class IntersectionDomainData
      : public Context
    {

    public:

      class IntersectionDomain
      {

        friend class IntersectionDomainData;

      public:

        using EntitySet    = typename Context::EntitySet;
        using Entity       = typename EntitySet::template Codim<0>::Entity;
        using Intersection = typename EntitySet::Intersection;
        using Field        = typename EntitySet::Field;

      private:

        using Index = typename EntitySet::IndexSet::Index;
        static constexpr Index invalid_index = EntitySet::IndexSet::invalidIndex();
        static constexpr auto invalid_type   = IntersectionType::invalid;

      public:

        struct Embedding {

          friend class IntersectionDomain;

          using Field   = typename IntersectionDomain::Field;
          using Global  = typename Intersection::Geometry;
          using Cell    = typename Intersection::LocalGeometry;
          using Inside  = typename Intersection::LocalGeometry;
          using Outside = typename Intersection::LocalGeometry;
          using LocalCoordinate = typename Global::LocalCoordinate;
          using GlobalCoordinate = typename Global::GlobalCoordinate;
          using CellCoordinate = typename Inside::GlobalCoordinate;

          const Global& global() const
          {
            return _data.intersectionGeometry();
          }

          const Inside& inside() const
          {
            return _data.intersectionGeometryInInside();
          }

          const Outside& outside() const
          {
            return _data.intersectionGeometryInOutside();
          }

        private:

          Embedding(const IntersectionDomainData& data)
            : _data(data)
          {}

          const IntersectionDomainData& _data;

        };

        Embedding embedding() const
        {
          assert(bound());
          return {_data};
        }

        EmbeddedPoint<Embedding> centroid() const
        {
          auto ref_el = referenceElement(embedding().global());
          return {ref_el.position(0,0),embedding()};
        }

        const Intersection& intersection() const
        {
          return _data.intersection();
        }

        typename Context::Inside& inside()
        {
          return _data.inside();
        }

        typename Context::Outside& outside()
        {
          return _data.outside();
        }

        const Intersection& entity() const
        {
          return intersection();
        }

        Index index() const
        {
          return _data.intersectionIndex();
        }

        IntersectionType type() const
        {
          return _data.intersectionType();
        }

        Field volume() const
        {
          return embedding().global().volume();
        }

        typename Embedding::GlobalCoordinate centerUnitOuterNormal() const
        {
          return intersection().centerUnitOuterNormal();
        }

        template<typename P>
        typename Embedding::GlobalCoordinate unitOuterNormal(const P& p) const
        {
          return intersection().unitOuterNormal(p.local());
        }

        bool bound() const
        {
          return _data._intersection;
        }

        auto quadratureRule(std::size_t order, QuadratureType::Enum quadrature_type = QuadratureType::GaussLegendre) const
        {
          auto& rule = QuadratureRules<typename Embedding::Global::ctype,Embedding::Global::mydimension>::rule(embedding().global().type(),order,quadrature_type);
          return QuadratureRule(rule,embedding());
        }

        IntersectionDomain(const IntersectionDomainData& data)
          : _data(data)
        {}

      private:

        const IntersectionDomainData& _data;

      };

      using Domain = IntersectionDomain;

      IntersectionDomain intersectionDomain() const
      {
        return {*this};
      }

      Domain domain() const
      {
        return {*this};
      }

      const typename IntersectionDomain::Intersection& intersection() const
      {
        assert(_intersection);
        return *_intersection;
      }

      IntersectionType intersectionType() const
      {
        assert(_intersection);
        return _intersection_type;
      }

      typename IntersectionDomain::Index intersectionIndex() const
      {
        assert(_intersection);
        return _intersection_index;
      }

      using Context::bind;
      using Context::unbind;

      template<typename IntersectionType>
      Context* bind(
        IntersectionType type,
        const typename IntersectionDomain::Intersection& intersection,
        typename IntersectionDomain::Index index,
        const typename IntersectionDomain::Entity&,
        typename IntersectionDomain::Index,
        typename IntersectionDomain::Index
        )
      {
        _intersection = &intersection;
        _intersection_index = index;
        _intersection_type = type;
        return this;
      }

      template<typename IntersectionType>
      Context* unbind(
        IntersectionType type,
        const typename IntersectionDomain::Intersection& intersection,
        typename IntersectionDomain::Index index,
        const typename IntersectionDomain::Entity&,
        typename IntersectionDomain::Index,
        typename IntersectionDomain::Index
        )
      {
        _intersection = nullptr;
        _intersection_index = IntersectionDomain::invalid_index;
        _intersection_type = IntersectionDomain::invalid_type;
        _geometry.reset();
        _geometry_in_inside.reset();
        _geometry_in_outside.reset();
        return this;
      }

      IntersectionDomainData(Context&& ctx)
        : Context(std::move(ctx))
        , _intersection(nullptr)
        , _intersection_index(IntersectionDomain::invalid_index)
        , _intersection_type(IntersectionDomain::invalid_type)
      {}

      const typename IntersectionDomain::Embedding::Global& intersectionGeometry() const
      {
        assert(_intersection);
        if (not _geometry)
          return _geometry.emplace(_intersection->geometry());
        else
          return *_geometry;
      }

      const typename IntersectionDomain::Embedding::Inside& intersectionGeometryInInside() const
      {
        assert(_intersection);
        if (not _geometry_in_inside)
          return _geometry_in_inside.emplace(_intersection->geometryInInside());
        else
          return *_geometry_in_inside;
      }

      const typename IntersectionDomain::Embedding::Outside& intersectionGeometryInOutside() const
      {
        assert(_intersection);
        if (not _geometry_in_outside)
          return _geometry_in_outside.emplace(_intersection->geometryInOutside());
        else
          return *_geometry_in_outside;
      }

    private:

      const typename IntersectionDomain::Intersection* _intersection;
      typename IntersectionDomain::Index _intersection_index;
      IntersectionType _intersection_type;
      mutable std::optional<typename IntersectionDomain::Embedding::Global> _geometry;
      mutable std::optional<typename IntersectionDomain::Embedding::Inside> _geometry_in_inside;
      mutable std::optional<typename IntersectionDomain::Embedding::Outside> _geometry_in_outside;

    };

    template<typename Context>
    auto intersectionDomainData(Context&& ctx)
    {
      return IntersectionDomainData<Context>{std::move(ctx)};
    }


    template<typename Context>
    class InsideCell
      : public Context
    {

      struct InsideTraits
      {
        using EntitySet = typename Context::EntitySet;
        using Element   = typename EntitySet::template Codim<0>::Entity;
        using Index     = typename EntitySet::IndexSet::Index;
      };

    public:

      using Inside    = InsideCell;
      using Cell      = InsideCell;

      InsideCell(Context&& ctx)
        : Context(std::move(ctx))
      {}

      Cell& cell()
      {
        return *this;
      }

      Inside& inside()
      {
        return *this;
      }

      const Inside& inside() const
      {
        return *this;
      }

    };

    template<typename Context>
    auto insideCell(Context&& ctx)
    {
      return InsideCell<Context>{std::move(ctx)};
    }


    template<typename Context, typename CellContext>
    class OutsideCell
      : public Context
    {

      struct OutsideTraits
      {
        using EntitySet    = typename CellContext::EntitySet;
        using Entity       = typename EntitySet::template Codim<0>::Entity;
        using Intersection = typename EntitySet::Intersection;
        using Index        = typename EntitySet::IndexSet::Index;
      };

    public:

      using Outside = CellContext;

    private:

      Outside _outside;

    public:

      OutsideCell(Context&& ctx, CellContext&& cell_ctx)
        : Context(std::move(ctx))
        , _outside(std::move(cell_ctx))
      {}

      Outside& outside()
      {
        return _outside;
      }

      const Outside& outside() const
      {
        return _outside;
      }

      Context* setup()
      {
        Dune::PDELab::Context::setup(_outside);
        return this;
      }

      using Context::bind;
      using Context::unbind;

      template<typename IntersectionType>
      Context* bind(
        IntersectionType type,
        const typename OutsideTraits::Intersection&,
        typename OutsideTraits::Index,
        const typename OutsideTraits::Entity& entity,
        typename OutsideTraits::Index entity_index,
        typename OutsideTraits::Index unique_index
        )
      {
        doBind(type,entity,entity_index,unique_index);
        return this;
      }

      template<typename IntersectionType>
      Context* unbind(
        IntersectionType type,
        const typename OutsideTraits::Intersection&,
        typename OutsideTraits::Index,
        const typename OutsideTraits::Entity& entity,
        typename OutsideTraits::Index entity_index,
        typename OutsideTraits::Index unique_index
        )
      {
        doUnbind(type,entity,entity_index,unique_index);
        return this;
      }

    private:

      template<typename IntersectionType>
      void doBind(
        IntersectionType,
        const typename OutsideTraits::Entity& outside_element,
        typename OutsideTraits::Index outside_entity_index,
        typename OutsideTraits::Index outside_unique_index
        )
      {}

      void doBind(
        std::integral_constant<IntersectionType,IntersectionType::skeleton>,
        const typename OutsideTraits::Entity& outside_element,
        typename OutsideTraits::Index outside_entity_index,
        typename OutsideTraits::Index outside_unique_index
        )
      {
        Dune::PDELab::Context::bind(_outside,outside_element,outside_entity_index,outside_unique_index);
      }

      void doBind(
        std::integral_constant<IntersectionType,IntersectionType::periodic>,
        const typename OutsideTraits::Entity& outside_element,
        typename OutsideTraits::Index outside_entity_index,
        typename OutsideTraits::Index outside_unique_index
        )
      {
        Dune::PDELab::Context::bind(_outside,outside_element,outside_entity_index,outside_unique_index);
      }

      template<typename IntersectionType>
      void doUnbind(
        IntersectionType,
        const typename OutsideTraits::Entity& outside_element,
        typename OutsideTraits::Index outside_entity_index,
        typename OutsideTraits::Index outside_unique_index
        )
      {}

      void doUnbind(
        std::integral_constant<IntersectionType,IntersectionType::skeleton>,
        const typename OutsideTraits::Entity& outside_element,
        typename OutsideTraits::Index outside_entity_index,
        typename OutsideTraits::Index outside_unique_index
        )
      {
        Dune::PDELab::Context::unbind(_outside,outside_element,outside_entity_index,outside_unique_index);
      }

      void doUnbind(
        std::integral_constant<IntersectionType,IntersectionType::periodic>,
        const typename OutsideTraits::Entity& outside_element,
        typename OutsideTraits::Index outside_entity_index,
        typename OutsideTraits::Index outside_unique_index
        )
      {
        Dune::PDELab::Context::unbind(_outside,outside_element,outside_entity_index,outside_unique_index);
      }

    };

    template<typename CellContext, typename Context>
    auto outsideCell(CellContext&& cell_ctx, Context&& ctx)
    {
      return OutsideCell<Context,CellContext>{std::move(ctx),std::move(cell_ctx)};
    }

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ASSEMBLER_CELLDATA_HH