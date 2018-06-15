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

    namespace CellType {

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

    template<typename ES>
    class CellGridData
    {

    public:

      using EntitySet = ES;
      using Entity    = typename ES::template Codim<0>::Entity;
      using Index     = typename ES::IndexSet::Index;
      using Geometry  = typename Entity::Geometry;

      const Entity& entity() const
      {
        assert(_entity);
        return *_entity;
      }

      Index entityIndex() const {
        return _entity_index;
      }

      Index uniqueIndex() const {
        return _unique_index;
      }

      bool bound() const
      {
        return _entity;
      }

      template<typename Context, typename CellContext>
      void bind(Context& ctx, CellContext& cell_ctx, const Entity& entity, Index entity_index, Index unique_index)
      {
        _entity = &entity;
        _entity_index = entity_index;
        _unique_index = unique_index;
      }

      template<typename Context, typename CellContext>
      void unbind(Context& ctx, CellContext& cell_ctx, const Entity& entity, Index entity_index, Index unique_index)
      {
        _entity = nullptr;
        _entity_index = ES::IndexSet::invalidIndex();
        _unique_index = ES::IndexSet::invalidIndex();
      }

      CellGridData()
        : _entity(nullptr)
        , _entity_index(ES::IndexSet::invalidIndex())
        , _unique_index(ES::IndexSet::invalidIndex())
      {}

    private:

      const Entity* _entity;
      Index _entity_index;
      Index _unique_index;

    };

    template<typename ES>
    auto cellGridData(const ES&)
    {
      return CellGridData<ES>{};
    }


    template<typename ES>
    class CellDomainData
    {

    public:

      class CellDomain
      {

        friend class CellDomainData;

      public:

        using EntitySet = ES;
        using Entity    = typename ES::template Codim<0>::Entity;

      private:

        using Index = typename EntitySet::IndexSet::Index;

      public:

        struct Embedding {

          using Global  = typename Entity::Geometry;
          using Local   = IdentityGeometry<typename ES::ctype,ES::dimension>;

          Global global() const
          {
            return _entity->geometry();
          }

          Local local() const
          {
            return Local(_entity->type());
          }

          Local inside() const
          {
            return local();
          }

        private:

          friend class CellDomain;

          Embedding(const Entity& e)
            : _entity(&e)
          {}

          const Entity* _entity;

        };

        Embedding embedding() const
        {
          assert(_entity);
          return {*_entity};
        }

        const Entity& entity() const
        {
          assert(_entity);
          return *_entity;
        }

        bool bound() const
        {
          return _entity;
        }

        auto quadratureRule(std::size_t order, QuadratureType::Enum quadrature_type = QuadratureType::GaussLegendre) const
        {
          auto& rule = QuadratureRules<typename Embedding::Global::ctype,Embedding::Global::mydimension>::rule(embedding().global().type(),order,quadrature_type);
          return QuadratureRule(rule,embedding());
        }

        CellDomain()
          : _entity(nullptr)
        {}

      private:

        const Entity* _entity;

      };

      CellDomain& cellDomain()
      {
        return _cell_domain;
      }

      const CellDomain& cellDomain() const
      {
        return _cell_domain;
      }

      CellDomain& domain()
      {
        return _cell_domain;
      }

      const CellDomain& domain() const
      {
        return _cell_domain;
      }

      template<typename Context>
      void bind(
        Context& ctx,
        const typename CellDomain::Entity& entity,
        typename CellDomain::Index entity_index,
        typename CellDomain::Index unique_index
        )
      {
        _cell_domain._entity = &entity;
      }

      template<typename Context>
      void unbind(
        Context& ctx,
        const typename CellDomain::Entity& entity,
        typename CellDomain::Index entity_index,
        typename CellDomain::Index unique_index
        )
      {
        _cell_domain._entity = nullptr;
      }

    private:

      CellDomain _cell_domain;

    };

    template<typename ES>
    auto cellDomainData(const ES&)
    {
      return CellDomainData<ES>{};
    }


    template<typename ES>
    class IntersectionDomainData
    {

    public:

      class IntersectionDomain
      {

        friend class IntersectionDomainData;

      public:

        using EntitySet    = ES;
        using Entity       = typename ES::template Codim<0>::Entity;
        using Intersection = typename ES::Intersection;

      private:

        using Cell  = typename ES::template Codim<0>::Entity;
        using Index = typename EntitySet::IndexSet::Index;
        static constexpr Index invalid_index = EntitySet::IndexSet::invalidIndex();
        static constexpr auto invalid_type   = IntersectionType::invalid;

      public:

        struct Embedding {

          using Global  = typename Intersection::Geometry;
          using Inside  = typename Intersection::LocalGeometry;
          using Outside = typename Intersection::LocalGeometry;

          Global global() const
          {
            return _intersection->geometry();
          }

          Inside inside() const
          {
            return _intersection->geometryInInside();
          }

          Outside outside() const
          {
            return _intersection->geometryInOutside();
          }

        private:

          friend class IntersectionDomain;

          Embedding(const Intersection& intersection)
            : _intersection(&intersection)
          {}

          const Intersection* _intersection;

        };

        Embedding embedding() const
        {
          assert(bound());
          return {*_intersection};
        }

        const Intersection& intersection() const
        {
          assert(bound());
          return *_intersection;
        }

        const Intersection& entity() const
        {
          return intersection();
        }

        Index index() const
        {
          assert(bound());
          return _index;
        }

        bool bound() const
        {
          return _intersection;
        }

        auto quadratureRule(std::size_t order, QuadratureType::Enum quadrature_type = QuadratureType::GaussLegendre) const
        {
          assert(bound());
          auto& rule = QuadratureRules<typename Embedding::Geometry::ctype,Embedding::Geometry::mydimension>::rule(embedding().global().type(),order,quadrature_type);
          return QuadratureRule(rule,embedding());
        }

        IntersectionDomain()
          : _intersection(nullptr)
          , _index(invalid_index)
          , _type(invalid_type)
        {}

      private:

        const Intersection* _intersection;
        Index _index;
        IntersectionType _type;

      };

      IntersectionDomain& intersectionDomain()
      {
        return _intersection_domain;
      }

      const IntersectionDomain& intersectionDomain() const
      {
        return _intersection_domain;
      }

      template<typename Context, typename IntersectionType>
      void bind(
        Context& ctx,
        IntersectionType type,
        const typename IntersectionDomain::Intersection& intersection,
        typename IntersectionDomain::Index index,
        const typename IntersectionDomain::Cell&,
        typename IntersectionDomain::Index,
        typename IntersectionDomain::Index
        )
      {
        _intersection_domain._intersection = &intersection;
        _intersection_domain._index = index;
        _intersection_domain._type = type;
      }

      template<typename Context, typename IntersectionType>
      void unbind(
        Context& ctx,
        IntersectionType type,
        const typename IntersectionDomain::Intersection& intersection,
        typename IntersectionDomain::Index index,
        const typename IntersectionDomain::Cell&,
        typename IntersectionDomain::Index,
        typename IntersectionDomain::Index
        )
      {
        _intersection_domain._intersection = nullptr;
        _intersection_domain._index = IntersectionDomain::invalid_index;
        _intersection_domain._type = IntersectionDomain::invalid_type;
      }

    private:

      IntersectionDomain _intersection_domain;

    };

    template<typename ES>
    auto intersectionDomainData(const ES&)
    {
      return IntersectionDomainData<ES>{};
    }


    template<typename ES, typename... Components>
    class InsideCell
      : public Components...
    {

      struct InsideTraits
      {
        using EntitySet = ES;
        using Element   = typename ES::template Codim<0>::Entity;
        using Index     = typename ES::IndexSet::Index;
      };

    public:

      using Inside    = InsideCell;

      InsideCell(Components&&... components)
        : Components(std::forward<Components>(components))...
      {}

      Inside& inside()
      {
        return *this;
      }

      const Inside& inside() const
      {
        return *this;
      }

      template<typename Context, typename Engine>
      void setup(Context& ctx, Engine& engine)
      {
        auto setup = [&](auto&& c)
          -> decltype(c.setup(ctx,*this,engine))
          { return c.setup(ctx,*this,engine); };
        applyToVariadicArguments{invoke_if_possible_discard_return(setup,static_cast<Components&>(*this))...};
      }

      template<typename Context>
      void bind(Context& ctx, const typename InsideTraits::Element& element, typename InsideTraits::Index entity_index, typename InsideTraits::Index unique_index)
      {
        auto bind = [&](auto&& c)
          -> decltype(c.bind(ctx,*this,element,entity_index,unique_index))
          { return c.bind(ctx,*this,element,entity_index,unique_index); };
        applyToVariadicArguments{invoke_if_possible_discard_return(bind,static_cast<Components&>(*this))...};
      }

      template<typename Context>
      void unbind(Context& ctx, const typename InsideTraits::Element& element, typename InsideTraits::Index entity_index, typename InsideTraits::Index unique_index)
      {
        auto unbind = [&](auto&& c)
          -> decltype(c.unbind(ctx,*this,element,entity_index,unique_index))
          { return c.unbind(ctx,*this,element,entity_index,unique_index); };
        auto apply = [&](auto&& c) { return invoke_if_possible_discard_return(unbind,std::forward<decltype(c)>(c)); };
        applyToVariadicArgumentsWithOrder(apply,std::forward_as_tuple(static_cast<Components&>(*this)...),reverse_index_sequence_for<Components...>{});
      }

    };

    template<typename EntitySet, typename... Components>
    auto insideCell(const EntitySet&, Components&&... components)
    {
      return InsideCell<EntitySet,std::decay_t<Components>...>{std::forward<Components>(components)...};
    }


    template<typename ES, typename... Components>
    class OutsideCell
    {

      struct OutsideTraits
      {
        using EntitySet    = ES;
        using Element      = typename ES::template Codim<0>::Entity;
        using Intersection = typename ES::Intersection;
        using Index        = typename ES::IndexSet::Index;
      };

    private:

      struct Outside
        : public Components...
      {

        Outside(Components&&... components)
          : Components(std::forward<Components>(components))...
        {}

      };

      Outside _outside;

      template<typename Context>
      void doBind(Context& ctx, const typename OutsideTraits::Element& outside_element, typename OutsideTraits::Index outside_entity_index, typename OutsideTraits::Index outside_unique_index)
      {
        auto bind = [&](auto&& c)
          -> decltype(c.bind(ctx,_outside,outside_element,outside_entity_index,outside_unique_index))
          { return c.bind(ctx,_outside,outside_element,outside_entity_index,outside_unique_index); };
        applyToVariadicArguments{invoke_if_possible_discard_return(bind,static_cast<Components&>(_outside))...};
      }

      template<typename Context>
      void doUnbind(Context& ctx, const typename OutsideTraits::Element& outside_element, typename OutsideTraits::Index outside_entity_index, typename OutsideTraits::Index outside_unique_index)
      {
        auto unbind = [&](auto&& c)
          -> decltype(c.unbind(ctx,_outside,outside_element,outside_entity_index,outside_unique_index))
          { return c.unbind(ctx,_outside,outside_element,outside_entity_index,outside_unique_index); };
        auto apply = [&](auto&& c) { return invoke_if_possible_discard_return(unbind,std::forward<decltype(c)>(c)); };
        applyToVariadicArgumentsWithOrder(apply,std::forward_as_tuple(static_cast<Components&>(_outside)...),reverse_index_sequence_for<Components...>{});
      }

    public:

      OutsideCell(Components&&... components)
        : _outside(std::forward<Components>(components)...)
      {}

      Outside& outside()
      {
        return _outside;
      }

      const Outside& outside() const
      {
        return _outside;
      }

      template<typename Context, typename Engine>
      void setup(Context& ctx, Engine& engine)
      {
        auto setup = [&](auto&& c)
          -> decltype(c.setup(ctx,_outside,engine))
          { return c.setup(ctx,_outside,engine); };
        applyToVariadicArguments{invoke_if_possible_discard_return(setup,static_cast<Components&>(_outside))...};
      }

      template<typename Context>
      void bind(
        Context& ctx,
        std::integral_constant<IntersectionType,IntersectionType::skeleton>,
        const typename OutsideTraits::Intersection&, typename OutsideTraits::Index,
        const typename OutsideTraits::Element& outside_element, typename OutsideTraits::Index outside_entity_index, typename OutsideTraits::Index outside_unique_index
        )
      {
        doBind(ctx,outside_element,outside_entity_index,outside_unique_index);
      }

      template<typename Context>
      void bind(
        Context& ctx,
        std::integral_constant<IntersectionType,IntersectionType::periodic>,
        const typename OutsideTraits::Intersection&, typename OutsideTraits::Index,
        const typename OutsideTraits::Element& outside_element, typename OutsideTraits::Index outside_entity_index, typename OutsideTraits::Index outside_unique_index
        )
      {
        doBind(ctx,outside_element,outside_entity_index,outside_unique_index);
      }

      template<typename Context>
      void unbind(
        Context& ctx,
        std::integral_constant<IntersectionType,IntersectionType::skeleton>,
        const typename OutsideTraits::Intersection&, typename OutsideTraits::Index,
        const typename OutsideTraits::Element& outside_element, typename OutsideTraits::Index outside_entity_index, typename OutsideTraits::Index outside_unique_index
        )
      {
        doUnbind(ctx,outside_element,outside_entity_index,outside_unique_index);
      }

      template<typename Context>
      void unbind(
        Context& ctx,
        std::integral_constant<IntersectionType,IntersectionType::periodic>,
        const typename OutsideTraits::Intersection&, typename OutsideTraits::Index,
        const typename OutsideTraits::Element& outside_element, typename OutsideTraits::Index outside_entity_index, typename OutsideTraits::Index outside_unique_index
        )
      {
        doUnbind(ctx,outside_element,outside_entity_index,outside_unique_index);
      }

    };

    template<typename EntitySet, typename... Components>
    auto outsideCell(const EntitySet&, Components&&... components)
    {
      return OutsideCell<EntitySet,std::decay_t<Components>...>{std::forward<Components>(components)...};
    }

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ASSEMBLER_CELLDATA_HH
