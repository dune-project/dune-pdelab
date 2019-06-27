// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_ASSEMBLER_ELEMENTDATA_HH
#define DUNE_PDELAB_ASSEMBLER_ELEMENTDATA_HH

#include <cassert>
#include <type_traits>
#include <optional>

#include <dune/common/reservedvector.hh>

#include <dune/geometry/identitygeometry.hh>

#include <dune/typetree/childextraction.hh>

#include <dune/pdelab/common/intersectiontype.hh>
#include <dune/pdelab/common/quadraturerules.hh>
#include <dune/pdelab/gridfunctionspace/flavor.hh>
#include <dune/pdelab/assembler/utility.hh>
#include <dune/pdelab/assembler/context.hh>
#include <dune/pdelab/assembler/finiteelementwrapper.hh>
#include <dune/pdelab/assembler/quadraturerule.hh>

namespace Dune::PDELab::Experimental {

  namespace ElementFlavor {

    template<typename Context>
    using TestLocalSpace  = typename Context::Engine::template TestLocalSpace<Context>;

    template<typename Context>
    using TestSpaceCache  = typename Context::Engine::template TestSpaceCache<Context>;

    template<typename Context>
    using TrialLocalSpace = typename Context::Engine::template TrialLocalSpace<Context>;

    template<typename Context>
    using TrialSpaceCache = typename Context::Engine::template TrialSpaceCache<Context>;

    template<bool enable_flavors>
    struct Inside
    {

      template<typename QP>
      static constexpr auto& quadratureCoordinate(const QP& qp)
      {
        return qp.inside();
      }

      template<typename QP>
      static constexpr auto& elementJacobianTransposed(const QP& qp)
      {
        return qp.insideJacobianTransposed();
      }

      template<typename QP>
      static constexpr auto& elementJacobianInverseTransposed(const QP& qp)
      {
        return qp.insideJacobianInverseTransposed();
      }

      template<typename Embedding>
      static constexpr auto embeddingDescriptor(const Embedding& embedding)
      {
        return embedding.insideDescriptor();
      }

      using Test  = std::conditional_t<enable_flavors,Flavor::InsideTest,Flavor::Generic>;
      using Trial = std::conditional_t<enable_flavors,Flavor::InsideTrial,Flavor::Generic>;

      struct Flavor {
        using Test  = Inside::Test;
        using Trial = Inside::Trial;
      };

    };

    template<bool enable_flavors>
    struct Outside
    {
      template<typename QP>
      static constexpr auto& quadratureCoordinate(const QP& qp)
      {
        return qp.outside();
      }

      template<typename QP>
      static constexpr auto& elementJacobianTransposed(const QP& qp)
      {
        return qp.outsideJacobianTransposed();
      }

      template<typename QP>
      static constexpr auto& elementJacobianInverseTransposed(const QP& qp)
      {
        return qp.outsideJacobianInverseTransposed();
      }

      template<typename Embedding>
      static constexpr auto embeddingDescriptor(const Embedding& embedding)
      {
        return embedding.outsideDescriptor();
      }

      using Test  = std::conditional_t<enable_flavors,Flavor::OutsideTest,Flavor::Generic>;
      using Trial = std::conditional_t<enable_flavors,Flavor::OutsideTrial,Flavor::Generic>;

      struct Flavor {
        using Test  = Outside::Test;
        using Trial = Outside::Trial;
      };

    };

  }

  template<typename Context>
  class ElementGridData
    : public Context
  {

  public:

    using EntitySet = typename Context::EntitySet;
    using Entity    = typename EntitySet::template Codim<0>::Entity;
    using Index     = typename EntitySet::IndexSet::Index;
    using Geometry  = typename Entity::Geometry;
    using Embedding = ElementEmbedding<Geometry>;
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

    ElementGridData(Context&& ctx)
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
  auto elementGridData(Context&& ctx)
  {
    return ElementGridData<Context>{std::move(ctx)};
  }


  template<typename Context>
  class ElementDomainData
    : public Context
  {

  public:

    class ElementDomain
    {

      friend class ElementDomainData;

    public:

      using EntitySet = typename Context::EntitySet;
      using Entity    = typename EntitySet::template Codim<0>::Entity;
      using Field     = typename EntitySet::ctype;

    private:

      using Index = typename EntitySet::IndexSet::Index;

    public:

      // re-use the element embedding
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
        return QuadratureRule(_data,rule,embedding());
      }

      ElementDomain(ElementDomainData& data)
        : _data(data)
      {}

    private:

      ElementDomainData& _data;

    };

    ElementDomain elementDomain()
    {
      return {*this};
    }

    ElementDomainData(Context&& ctx)
      : Context(std::move(ctx))
    {}

  };

  template<typename Context>
  auto elementDomainData(Context&& ctx)
  {
    return ElementDomainData<Context>{std::move(ctx)};
  }


  template<typename Context>
  class IntersectionDomainData
    : public Context
  {

  public:

    // TODO: The upper bound here is kind of a hack!
    using EmbeddingDescriptor = Dune::ReservedVector<
      typename Context::EntitySet::Intersection::Geometry::GlobalCoordinate,
      (1 << Context::EntitySet::Intersection::Geometry::mydimension)
      >;

    class IntersectionDomain
    {

      friend class IntersectionDomainData;

    public:

      using EntitySet    = typename Context::EntitySet;
      using Entity       = typename EntitySet::template Codim<0>::Entity;
      using Intersection = typename EntitySet::Intersection;
      using Field        = typename EntitySet::Field;
      using Index        = typename EntitySet::IndexSet::Index;

    private:

      static constexpr Index invalid_index = EntitySet::IndexSet::invalidIndex();
      static constexpr auto invalid_type   = IntersectionType::invalid;

    public:

      struct Embedding {

        friend class IntersectionDomain;

        using size_type                        = std::size_t;
        using Field                            = typename IntersectionDomain::Field;
        using Intersection                     = typename Context::EntitySet::Intersection;
        using Global                           = typename Intersection::Geometry;
        using Element                          = typename Intersection::LocalGeometry;
        using Inside                           = typename Intersection::LocalGeometry;
        using Outside                          = typename Intersection::LocalGeometry;
        using LocalCoordinate                  = typename Global::LocalCoordinate;
        using GlobalCoordinate                 = typename Global::GlobalCoordinate;
        using ElementCoordinate                = typename Inside::GlobalCoordinate;
        using JacobianTransposed               = typename Global::JacobianTransposed;
        using JacobianInverseTransposed        = typename Global::JacobianInverseTransposed;
        using InsideJacobianTransposed         = typename Context::Inside::Embedding::JacobianTransposed;
        using InsideJacobianInverseTransposed  = typename Context::Inside::Embedding::JacobianInverseTransposed;
        using OutsideJacobianTransposed        = typename Context::Outside::Embedding::JacobianTransposed;
        using OutsideJacobianInverseTransposed = typename Context::Outside::Embedding::JacobianInverseTransposed;
        using EmbeddingDescriptor              = typename IntersectionDomainData::EmbeddingDescriptor;

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

        template<typename P>
        GlobalCoordinate unitOuterNormal(const P& p) const
        {
          return _data.intersection().unitOuterNormal(p.local());
        }

        template<typename P>
        GlobalCoordinate unitOuterNormal(const LocalCoordinate& local) const
        {
          return _data.intersection().unitOuterNormal(local);
        }

        template<typename P>
        InsideJacobianTransposed insideJacobianTransposed(const P& p) const
        {
          return _data.inside().embedding().global().jacobianTransposed(p.inside());
        }

        template<typename P>
        InsideJacobianInverseTransposed insideJacobianInverseTransposed(const P& p) const
        {
          return _data.inside().embedding().global().jacobianInverseTransposed(p.inside());
        }

        template<typename P>
        OutsideJacobianTransposed outsideJacobianTransposed(const P& p) const
        {
          return _data.outside().embedding().global().jacobianTransposed(p.outside());
        }

        template<typename P>
        OutsideJacobianInverseTransposed outsideJacobianInverseTransposed(const P& p) const
        {
          return _data.outside().embedding().global().jacobianInverseTransposed(p.outside());
        }

        EmbeddingDescriptor insideDescriptor() const
        {
          return _data.insideDescriptor();
        }

        EmbeddingDescriptor outsideDescriptor() const
        {
          return _data.outsideDescriptor();
        }

        size_type insideIndex() const
        {
          return _data.intersection().indexInInside();
        }

        size_type outsideIndex() const
        {
          return _data.intersection().indexInOutside();
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

      typename Embedding::GlobalCoordinate unitOuterNormal(const typename Embedding::LocalCoordinate& local) const
      {
        return intersection().unitOuterNormal(local);
      }

      bool bound() const
      {
        return _data._intersection;
      }

      auto quadratureRule(std::size_t order, QuadratureType::Enum quadrature_type = QuadratureType::GaussLegendre) const
      {
        auto& rule = QuadratureRules<typename Embedding::Global::ctype,Embedding::Global::mydimension>::rule(embedding().global().type(),order,quadrature_type);
        return QuadratureRule(_data,rule,embedding());
      }

      IntersectionDomain(IntersectionDomainData& data)
        : _data(data)
      {}

    private:

      IntersectionDomainData& _data;

    };

    using Domain = IntersectionDomain;

    IntersectionDomain intersectionDomain()
    {
      return {*this};
    }

    Domain domain()
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
      _inside_descriptor.reset();
      _outside_descriptor.reset();
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

    const EmbeddingDescriptor& insideDescriptor() const
    {
      assert(_intersection);
      if (not _inside_descriptor)
      {
        auto& geo = intersectionGeometryInInside();
        _inside_descriptor.emplace();
        for (auto i = 0, corners = geo.corners() ; i < corners ; ++i)
          _inside_descriptor->push_back(geo.corner(i));
      }
      return *_inside_descriptor;
    }

    const EmbeddingDescriptor& outsideDescriptor() const
    {
      assert(_intersection);
      if (not _outside_descriptor)
      {
        auto& geo = intersectionGeometryInOutside();
        _outside_descriptor.emplace();
        for (auto i = 0, corners = geo.corners() ; i < corners ; ++i)
          _outside_descriptor->push_back(geo.corner(i));
      }
      return *_outside_descriptor;
    }

  private:

    const typename IntersectionDomain::Intersection* _intersection;
    typename IntersectionDomain::Index _intersection_index;
    IntersectionType _intersection_type;
    mutable std::optional<typename IntersectionDomain::Embedding::Global> _geometry;
    mutable std::optional<typename IntersectionDomain::Embedding::Inside> _geometry_in_inside;
    mutable std::optional<typename IntersectionDomain::Embedding::Outside> _geometry_in_outside;
    mutable std::optional<EmbeddingDescriptor> _inside_descriptor;
    mutable std::optional<EmbeddingDescriptor> _outside_descriptor;

  };

  template<typename Context>
  auto intersectionDomainData(Context&& ctx)
  {
    return IntersectionDomainData<Context>{std::move(ctx)};
  }


  template<typename Context>
  class InsideElement
    : public Context
  {

    struct InsideTraits
    {
      using EntitySet = typename Context::EntitySet;
      using Element   = typename EntitySet::template Codim<0>::Entity;
      using Index     = typename EntitySet::IndexSet::Index;
    };

  public:

    using Inside    = InsideElement;
    using Element   = InsideElement;

    InsideElement(Context&& ctx)
      : Context(std::move(ctx))
    {}

    Element& element()
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
  auto insideElement(Context&& ctx)
  {
    return InsideElement<Context>{std::move(ctx)};
  }


  template<typename Context, typename ElementContext>
  class OutsideElement
    : public Context
  {

    struct OutsideTraits
    {
      using EntitySet    = typename ElementContext::EntitySet;
      using Entity       = typename EntitySet::template Codim<0>::Entity;
      using Intersection = typename EntitySet::Intersection;
      using Index        = typename EntitySet::IndexSet::Index;
    };

  public:

    using Outside = ElementContext;

  private:

    Outside _outside;

  public:

    OutsideElement(Context&& ctx, ElementContext&& element_ctx)
      : Context(std::move(ctx))
      , _outside(std::move(element_ctx))
    {}

    Outside& outside()
    {
      return _outside;
    }

    const Outside& outside() const
    {
      return _outside;
    }

    void setup()
    {
      Context::setup();
      outside().setup();
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

    template<typename QuadratureRule>
    std::enable_if_t<
      models<Concept::IntersectionEmbedding,typename QuadratureRule::Embedding>()
      >
    beginQuadrature(QuadratureRule& qr)
    {
      Context::beginQuadrature(qr);
      if (_outside.bound())
        _outside.beginQuadrature(qr);
    }

    template<typename QuadratureRule>
    std::enable_if_t<
      models<Concept::IntersectionEmbedding,typename QuadratureRule::Embedding>()
      >
    endQuadrature(QuadratureRule& qr)
    {
      if (_outside.bound())
       _outside.endQuadrature(qr);
      Context::endQuadrature(qr);
    }

    template<typename QuadratureRule>
    std::enable_if_t<
      not models<Concept::IntersectionEmbedding,typename QuadratureRule::Embedding>()
      >
    beginQuadrature(QuadratureRule& qr)
    {
      Context::beginQuadrature(qr);
    }

    template<typename QuadratureRule>
    std::enable_if_t<
      not models<Concept::IntersectionEmbedding,typename QuadratureRule::Embedding>()
      >
    endQuadrature(QuadratureRule& qr)
    {
      Context::endQuadrature(qr);
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
      Dune::PDELab::Experimental::Context::bind(_outside,outside_element,outside_entity_index,outside_unique_index);
    }

    void doBind(
      std::integral_constant<IntersectionType,IntersectionType::periodic>,
      const typename OutsideTraits::Entity& outside_element,
      typename OutsideTraits::Index outside_entity_index,
      typename OutsideTraits::Index outside_unique_index
      )
    {
      Dune::PDELab::Experimental::Context::bind(_outside,outside_element,outside_entity_index,outside_unique_index);
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
      Dune::PDELab::Experimental::Context::unbind(_outside,outside_element,outside_entity_index,outside_unique_index);
    }

    void doUnbind(
      std::integral_constant<IntersectionType,IntersectionType::periodic>,
      const typename OutsideTraits::Entity& outside_element,
      typename OutsideTraits::Index outside_entity_index,
      typename OutsideTraits::Index outside_unique_index
      )
    {
      Dune::PDELab::Experimental::Context::unbind(_outside,outside_element,outside_entity_index,outside_unique_index);
    }

  };

  template<typename ElementContext, typename Context>
  auto outsideElement(ElementContext&& element_ctx, Context&& ctx)
  {
    return OutsideElement<Context,ElementContext>{std::move(ctx),std::move(element_ctx)};
  }

} // namespace Dune::PDELab::Experimental

#endif // DUNE_PDELAB_ASSEMBLER_ELEMENTDATA_HH
