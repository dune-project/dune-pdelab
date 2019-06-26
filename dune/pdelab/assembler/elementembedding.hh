// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_ASSEMBLER_ELEMENTEMBEDDING_HH
#define DUNE_PDELAB_ASSEMBLER_ELEMENTEMBEDDING_HH

#include <cstddef>

#include <dune/common/reservedvector.hh>

#include <dune/geometry/identitygeometry.hh>

namespace Dune::PDELab::Experimental {

  template<typename Geometry>
  class ElementEmbedding
  {

  public:

    using size_type                        = std::size_t;
    using Global                           = Geometry;
    using Field                            = typename Geometry::ctype;
    using Element                          = IdentityGeometry<Field,Geometry::mydimension>;
    using Local                            = Element;
    using Inside                           = Element;
    using LocalCoordinate                  = typename Geometry::LocalCoordinate;
    using ElementCoordinate                = LocalCoordinate;
    using GlobalCoordinate                 = typename Geometry::GlobalCoordinate;
    using JacobianTransposed               = typename Geometry::JacobianTransposed;
    using JacobianInverseTransposed        = typename Geometry::JacobianInverseTransposed;
    using InsideJacobianTransposed         = JacobianTransposed;
    using InsideJacobianInverseTransposed  = JacobianInverseTransposed;
    using OutsideJacobianTransposed        = int;
    using OutsideJacobianInverseTransposed = int;
    using EmbeddingDescriptor              = ReservedVector<
      LocalCoordinate,
      (1 << (Geometry::mydimension-1))
      >;

    static constexpr int dimLocal = Geometry::mydimension;
    static constexpr int dimWorld = Geometry::coorddimension;

    const Global& global() const
    {
      return *_global;
    }

    Element local() const
    {
      return Local(global().type());
    }

    Element inside() const
    {
      return local();
    }

    LocalCoordinate inside(const LocalCoordinate& local) const
    {
      return local;
    }

    template<typename P>
    InsideJacobianTransposed jacobianTransposed(const P& p) const
    {
      return global().jacobianTransposed(p.inside());
    }

    template<typename P>
    InsideJacobianInverseTransposed jacobianInverseTransposed(const P& p) const
    {
      return global().jacobianInverseTransposed(p.inside());
    }

    template<typename P>
    InsideJacobianTransposed insideJacobianTransposed(const P& p) const
    {
      return global().jacobianTransposed(p.inside());
    }

    template<typename P>
    InsideJacobianInverseTransposed insideJacobianInverseTransposed(const P& p) const
    {
      return global().jacobianInverseTransposed(p.inside());
    }

    template<typename P>
    Field integrationElement(const P& p) const
    {
      return global().integrationElement(p.local());
    }

    ElementEmbedding(const Geometry& geo)
      : _global(&geo)
    {}

    EmbeddingDescriptor insideDescriptor() const
    {
      return {};
    }

  private:

    const Global* _global;

  };

} // namespace Dune::PDELab::Experimental

#endif // DUNE_PDELAB_ASSEMBLER_ELEMENTEMBEDDING_HH
