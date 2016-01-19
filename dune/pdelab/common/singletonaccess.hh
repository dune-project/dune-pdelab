#ifndef DUNE_PDELAB_COMMON_SINGLETONACCESS_HH
#define DUNE_PDELAB_COMMON_SINGLETONACCESS_HH

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>

namespace Dune {
  namespace PDELab {

    //! Wrapper for Dune::QuadratureRule with value semantics.
    /**
     * This class wraps a Dune::QuadratureRule and exposes the relevant parts
     * of its interface (iteration, size, order and geometry information).
     * It does, however, not leak internal information by inheriting from
     * std::vector, so you should not rely on the additional types and methods
     * inherited from std::vector.
     *
     * In contrast to Dune::QuadratureRule, QuadratureRuleWrapper should be
     * used with value semantics, i.e. you should store it by value (not by
     * reference) and create copies as needed. Copies are very cheap, as the
     * class internally only stores a reference to the underlying QuadratureRule.
     *
     * \note Users will normally not construct a QuadratureRuleWrapper directly, but
     *       obtain the rule by calling Dune::PDELab::quadratureRule().
     */
    template<typename QR>
    class QuadratureRuleWrapper
    {

    public:

      //! The coordinate type of the local coordinates of the rule.
      using CoordType = typename QR::CoordType;

      //! The size type used by the container.
      using size_type = typename QR::size_type;

      //! A const iterator over the quadrature points.
      using const_iterator = typename QR::const_iterator;

      //! An iterator over the quadrature points (always const, as the container is read-only).
      using iterator = const_iterator;

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
        return _quadrature_rule->begin();
      }

      //! Returns an iterator pointing after the last quadrature point.
      const_iterator end() const
      {
        return _quadrature_rule->end();
      }

#ifndef DOXYGEN

      QuadratureRuleWrapper(const QR& quadrature_rule)
        : _quadrature_rule(&quadrature_rule)
      {}

#endif

    private:

      const QR* _quadrature_rule;

    };

    //! Returns a quadrature rule for the given geometry.
    /**
     * This function returns a quadrature rule for the given geometry `geo` which can
     * be used for exact integration of polynomial functions up to the order `order`.
     * Optionally, the type of quadrature rule can be passed as a third parameter. It
     * defaults to Gauss-Legendre quadrature.
     *
     * This function should be used unqualified, e.g.
     *
     * ```c++
     * auto rule = quadratureRule(geo,4);
     * ```
     *
     * \warning Do *not* add `Dune::PDELab` in front of it, as that might cause compiler
     * errors with future versions of the core modules.
     *
     * The function wraps the quadrature rule in a QuadratureRuleWrapper to provide value
     * semantics.
     *
     * \param geo             The geometry for which to obtain a quadrature rule.
     * \param order           The maximum polynomial order up to which the rule should be exact.
     * \param quadrature_type The type of quadrature (default: Gauss-Legendre).
     */
    template<typename Geometry>
    QuadratureRuleWrapper<
      QuadratureRule<
        typename Geometry::ctype,
        Geometry::mydimension
        >
      >
    quadratureRule(const Geometry& geo, std::size_t order, QuadratureType::Enum quadrature_type = QuadratureType::GaussLegendre)
    {
      return { QuadratureRules<typename Geometry::ctype,Geometry::mydimension>::rule(geo.type(),order,quadrature_type) };
    }



    //! Wrapper for Dune::ReferenceElement with value semantics.
    /**
     * This class wraps a Dune::ReferenceElement and exposes the relevant parts
     * of its interface. Note, however, that its signature is different, so code
     * that relies on partial specialization of the reference element type will
     * not work with the wrapper.
     *
     * For the exact documentation of the individual members, please see the documentation
     * of Dune::ReferenceElement.
     *
     * This class exports three additional items of static information:
     *
     * - `CoordinateField` denotes the field type of the spatial coordinates of the
     *   reference element.
     * - `Coordinate` denotes the type of spatial coordinate vectors.
     * - `dimension` denotes the dimension of the element.
     *
     * \note Users will normally not construct a ReferenceWrapper directly, but
     *       obtain the rule by calling Dune::PDELab::referenceElement().
     *
     * \sa Dune::ReferenceElement.
     *
     */
    template<typename RE>
    class ReferenceElementWrapper
    {

    public:

      template<int codim>
      using Codim = typename RE::template Codim<codim>;

      //! The coordinate field type.
      using ctype = typename Codim<0>::Geometry::ctype;
      //! The coordinate field type.
      using CoordinateField = ctype;

      //! The dimension of the reference element.
      static const std::size_t dimension = Codim<0>::Geometry::coorddimension;

      //! The coordinate type of the reference element.
      using Coordinate = FieldVector<CoordinateField,dimension>;

      int size(int c) const
      {
        return _ref_el->size(c);
      }

      int size(int i, int c, int cc) const
      {
        return _ref_el->size(i,c,cc);
      }

      int subEntity(int i, int c, int ii, int cc) const
      {
        return _ref_el->subEntity(i,c,ii,cc);
      }

      GeometryType type(int i, int c) const
      {
        return _ref_el->type(i,c);
      }

      GeometryType type() const
      {
        return _ref_el->type();
      }

      const Coordinate& position(int i, int c) const
      {
        return _ref_el->position(i,c);
      }

      bool checkInside(const Coordinate& local) const
      {
        return _ref_el->checkInside(local);
      }

      template<int codim>
      typename Codim<codim>::Geometry geometry(int i) const
      {
        return _ref_el->geometry(i);
      }

      CoordinateField volume() const
      {
        return _ref_el->volume();
      }

      const Coordinate& integrationOuterNormal(int face) const
      {
        return _ref_el->integrationOuterNormal(face);
      }

#ifndef DOXYGEN

      ReferenceElementWrapper(const RE& ref_el)
        : _ref_el(&ref_el)
      {}

#endif // DOXYGEN

    private:

      const RE* _ref_el;

    };


    //! Returns the reference element for the given geometry.
    /**
     * This function returns the reference element for the given geometry `geo`.
     *
     * This function should be used unqualified, e.g.
     *
     * ```c++
     * auto rule = quadratureRule(geo,4);
     * ```
     *
     * \warning Do *not* add `Dune::PDELab` in front of it, as that might cause compiler
     * errors with future versions of the core modules.
     *
     * The function wraps the reference element in a ReferenceElementWrapper to provide value
     * semantics.
     *
     * \param geo             The geometry for which to obtain the reference element.
     */
    template<typename Geometry>
    ReferenceElementWrapper<
      ReferenceElement<
        typename Geometry::ctype,
        Geometry::mydimension
        >
      >
    referenceElement(const Geometry& geo)
    {
      return { ReferenceElements<typename Geometry::ctype,Geometry::mydimension>::general(geo.type()) };
    }

  } // namespace PDELab

  // inject the functions into the Dune namespace to enable ADL.
  // TODO: Remove this injection once the core modules gain their own versions of these functions!
  using Dune::PDELab::quadratureRule;
  using Dune::PDELab::referenceElement;

} // namespace Dune

#endif // DUNE_PDELAB_COMMON_SINGLETONACCESS_HH
