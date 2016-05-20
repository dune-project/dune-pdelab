#ifndef DUNE_PDELAB_COMMON_QUADRATURERULES_HH
#define DUNE_PDELAB_COMMON_QUADRATURERULES_HH

#include <dune/geometry/quadraturerules.hh>

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

  } // namespace PDELab

  // inject the function into the Dune namespace to enable ADL.
  // TODO: Remove this injection once the core modules gain their own version of this function!
  using Dune::PDELab::quadratureRule;

} // namespace Dune

#endif // DUNE_PDELAB_COMMON_QUADRATURERULES_HH
