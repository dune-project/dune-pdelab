#ifndef DUNE_PDELAB_COMMON_REFERENCEELEMENTS_HH
#define DUNE_PDELAB_COMMON_REFERENCEELEMENTS_HH

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>

namespace Dune {
  namespace PDELab {

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

  // inject the function into the Dune namespace to enable ADL.
  // TODO: Remove this injection once the core modules gain their own version of this function!
  using Dune::PDELab::referenceElement;

} // namespace Dune

#endif // DUNE_PDELAB_COMMON_REFERENCEELEMENTS_HH
