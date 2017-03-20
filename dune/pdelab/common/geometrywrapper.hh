// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_COMMON_GEOMETRYWRAPPER_HH
#define DUNE_PDELAB_COMMON_GEOMETRYWRAPPER_HH

#include <dune/common/fvector.hh>

namespace Dune {
  namespace PDELab {

    //! Wrap element
    /**
     * \todo Please doc me!
     */
    template<typename E>
    class ElementGeometry
    {
    public:
      //! \todo Please doc me!
      typedef typename E::Geometry Geometry;
      //! \todo Please doc me!
      typedef E Entity;

      //! \todo Please doc me!
      ElementGeometry (const E& e_)
        : e(e_)
      {}

      //! \todo Please doc me!
      Geometry geometry () const
      {
        return e.geometry();
      }

      //! \todo Please doc me!
      const Entity& entity () const
      {
        return e;
      }

      //! \todo Please doc me!
      const Entity& hostEntity () const
      {
        return e;
      }

    private:
      const E& e;
    };


    //! Wrap intersection
    /**
     * \todo Please doc me!
     */
    template<typename I>
    class IntersectionGeometry
    {
    private:
      const I& i;
      const unsigned int index;
    public:
      //! \todo Please doc me!
      typedef typename I::Geometry Geometry;
      //! \todo Please doc me!
      typedef typename I::LocalGeometry LocalGeometry;
      //! \todo Please doc me!
      typedef typename I::Entity Entity;
      //! \todo Please doc me!
      typedef typename Geometry::ctype ctype;
      //! \todo Please doc me!
      enum DUNE_DEPRECATED_MSG("Deprecated: get this dimension from the grid itself, or from an element!") { dimension=Entity::dimension };

      /** \brief Dimension of the domain space of the geometry */
      enum { mydimension=I::mydimension };

      /** \brief Dimension of the image space of the geometry */
      enum { coorddimension=Geometry::coorddimension };

      //! \todo Please doc me!
      IntersectionGeometry (const I& i_, unsigned int index_)
        : i(i_), index(index_)
      {}

      //! \todo Please doc me!
      int insideDomainIndex() const
      {
        return 0;
      }

      //! \todo Please doc me!
      int outsideDomainIndex() const
      {
        const bool is_boundary = i.boundary();
        return 0 - int(is_boundary);
      }

      //! return true if intersection is with interior or exterior boundary (see the cases above)
      bool boundary () const
      {
        return i.boundary();
      }

      //! @brief return true if intersection is shared with another element.
      bool neighbor () const
      {
        return i.neighbor();
      }

      /*! @brief geometrical information about this intersection in local
        coordinates of the inside() entity.

        This method returns a Geometry object that provides a mapping from
        local coordinates of the intersection to local coordinates of the
        inside() entity.
      */
      LocalGeometry geometryInInside () const
      {
        return i.geometryInInside();
      }

      /*! @brief geometrical information about this intersection in local
        coordinates of the outside() entity.

        This method returns a Geometry object that provides a mapping from
        local coordinates of the intersection to local coordinates of the
        outside() entity.
      */
      LocalGeometry geometryInOutside () const
      {
        return i.geometryInOutside();
      }

      /*! @brief geometrical information about this intersection in global coordinates.

        This method returns a Geometry object that provides a mapping from
        local coordinates of the intersection to global (world) coordinates.
      */
      Geometry geometry () const
      {
        return i.geometry();
      }

      //! Local number of codim 1 entity in the inside() Entity where intersection is contained in
      int indexInInside () const
      {
        return i.indexInInside ();
      }

      //! Local number of codim 1 entity in outside() Entity where intersection is contained in
      int indexInOutside () const
      {
        return i.indexInOutside ();
      }

      /*! @brief Return an outer normal (length not necessarily 1)

        The returned vector may depend on local position within the intersection.
      */
      Dune::FieldVector<ctype, coorddimension> outerNormal (const Dune::FieldVector<ctype, mydimension>& local) const
      {
        return i.outerNormal(local);
      }

      /*! @brief return outer normal scaled with the integration element
        @copydoc outerNormal
        The normal is scaled with the integration element of the intersection. This
        method is redundant but it may be more efficent to use this function
        rather than computing the integration element via intersectionGlobal().
      */
      Dune::FieldVector<ctype, coorddimension> integrationOuterNormal (const Dune::FieldVector<ctype, mydimension>& local) const
      {
        return i.integrationOuterNormal(local);
      }

      /*! @brief Return unit outer normal (length == 1)

        The returned vector may depend on the local position within the intersection.
        It is scaled to have unit length.
      */
      Dune::FieldVector<ctype, coorddimension> unitOuterNormal (const Dune::FieldVector<ctype, mydimension>& local) const
      {
        return i.unitOuterNormal(local);
      }

      /*! @brief Return unit outer normal (length == 1)

        The returned vector may depend on the local position within the intersection.
        It is scaled to have unit length.
      */
      Dune::FieldVector<ctype, coorddimension> centerUnitOuterNormal () const
      {
        return i.centerUnitOuterNormal();
      }

      /*! @brief return Entity on the inside of this intersection.
          That is the Entity where we started this.
      */
      Entity inside() const
      {
        return i.inside();
      }

      /*! @brief return Entity on the inside of this intersection.
          That is the Entity where we started this.
      */
      Entity insideHostEntity() const
      {
        DUNE_THROW(Dune::Exception,"This should never be called.");
        return i.inside();
      }

      /*! @brief return Entity on the outside of this intersection.
          That is the neighboring Entity.

        @warning Don't call this method if there is no neighboring Entity
        (neighbor() returns false). In this case the result is undefined.
      */
      Entity outside() const
      {
        return i.outside();
      }

      //! \todo Please doc me!
      const I& intersection () const
      {
        return i;
      }

      unsigned int intersectionIndex() const
      {
        return index;
      }
    };

  }
}

#endif // DUNE_PDELAB_COMMON_GEOMETRYWRAPPER_HH
