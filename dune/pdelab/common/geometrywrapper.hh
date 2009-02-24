// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_GEOMETRYWRAPPER_HH
#define DUNE_PDELAB_GEOMETRYWRAPPER_HH

namespace Dune {
  namespace PDELab {

	template<typename E>
	class ElementGeometry
	{
	public:
	  typedef typename E::Geometry Geometry;
	  typedef E Entity;

	  ElementGeometry (const E& e_)
		: e(e_)
	  {}

	  const Geometry& geometry () const
	  {
		return e.geometry();
	  }

	  const Entity& entity () const
	  {
		return e;
	  }
  
	private:
	  const E& e;
	};


	template<typename I>
	class IntersectionGeometry
	{
	public:
	  typedef typename I::Geometry Geometry;
	  typedef typename I::LocalGeometry LocalGeometry;
	  typedef typename I::Entity Entity;
	  typedef typename I::EntityPointer EntityPointer;
	  typedef typename Entity::ctype ctype;
	  enum { dimension=Entity::dimension };
	  enum { dimensionworld=Entity::dimensionworld };

	  IntersectionGeometry (const I& i_)
		: i(i_)
	  {}


	  //! return true if intersection is with interior or exterior boundary (see the cases above)
	  bool boundary () const
	  {
		return i.boundary();
	  }

	  /**
		 \brief Identifier for boundary segment from macro grid.

		 One can attach a boundary Id to a boundary segment on the macro
		 grid. This Id will also be used for all fragments of these
		 boundary segments.

		 The numbering is defined as:
		 - Id==0 for all intersections without boundary()==false
		 - Id>=0 for all intersections without boundary()==true
     
		 The way the Identifiers are attached to the grid may differ
		 between the different grid implementations.

	  */
	  int boundaryId () const
	  {
		return i.boundaryId();
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
	  const LocalGeometry& intersectionSelfLocal () const
	  {
		return i.intersectionSelfLocal();
	  }
	  /*! @brief geometrical information about this intersection in local
		coordinates of the outside() entity.

		This method returns a Geometry object that provides a mapping from
		local coordinates of the intersection to local coordinates of the
		outside() entity.
	  */
	  const LocalGeometry& intersectionNeighborLocal () const
	  {
		return i.intersectionNeighborLocal();
	  }

	  /*! @brief geometrical information about this intersection in global coordinates.

		This method returns a Geometry object that provides a mapping from
		local coordinates of the intersection to global (world) coordinates.
	  */
	  const Geometry& intersectionGlobal () const
	  {
		return i.intersectionGlobal();
	  }

	  //! Local number of codim 1 entity in the inside() Entity where intersection is contained in
	  int numberInSelf () const
	  {
		return i.numberInSelf ();
	  }

	  //! Local number of codim 1 entity in outside() Entity where intersection is contained in
	  int numberInNeighbor () const
	  {
		return i.numberInNeighbor ();
	  }

	  /*! @brief Return an outer normal (length not necessarily 1)

		The returned vector may depend on local position within the intersection.
	  */
	  Dune::FieldVector<ctype, dimensionworld> outerNormal (const Dune::FieldVector<ctype, dimension-1>& local) const
	  {
		return i.outerNormal(local);
	  }

	  /*! @brief return outer normal scaled with the integration element
		@copydoc outerNormal
		The normal is scaled with the integration element of the intersection. This
		method is redundant but it may be more efficent to use this function
		rather than computing the integration element via intersectionGlobal().
	  */
	  Dune::FieldVector<ctype, dimensionworld> integrationOuterNormal (const Dune::FieldVector<ctype, dimension-1>& local) const
	  {
		return i.integrationOuterNormal(local);
	  }

	  /*! @brief Return unit outer normal (length == 1)

		The returned vector may depend on the local position within the intersection.
		It is scaled to have unit length.
	  */
	  Dune::FieldVector<ctype, dimensionworld> unitOuterNormal (const Dune::FieldVector<ctype, dimension-1>& local) const
	  {
		return i.unitOuterNormal(local);
	  }

	  /*! @brief return EntityPointer to the Entity on the inside of this
		intersection. That is the Entity where we started this .
	  */
	  EntityPointer inside() const
	  {
		return i.inside();
	  }

	  /*! @brief return EntityPointer to the Entity on the outside of this
		intersection. That is the neighboring Entity.

		@warning Don't call this method if there is no neighboring Entity
		(neighbor() returns false). In this case the result is undefined.
	  */
	  EntityPointer outside() const
	  {
		return i.outside();
	  }

	private:
	  const I& i;
	};

  }
}

#endif
