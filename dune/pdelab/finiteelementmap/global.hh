// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_FINITEELEMENTMAP_GLOBAL_HH
#define DUNE_PDELAB_FINITEELEMENTMAP_GLOBAL_HH

#include <dune/pdelab/common/countingptr.hh>
#include <dune/pdelab/finiteelementmap/finiteelementmap.hh>

namespace Dune {
  namespace PDELab {

    //! \brief Generic finite element map for global finite elements created
    //!        with a geometry
    /**
     * \tparam Factory Finite element factory class that supports make(const
     *                 Geometry&).
     */
    template<class Factory>
    class GeometryFiniteElementMap :
      public Countable
    {
      Factory& factory;

    public:
      typedef FiniteElementMapTraits<typename Factory::FiniteElement> Traits;

      //! construct GeometryFiniteElementMap
      /**
       * \param factory_ Reference to a factory object.
       *
       * \note The constructed finite element map object stores a reference to
       *       to factory object.  The finite element map object becomes
       *       invalid as soon the reference to the factory becomes invalid.
       *       The only allowed operation on an invalid finite element map is
       *       calling the destructor, all other operations result in
       *       undefined behaviour.
       */
      GeometryFiniteElementMap(Factory& factory_) : factory(factory_) {}

      //! Return finite element for the given entity.
      /**
       * \param e Grid element to create a finite element for.
       *
       * The returned value is valid for as long as both this finite element
       * map as well as the reference to the grid element are valid.
       */
      template<class Element>
      typename Traits::FiniteElementType find(const Element& e) const {
        return factory.make(e.geometry());
      }
    };

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_FINITEELEMENTMAP_GLOBAL_HH
