// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_FINITEELEMENTMAP_POWERFEM_HH
#define DUNE_PDELAB_FINITEELEMENTMAP_POWERFEM_HH

#include <cstddef>

#include <dune/localfunctions/meta/power.hh>

#include <dune/pdelab/finiteelementmap/finiteelementmap.hh>

namespace Dune {
  namespace PDELab {

    //! FiniteElementMap for PowerFiniteElements
    /**
     * \tparam BackendFEM Map for finite elements that should be raised to a
     *                    power.
     * \tparam dimR       Power to raise the backend FiniteElements to.
     */
    template<class BackendFEM, std::size_t dimR>
    class PowerFiniteElementMap
    {
      typedef typename BackendFEM::Traits::FiniteElementType BackendFE;
      typedef PowerFiniteElementFactory<BackendFE, dimR> Factory;

      const BackendFEM &backend;
      static const Factory factory;

    public:
      //! export Traits
      typedef FiniteElementMapTraits<typename Factory::FiniteElement> Traits;

      //! construct PowerFiniteElementMap
      /**
       * \param backend_ Reference to a finite element map for the underlying
       *                 finite elements.
       *
       * \note The constructed finite element map object stores a reference to
       *       to backend object.  The finite element map object becomes
       *       invalid as soon the reference to the backend becomes invalid.
       *       The only allowed operation on an invalid finite element map is
       *       calling the destructor, all other operations result in
       *       undefined behaviour.
       */
      PowerFiniteElementMap(const BackendFEM& backend_) : backend(backend_) { }

      //! Return finite element for the given entity.
      /**
       * \param e Grid element to create a finite element for.
       *
       * The returned value is valid for as long as both this finite element
       * map as well as the reference to the grid element are valid.
       */
      template<class Element>
      typename Traits::FiniteElementType find(const Element& e) const
      { return factory.make(backend.find(e)); }

    };

    template<class BackendFEM, std::size_t dimR>
    const typename PowerFiniteElementMap<BackendFEM, dimR>::Factory
    PowerFiniteElementMap<BackendFEM, dimR>::factory = Factory();

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_FINITEELEMENTMAP_POWERFEM_HH
