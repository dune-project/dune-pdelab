// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_CONSTRAINTS_COMMON_CONSTRAINTSPARAMETERS_HH
#define DUNE_PDELAB_CONSTRAINTS_COMMON_CONSTRAINTSPARAMETERS_HH

#include <dune/common/fvector.hh>
#include <dune/typetree/typetree.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup GridFunctionSpace
    //! \ingroup PDELab
    //! \{

    /**
     * Interface for the constraints parameters describing Dirichlet constraints.
     *
     * \note This class not only describes the required interface of the parameter class,
     *       but can also be used as a convenient standard implementation that will add
     *       Dirichlet constraints to all locations it is queried about.
     */
    struct DirichletConstraintsParameters :
      public TypeTree::LeafNode
    {

      /**
       * Indicates whether the given position should be Dirichlet-constrained.
       *
       * \param intersection The grid intersection containing the queried location.
       * \param coord        The position of the queried location in local coordinates
       *                     of the intersection.
       *
       * \returns            true iff the given location should have a Dirichlet constraint.
       */
      template<typename I>
      bool isDirichlet(const I & intersection, const FieldVector<typename I::ctype, I::dimension-1> & coord) const
      {
        return true;
      }

      /**
       * Indicates whether the given position should be Neumann-constrained.
       *
       * Most of the time, this method will be equivalent to !isDirichlet(...), but
       * sometimes (in particular in multi-domain scenarios), both methods may return false.
       *
       * \param intersection The grid intersection containing the queried location.
       * \param coord        The position of the queried location in local coordinates
       *                     of the intersection.
       *
       * \returns            true iff the given location should have a Neumann constraint.
       */
      template<typename I>
      bool isNeumann(const I & intersection, const FieldVector<typename I::ctype, I::dimension-1> & coord) const
      {
        return false;
      }


      /**
       * Sets the current time (only required for instationary problems).
       *
       * \note This method only needs to be implement for instationary problems.
       */
      template<typename T>
      void setTime(const T& time)
      {}

    };

    /**
     * Dirichlet constraints parameters convenience implementation that does not apply
     * Dirichlet constraints anywhere.
     *
     * \sa DirichletConstraintsParameters
     */
    struct NoDirichletConstraintsParameters :
      public TypeTree::LeafNode
    {

      /**
       * Predicate implementation that will always return false.
       */
      template<typename I>
      bool isDirichlet(const I & intersection, const FieldVector<typename I::ctype, I::dimension-1> & coord) const
      {
        return false;
      }

      /**
       * Predicate implementation that will always return true.
       */
      template<typename I>
      bool isNeumann(const I & intersection, const FieldVector<typename I::ctype, I::dimension-1> & coord) const
      {
        return true;
      }

      /**
       * Sets the current time (only required for instationary problems).
       *
       * \note This method only needs to be implement for instationary problems.
       */
      template<typename T>
      void setTime(const T& time)
      {}

    };


    /**
     * Interface for the constraints parameters describing flux (Neumann) constraints.
     *
     * \note This class not only describes the required interface of the parameter class,
     *       but can also be used as a convenient standard implementation that will add
     *       Neumann constraints to all locations it is queried about.
     */
    struct FluxConstraintsParameters :
      public TypeTree::LeafNode
    {

      /**
       * Indicates whether the given position should be Neumann-constrained.
       *
       * \param intersection The grid intersection containing the queried location.
       * \param coord        The position of the queried location in local coordinates
       *                     of the intersection.
       *
       * \returns            true iff the given location should have a Neumann constraint.
       */
      template<typename I>
      bool isNeumann(const I & intersection, const FieldVector<typename I::ctype, I::dimension-1> & coord) const
      {
        return true;
      }

      /**
       * Sets the current time (only required for instationary problems).
       *
       * \note This method only needs to be implement for instationary problems.
       */
      template<typename T>
      void setTime(const T& time)
      {}

    };

    /**
     * Flux (Neumann) constraints parameters convenience implementation that does not apply
     * Neumann constraints anywhere.
     *
     * \sa FluxConstraintsParameters
     */
    struct NoFluxConstraintsParameters :
      public TypeTree::LeafNode
    {

      /**
       * Predicate implementation that will always return false.
       */
      template<typename I>
      bool isNeumann(const I & intersection, const FieldVector<typename I::ctype, I::dimension-1> & coord) const
      {
        return false;
      }

      /**
       * Sets the current time (only required for instationary problems).
       *
       * \note This method only needs to be implement for instationary problems.
       */
      template<typename T>
      void setTime(const T& time)
      {}

    };

    /**
     * Adapter class that deduces Neumann constraints from a model of DirichletConstraintsParameters
     * by applying the standard logic that every location on the boundary is either Neumann or
     * Dirichlet constrained.
     *
     * \warning This adapter will only work for the standard case of \f$\partial \Omega = \Gamma_D \cup \Gamma_N\f$.
     *          Do not try to use it for more advanced problems such as flow with flux, outflow and Dirichlet constraints!
     *
     * \tparam DirichletConstraintsParameters The type of the underlying Dirichlet constraints parameters implementation.
     */
    template<typename DirichletConstraintsParameters>
    struct FluxFromDirichletConstraintsAdapter :
      public TypeTree::LeafNode
    {

      /**
       * Forwards to the underlying DirichletConstraintsParameters implementation.
       */
      template<typename I>
      bool isDirichlet(const I & intersection, const FieldVector<typename I::ctype, I::dimension-1> & coord) const
      {
        return _dirichletConstraintsParameters.isDirichlet(intersection,coord);
      }

      /**
       * Returns the opposite of isDirichlet().
       */
      template<typename I>
      bool isNeumann(const I & intersection, const FieldVector<typename I::ctype, I::dimension-1> & coord) const
      {
        return !_dirichletConstraintsParameters.isDirichlet(intersection,coord);
      }

      /**
       * Forwards the new time to the underlying DirichletConstraintsParameters object.
       */
      template<typename T>
      void setTime(const T& time)
      {
        _dirichletConstraintsParameters.setTime(time);
      }

      FluxFromDirichletConstraintsAdapter(DirichletConstraintsParameters& dirichletConstraintsParameters)
        : _dirichletConstraintsParameters(dirichletConstraintsParameters)
      {}

    private:

      DirichletConstraintsParameters& _dirichletConstraintsParameters;

    };

    /**
     * Adapter class that deduces Dirichlet constraints from a model of FluxConstraintsParameters
     * by applying the standard logic that every location on the boundary is either Neumann or
     * Dirichlet constrained.
     *
     * \warning This adapter will only work for the standard case of \f$\partial \Omega = \Gamma_D \cup \Gamma_N\f$.
     *          Do not try to use it for more advanced problems such as flow with flux, outflow and Dirichlet constraints!
     *
     * \tparam FluxConstraintsParameters The type of the underlying Neumann constraints parameters implementation.
     */
    template<typename FluxConstraintsParameters>
    struct DirichletFromFluxConstraintsAdapter :
      public TypeTree::LeafNode
    {

      /**
       * Returns the opposite of isNeumann().
       */
      template<typename I>
      bool isDirichlet(const I & intersection, const FieldVector<typename I::ctype, I::dimension-1> & coord) const
      {
        return !_fluxConstraintsParameters.isNeumann(intersection,coord);
      }

      /**
       * Forwards to the underlying FluxConstraintsParameters implementation
       */
      template<typename I>
      bool isNeumann(const I & intersection, const FieldVector<typename I::ctype, I::dimension-1> & coord) const
      {
        return _fluxConstraintsParameters.isNeumann(intersection,coord);
      }

      /**
       * Forwards the new time to the underlying FluxConstraintsParameters object.
       */
      template<typename T>
      void setTime(const T& time)
      {
        _fluxConstraintsParameters.setTime(time);
      }

      DirichletFromFluxConstraintsAdapter(FluxConstraintsParameters& fluxConstraintsParameters)
        : _fluxConstraintsParameters(fluxConstraintsParameters)
      {}

    private:

      FluxConstraintsParameters& _fluxConstraintsParameters;

    };

    /**
     * Empty parameters to allow running the assembler for NoConstraints
     *
     * \note This class does not do anything. It only works, because the
     *       NoConstraints don't even try to assemble anything.
     *
     * \remark it is also used internally for the
     *       Dune::PDELab::constraints assembler TMP without any parameters
     */
    class NoConstraintsParameters : public TypeTree::LeafNode {};

    //! \}

  }
}

#endif // DUNE_PDELAB_CONSTRAINTS_COMMON_CONSTRAINTSPARAMETERS_HH
