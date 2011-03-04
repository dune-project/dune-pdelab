// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_GRIDOPERATORUTILITIES_HH
#define DUNE_PDELAB_GRIDOPERATORUTILITIES_HH

#include <dune/pdelab/backend/backendselector.hh>

namespace Dune{
  namespace PDELab{

    //! Traits class for the grid operator.
    /**
     * This class collects types and auxilliary information about the
     * grid operator.
     *
     * \tparam GO   The grid operator.
     * \tparam GFSU The trial function space.
     * \tparam GFSV The test function space.
     * \tparam MB   The matrix backend.
     * \tparam DF   The domain (solution) field type.
     * \tparam RF   The range (residual) field type.
     * \tparam JF   The jacobian field type.
     * \tparam A    The global assembler.
     * \tparam LA   The local assembler.
     */
    template<typename GO,
             typename GFSU, typename GFSV,
             typename MB,
             typename DF, typename RF, typename JF,
             typename A, typename LA>
    struct GridOperatorTraits
    {

      //! The type of the grid operator.
      typedef GO GridOperator;

      //! The trial grid function space.
      typedef GFSU TrialGridFunctionSpace;

      //! The test grid function space.
      typedef GFSV TestGridFunctionSpace;


      //! The matrix backend of the grid operator.
      typedef MB MatrixBackend;


      //! The field type of the domain (solution).
      typedef DF DomainField;

      //! The type of the domain (solution).
      typedef typename Dune::PDELab::BackendVectorSelector<GFSU,DF>::Type Domain;


      //! The field type of the range (residual).
      typedef RF RangeField;

      //! The type of the range (residual).
      typedef typename Dune::PDELab::BackendVectorSelector<GFSV,RF>::Type Range;


      //! The field type of the jacobian.
      typedef JF JacobianField;

      //! The type of the jacobian.
      typedef typename MatrixBackend::template Matrix<GridOperator,JacobianField> Jacobian;


      //! The global assembler of the grid operator.
      typedef A Assembler;

      //! The local assembler of the grid operator.
      typedef LA LocalAssembler;

    };

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDOPERATORUTILITIES_HH
