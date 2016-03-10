// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_GRIDOPERATOR_COMMON_GRIDOPERATORUTILITIES_HH
#define DUNE_PDELAB_GRIDOPERATOR_COMMON_GRIDOPERATORUTILITIES_HH

#include <dune/pdelab/backend/interface.hh>

namespace Dune{
  namespace PDELab{

    //! Traits class for the grid operator.
    /**
     * This class collects types and auxilliary information about the
     * grid operator.
     *
     * \tparam GFSU The trial function space.
     * \tparam GFSV The test function space.
     * \tparam MB   The matrix backend.
     * \tparam DF   The domain (solution) field type.
     * \tparam RF   The range (residual) field type.
     * \tparam JF   The jacobian field type.
     * \tparam CU   The type of the trial grid function space constraints.
     * \tparam CV   The type of the test grid function space constraints.
     * \tparam A    The global assembler.
     * \tparam LA   The local assembler.
     */
    template<typename GFSU, typename GFSV,
             typename MB,
             typename DF, typename RF, typename JF,
             typename CU, typename CV,
             typename A, typename LA>
    struct GridOperatorTraits
    {

      //! The trial grid function space.
      typedef GFSU TrialGridFunctionSpace;

      //! The test grid function space.
      typedef GFSV TestGridFunctionSpace;


      //! The type of the trial grid function space constraints.
      typedef CU TrialGridFunctionSpaceConstraints;

      //! The type of the test grid function space constraints.
      typedef CV TestGridFunctionSpaceConstraints;


      //! The matrix backend of the grid operator.
      typedef MB MatrixBackend;


      //! The field type of the domain (solution).
      typedef DF DomainField;

      //! The type of the domain (solution).
      using Domain = Dune::PDELab::Backend::Vector<GFSU,DF>;


      //! The field type of the range (residual).
      typedef RF RangeField;

      //! The type of the range (residual).
      using Range = Dune::PDELab::Backend::Vector<GFSV,RF>;


      //! The field type of the jacobian.
      typedef JF JacobianField;

      //! The type of the jacobian.
      using Jacobian = Dune::PDELab::Backend::Matrix<MB,Domain,Range,JF>;


      //! The global assembler of the grid operator.
      typedef A Assembler;

      //! The local assembler of the grid operator.
      typedef LA LocalAssembler;

    };

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDOPERATOR_COMMON_GRIDOPERATORUTILITIES_HH
