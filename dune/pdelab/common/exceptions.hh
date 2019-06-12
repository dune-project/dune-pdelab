// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_COMMON_EXCEPTIONS_HH
#define DUNE_PDELAB_COMMON_EXCEPTIONS_HH

#include <dune/common/exceptions.hh>

/**
 * \file
 * \brief PDELab-specific exceptions.
 */

namespace Dune {
  namespace PDELab {

    //! Base class for all PDELab exceptions.
    class Exception
      : public Dune::Exception
    {};

    //! Error when an invalid argument was passed to a function.
    class InvalidArgument
      : public Exception
    {};

    //! Assembly-related error.
    class AssemblyError
      : public Exception
    {};

    class OneStepError
      : public Exception
    {};

    class UnknownOneStepMethod
      : public OneStepError
    {};

    class LinearAlgebraError
      : public Exception
    {};

    //! GridFunctionSpace-related error.
    class GridFunctionSpaceError
      : public Exception
    {};

    //! Called a GridFunctionSpace method that requires initialization of the space.
    class UninitializedGridFunctionSpaceError
      : public GridFunctionSpaceError
    {};

    //! Called a method on a GridFunctionSpace that is not valid
    //! at its current place in the function space tree.
    class GridFunctionSpaceHierarchyError
      : public GridFunctionSpaceError
    {};

    //! Ordering-related error.
    class OrderingError
      : public Exception
    {};

    //! Error related to the logical structure of an Ordering.
    class OrderingStructureError
      : public OrderingError
    {};

    //! A PermutedOrdering got a permutation vector of the wrong size.
    class PermutedOrderingSizeError
      : public OrderingError
    {};

    //! The block size of a ChunkedBlockOrdering does not divide the block count of the underlying ordering.
    class ChunkedBlockOrderingSizeError
      : public OrderingError
    {};

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_COMMON_EXCEPTIONS_HH
