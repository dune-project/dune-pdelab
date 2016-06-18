// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_LOCALOPERATOR_VARIABLEFACTORIES_HH
#define DUNE_PDELAB_LOCALOPERATOR_VARIABLEFACTORIES_HH

#include <vector>

namespace Dune {
  namespace PDELab {

    //! return a container for basis evaluations
    template<typename LFS>
    std::vector<typename LFS::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType>
    makeValueContainer (const LFS& lfs)
    {
      return std::vector<typename LFS::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType>(lfs.size());
    }

    //! return a container for Jacobian evaluations
    template<typename LFS>
    std::vector<typename LFS::Traits::FiniteElementType::Traits::LocalBasisType::Traits::JacobianType>
    makeJacobianContainer (const LFS& lfs)
    {
      return std::vector<typename LFS::Traits::FiniteElementType::Traits::LocalBasisType::Traits::JacobianType>(lfs.size());
    }

    //! return a zero value of RangeFieldType of the basis
    template<typename LFS>
    typename LFS::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType
    makeZeroBasisFieldValue (const LFS& lfs)
    {
      typename LFS::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType x(0.0);
      return x;
    }

    //! return a zero value of RangeType of the basis
    template<typename LFS>
    typename LFS::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType
    makeZeroBasisValue (const LFS& lfs)
    {
      typename LFS::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType x(0.0);
      return x;
    }

    //! return a zero value of JacobianType of the basis
    template<typename LFS>
    typename LFS::Traits::FiniteElementType::Traits::LocalBasisType::Traits::JacobianType
    makeZeroJacobianValue (const LFS& lfs)
    {
      typename LFS::Traits::FiniteElementType::Traits::LocalBasisType::Traits::JacobianType x(0.0);
      return x;
    }
  }
}

#endif // DUNE_PDELAB_LOCALOPERATOR_VARIABLEFACTORIES_HH
