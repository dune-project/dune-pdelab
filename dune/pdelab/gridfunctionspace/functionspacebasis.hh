// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_FUNCTIONSPACEBASIS_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_FUNCTIONSPACEBASIS_HH

/** \file helper code to unify the handling of GridFunctionSpace and Dune::Functions::GlobalBasis

    \internal this header is only intended for internal use in PDELab. Interface may change without prior notice.
 */

#include <dune/pdelab/common/concepts.hh>

namespace Dune {
namespace PDELab {

    //! \brief thin wrapper around a dune-functions basis to add information about the backend
    template<typename FSB, typename B>
    struct FunctionsBasisInfo
    {
        FunctionsBasisInfo(FSB & b) :
            _basis(b) {}

        using Basis = FSB;
        using GridView = typename FSB::GridView;
        using LocalView = typename Basis::LocalView;
        using ContainerIndex = typename LocalView::MultiIndex;

        // try to get around without Traits
        struct Traits {
            using Backend = B;
        };

        FSB & basis() { return _basis; }
        const FSB & basis() const { return _basis; }

        FSB & _basis;
    };

    /** \brief Traits class to handle GridFunctionSpace and Dune::Functions::GlobalBasis consistently

        This traits class should allow us to abstract away differences between a
        GridFunctionSpace and Dune::Functions::GlobalBasis and slowly transition
        to dune-functions basis instead of GridFunctionSpace.
    */
    template<typename GFS, typename Enable = void>
    struct BasisTraits;

#ifndef DOXYGEN
    // Specialization for GridFunctionSpace
    template<typename GFS>
    struct BasisTraits<GFS,
                     std::enable_if_t<isGridFunctionSpace<GFS>()>>
    {
        // Multiindex with which we can access entries in a DOF container
        using ContainerIndex = typename GFS::Ordering::Traits::ContainerIndex;
    };

    // Specialization for Dune::Functions::GlobalBasis
    template<typename FSB>
    struct BasisTraits<FSB,
                     std::enable_if_t<isBasisInfo<FSB>()>>
    {
        // Multiindex with which we can access entries in a DOF container
        using Basis = typename FSB::Basis;
        using LocalView = typename Basis::LocalView;
        using ContainerIndex = typename LocalView::MultiIndex;
    };
#endif

} // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_FUNCTIONSPACEBASIS_HH
