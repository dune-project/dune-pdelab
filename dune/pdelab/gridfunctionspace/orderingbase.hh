// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_ORDERINGBASE_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_ORDERINGBASE_HH

#include <cstddef>

#include <dune/common/documentation.hh>
#include <dune/common/static_assert.hh>

#include <dune/pdelab/common/typetree/compositenodemacros.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup GridFunctionSpace
    //! \ingroup PDELab
    //! \{

    //! Transform PowerGridFunctionSpaces to ordering
    template<class OrderingTag>
    class TransformPowerGFSToOrdering {
      dune_static_assert(AlwaysFalse<OrderingTag>::value, "Can't find an "
                         "ordering for this particular ordering tag.");

    public:
      template<class GFSTraits, class TransformedChild, std::size_t k>
      struct result {
        typedef ImplementationDefined type;
      };
    };

    //! Transform CompositeGridFunctionSpaces to ordering
    template<class OrderingTag>
    class TransformCompositeGFSToOrdering {
      dune_static_assert(AlwaysFalse<OrderingTag>::value, "Can't find an "
                         "ordering for this particular ordering tag.");

    public:
      template<class GFSTraits, DUNE_TYPETREE_COMPOSITENODE_TEMPLATE_CHILDREN>
      struct result {
        typedef ImplementationDefined type;
      };
    };

   //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_ORDERINGBASE_HH
