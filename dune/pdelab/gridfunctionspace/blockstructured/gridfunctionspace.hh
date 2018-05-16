// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_BLOCKSTRUCTURED_GRIDFUNCTIONSPACE_GRIDFUNCTIONSPACE_HH
#define DUNE_PDELAB_BLOCKSTRUCTURED_GRIDFUNCTIONSPACE_GRIDFUNCTIONSPACE_HH

#include <cstddef>
#include <map>
#include <ostream>
#include <set>
#include <vector>

#include <dune/common/deprecated.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/common/stdstreams.hh>
#include <dune/common/typetraits.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/interfaceswitch.hh>
#include <dune/localfunctions/common/localkey.hh>

#include <dune/typetree/typetree.hh>

#include <dune/pdelab/common/partitionviewentityset.hh>
#include <dune/pdelab/backend/interface.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/constraints/noconstraints.hh>
#include <dune/pdelab/backend/istl/descriptors.hh>
#include <dune/pdelab/gridoperator/blockstructured/tags.hh>

namespace Dune {
  namespace Blockstructured {
    /** \brief A grid function space.
     *
     *  \tparam GV   Type implementing GridView
     *  \tparam FEM  Type implementing FiniteElementMapInterface
     *  \tparam CE   Type for constraints assembler
     *  \tparam B    Backend type
     *  \tparam O    Ordering tag
     */
    template<typename GV, typename FEM, typename CE=Dune::PDELab::NoConstraints,
             typename B=Dune::PDELab::ISTL::VectorBackend<>, typename O=Dune::PDELab::DefaultLeafOrderingTag>
    class GridFunctionSpace
      : public Dune::PDELab::GridFunctionSpace<GV,FEM,CE,B,O> {
      using Base = Dune::PDELab::GridFunctionSpace<GV, FEM, CE, B, O>;

      typedef TypeTree::TransformTree<GridFunctionSpace, Dune::PDELab::gfs_to_ordering<GridFunctionSpace> > ordering_transformation;

      template<typename, typename>
      friend
      class GridFunctionSpaceBase;

    public:
      //! export Traits class
      using Traits = typename Base::Traits;

      typedef typename GV::Traits::template Codim<0>::Entity Element;
      typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;

      DUNE_DEPRECATED
      typedef O SizeTag;

      typedef O OrderingTag;

      using ImplementationTag = BlockstructuredLeafGridFunctionSpaceTag;

      typedef typename ordering_transformation::Type Ordering;

      using Base::Base;
    };
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_GRIDFUNCTIONSPACE_HH
