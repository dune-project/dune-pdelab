//
// Created by marckoch on 14/05/18.
//

#ifndef DUNE_PDELAB_BLOCKSTRUCTURED_TAGS_HH
#define DUNE_PDELAB_BLOCKSTRUCTURED_TAGS_HH

#include <dune/pdelab/gridfunctionspace/tags.hh>
#include <dune/pdelab/ordering/transformations.hh>

namespace Dune {
  namespace Blockstructured {

    struct BlockstructuredLeafGridFunctionSpaceTag : public Dune::PDELab::LeafGridFunctionSpaceTag {};

    template<typename GridFunctionSpace, typename Params>
    Dune::PDELab::leaf_gfs_to_ordering_descriptor<
        GridFunctionSpace,
        Dune::PDELab::gfs_to_ordering<Params>,
        typename GridFunctionSpace::Traits::OrderingTag
    >
    registerNodeTransformation(GridFunctionSpace*, Dune::PDELab::gfs_to_ordering<Params>*, BlockstructuredLeafGridFunctionSpaceTag*);
  }
}
#endif //DUNE_PDELAB_TAGS_HH
