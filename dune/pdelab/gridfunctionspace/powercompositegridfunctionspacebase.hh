// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_POWERCOMPOSITEGRIDFUNCTIONSPACEBASE_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_POWERCOMPOSITEGRIDFUNCTIONSPACEBASE_HH

#include <algorithm>
#include <cstddef>
#include <functional>
#include <vector>

#include <dune/typetree/fixedcapacitystack.hh>
#include <dune/typetree/traversal.hh>
#include <dune/typetree/traversalutilities.hh>
#include <dune/typetree/utility.hh>
#include <dune/typetree/transformation.hh>
#include <dune/typetree/visitor.hh>

#include <dune/pdelab/constraints/common/constraintstransformation.hh>
#include <dune/pdelab/constraints/noconstraints.hh>
#include <dune/pdelab/gridfunctionspace/tags.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspacebase.hh>
#include <dune/pdelab/ordering/lexicographicordering.hh>
#include <dune/pdelab/ordering/entityblockedlocalordering.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup GridFunctionSpace grid function space
    //! \ingroup PDELab
    //! \{

    //! Trait class for the multi component grid function spaces
    template<typename G, typename B, typename O, std::size_t k>
    struct PowerCompositeGridFunctionSpaceTraits
    {
      enum{
        //! \brief True if this grid function space is composed of others.
        isComposite = 1,
        //! \brief number of child spaces
        noChilds = k
      };

      const static std::size_t CHILDREN = k;

      using EntitySet = G;

      using GridView = typename EntitySet::GridView;

      //! \brief the grid view where grid function is defined upon
      using GridViewType = GridView;

      //! \brief vector backend
      typedef B BackendType;

      typedef B Backend;

      //! \brief mapper
      typedef O MapperType;

      typedef O OrderingTag;

      //! \brief short cut for size type exported by Backend
      typedef typename B::size_type SizeType;
    };

    //! Mixin class providing common functionality of PowerGridFunctionSpace and CompositeGridFunctionSpace
    template<typename GridFunctionSpace, typename GV, typename B, typename O, std::size_t k>
    class PowerCompositeGridFunctionSpaceBase
      : public GridFunctionSpaceBase<
                 GridFunctionSpace,
                 PowerCompositeGridFunctionSpaceTraits<GV,B,O,k>
                 >
    {

#ifndef DOXYGEN

      const GridFunctionSpace& gfs() const
      {
        return static_cast<const GridFunctionSpace&>(*this);
      }

      GridFunctionSpace& gfs()
      {
        return static_cast<GridFunctionSpace&>(*this);
      }

#endif // DOXYGEN

    public:

      //! export traits class
      typedef PowerCompositeGridFunctionSpaceTraits<GV,B,O,k> Traits;

    private:

      typedef GridFunctionSpaceBase<GridFunctionSpace,Traits> BaseT;

    public:

      typedef O OrderingTag;

      // TODO: Do not just use constraints from child 0!
      //! extract type for storing constraints
      template<typename E>
      struct ConstraintsContainer
      {
        typedef typename std::conditional<
          std::is_same<
            typename GridFunctionSpace::template Child<0>::type::template ConstraintsContainer<E>::Type,
            EmptyTransformation
            >::value,
          EmptyTransformation,
          ConstraintsTransformation<
            typename GridFunctionSpace::Ordering::Traits::DOFIndex,
            typename GridFunctionSpace::Ordering::Traits::ContainerIndex,
            E
            >
          >::type Type;
      };

      //! get grid view
      const typename Traits::GridView& gridView () const
      {
        return gfs().template child<0>().gridView();
      }

      //! get grid view partition
      const typename Traits::EntitySet& entitySet () const
      {
        return gfs().template child<0>().entitySet();
      }

      PowerCompositeGridFunctionSpaceBase(const B& backend, const OrderingTag& ordering_tag)
        : BaseT(backend,ordering_tag)
      {}

    };

  }

}
//! \}
#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_POWERCOMPOSITEGRIDFUNCTIONSPACEBASE_HH
