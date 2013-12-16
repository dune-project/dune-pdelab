// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_HH

#include <cstddef>
#include <map>
#include <ostream>
#include <set>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/common/static_assert.hh>
#include <dune/common/stdstreams.hh>
#include <dune/common/typetraits.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/interfaceswitch.hh>
#include <dune/localfunctions/common/localkey.hh>

#include <dune/typetree/typetree.hh>

#include <dune/pdelab/backend/backendselector.hh>
#include <dune/pdelab/backend/istl/descriptors.hh>
#include <dune/pdelab/common/geometrywrapper.hh>
#include <dune/pdelab/gridfunctionspace/compositegridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/datahandleprovider.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspacebase.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/pdelab/gridfunctionspace/powergridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/tags.hh>
#include <dune/pdelab/gridfunctionspace/utility.hh>
#include <dune/pdelab/ordering/gridviewordering.hh>
#include <dune/pdelab/ordering/lexicographicordering.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup GridFunctionSpace grid function space
    //! \ingroup PDELab
    //! \{

    //=======================================
    // grid function space : single component case
    //=======================================

    //! collect types exported by a leaf grid function space
    /**
     * This is based on a global FiniteElementMap
     */
    template<typename G, typename L, typename C, typename B, typename O>
    struct GridFunctionSpaceTraits
    {
      //! True if this grid function space is composed of others.
      static const bool isComposite = false;

      //! the grid view where grid function is defined upon
      typedef G GridViewType;

      typedef G GridView;

      //! vector backend
      typedef B BackendType;

      typedef B Backend;

      //! short cut for size type exported by Backend
      typedef typename B::size_type SizeType;

      //! finite element map
      typedef L FiniteElementMapType;

      //! finite element map
      typedef L FiniteElementMap;

      //! finite element
      typedef typename L::Traits::FiniteElementType FiniteElementType;

      typedef typename L::Traits::FiniteElementType FiniteElement;

      //! type representing constraints
      typedef C ConstraintsType;

      //! tag describing the ordering.
      /**
       * The tag type may contain additional constants and typedefs to
       * control the behavior of the created ordering.
       */
      typedef O OrderingTag;

    };

    /** \brief A grid function space.
     *
     *  \tparam GV   Type implementing GridView
     *  \tparam FEM  Type implementing FiniteElementMapInterface
     *  \tparam CE   Type for constraints assembler
     *  \tparam B    Backend type
     *  \tparam P    Parameter type. Possible types are
     * \link GridFunctionGeneralMapper \endlink (arbitrary number of unknowns per
     * entity) or \link GridFunctionRestrictedMapper \endlink (fixed number of unknowns per
     * entity) or \link GridFunctionStaticSize \endlink (number of unknowns per
     * entity, known at compile-time)
     */
    template<typename GV, typename FEM, typename CE=NoConstraints,
             typename B=ISTLVectorBackend<>, typename P=DefaultLeafOrderingTag>
    class GridFunctionSpace
      : public TypeTree::LeafNode
      , public GridFunctionSpaceBase<
                 GridFunctionSpace<GV,FEM,CE,B,P>,
                 GridFunctionSpaceTraits<GV,FEM,CE,B,P>
                 >
      , public GridFunctionOutputParameters
      , public DataHandleProvider<GridFunctionSpace<GV,FEM,CE,B,P> >
    {

      typedef TypeTree::TransformTree<GridFunctionSpace,gfs_to_ordering<GridFunctionSpace> > ordering_transformation;

      template<typename,typename>
      friend class GridFunctionSpaceBase;

    public:
      //! export Traits class
      typedef GridFunctionSpaceTraits<GV,FEM,CE,B,P> Traits;

    private:

      typedef GridFunctionSpaceBase<GridFunctionSpace,Traits> BaseT;

    public:

      typedef typename GV::Traits::template Codim<0>::Entity Element;
      typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;

      typedef P SizeTag;

      typedef P OrderingTag;

      typedef LeafGridFunctionSpaceTag ImplementationTag;

      typedef typename ordering_transformation::Type Ordering;

       //! extract type for storing constraints
      template<typename E>
      struct ConstraintsContainer
      {

        //! \brief define Type as the Type of a container of E's
        typedef typename conditional<
          is_same<
            CE,
            NoConstraints
            >::value,
          EmptyTransformation,
          ConstraintsTransformation<typename Ordering::Traits::DOFIndex,typename Ordering::Traits::ContainerIndex,E>
          >::type Type;

      private:
        ConstraintsContainer () {}
      };

      //! constructor
      GridFunctionSpace (const GV& gridview, const FEM& fem, const CE& ce_, const B& backend = B(), const OrderingTag& ordering_tag = OrderingTag())
        : BaseT(backend,ordering_tag)
        , defaultce(ce_)
        , gv(gridview)
        , pfem(stackobject_to_shared_ptr(fem))
        , ce(ce_)
      {
      }

      //! constructor
      GridFunctionSpace (const GV& gridview, const shared_ptr<const FEM>& fem, const CE& ce_, const B& backend = B(), const OrderingTag& ordering_tag = OrderingTag())
        : BaseT(backend,ordering_tag)
        , defaultce(ce_)
        , gv(gridview)
        , pfem(fem)
        , ce(ce_)
      {}

      //! constructor
      GridFunctionSpace (const GV& gridview, const FEM& fem, const B& backend = B(), const OrderingTag& ordering_tag = OrderingTag())
        : BaseT(backend,ordering_tag)
        , gv(gridview)
        , pfem(stackobject_to_shared_ptr(fem))
        , ce(defaultce)
      {}

      //! constructor
      GridFunctionSpace (const GV& gridview, const shared_ptr<const FEM>& fem, const B& backend = B(), const OrderingTag& ordering_tag = OrderingTag())
        : BaseT(backend,ordering_tag)
        , gv(gridview)
        , pfem(fem)
        , ce(defaultce)
      {}

      //! get grid view
      const GV& gridView () const
      {
        return gv;
      }

      //! get finite element map
      const FEM& finiteElementMap () const
      {
        return *pfem;
      }

      //! get finite element map
      shared_ptr<const FEM> finiteElementMapStorage () const
      {
        return pfem;
      }

      // return constraints engine
      const typename Traits::ConstraintsType& constraints () const
      {
        return ce;
      }

      //------------------------------

      //! Direct access to the DOF ordering.
      const Ordering &ordering() const
      {
        if (!this->isRootSpace())
          {
            DUNE_THROW(GridFunctionSpaceHierarchyError,
                       "Ordering can only be obtained for root space in GridFunctionSpace tree.");
          }
        if (!_ordering)
          {
            create_ordering();
            this->update(*_ordering);
          }
        return *_ordering;
      }

      //! Direct access to the DOF ordering.
      Ordering &ordering()
      {
        if (!this->isRootSpace())
          {
            DUNE_THROW(GridFunctionSpaceHierarchyError,
                       "Ordering can only be obtained for root space in GridFunctionSpace tree.");
          }
        if (!_ordering)
          {
            create_ordering();
            this->update(*_ordering);
          }
        return *_ordering;
      }

      //! Direct access to the storage of the DOF ordering.
      shared_ptr<const Ordering> orderingStorage() const
      {
        if (!this->isRootSpace())
          {
            DUNE_THROW(GridFunctionSpaceHierarchyError,
                       "Ordering can only be obtained for root space in GridFunctionSpace tree.");
          }
        if (!_ordering)
          {
            create_ordering();
            this->update(*_ordering);
          }
        return _ordering;
      }

      //! Direct access to the storage of the DOF ordering.
      shared_ptr<Ordering> orderingStorage()
      {
        if (!this->isRootSpace())
          {
            DUNE_THROW(GridFunctionSpaceHierarchyError,
                       "Ordering can only be obtained for root space in GridFunctionSpace tree.");
          }
        if (!_ordering)
          {
            create_ordering();
            this->update(*_ordering);
          }
        return _ordering;
      }

    private:

      // This method here is to avoid a double update of the Ordering when the user calls
      // GFS::update() before GFS::ordering().
      void create_ordering() const
      {
        _ordering = make_shared<Ordering>(ordering_transformation::transform(*this));
      }

      CE defaultce;
      const GV gv;
      shared_ptr<FEM const> pfem;
      const CE& ce;

      mutable shared_ptr<Ordering> _ordering;
    };


  } // namespace PDELab
} // namespace Dune

#endif
