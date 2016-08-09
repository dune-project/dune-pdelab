// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_DUNEFUNCTIONSGRIDFUNCTIONSPACE_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_DUNEFUNCTIONSGRIDFUNCTIONSPACE_HH

#include <cstddef>
#include <map>
#include <ostream>
#include <set>
#include <vector>

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
#include <dune/pdelab/backend/istl/descriptors.hh>
#include <dune/pdelab/constraints/noconstraints.hh>
#include <dune/pdelab/gridfunctionspace/compositegridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/datahandleprovider.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/pdelab/gridfunctionspace/powergridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/tags.hh>
#include <dune/pdelab/ordering/gridviewordering.hh>
#include <dune/pdelab/ordering/lexicographicordering.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup GridFunctionSpace grid function space
    //! \ingroup PDELab
    //! \{

    //! \brief collect types exported by a finite element map
    template<class T>
    struct DuneFunctionsFiniteElementMapTraits
    {
      //! Type of finite element from local functions
      typedef T FiniteElementType;

      //! Type of finite element from local functions
      typedef T FiniteElement;
    };

    //! wrap up element from local functions
    template<typename DFBasis>
    class DuneFunctionsFiniteElementMap
    {
    public:

      using Traits = DuneFunctionsFiniteElementMapTraits<typename DFBasis::LocalView::Tree::FiniteElement>;

      DuneFunctionsFiniteElementMap(const DFBasis& basis)
      : basis_(basis),
        elementMapper_(basis.gridView()),
        localView_(elementMapper_.size())
      {}

      bool fixedSize() const
      {
        return false;
      }

      bool hasDOFs(int codim) const
      {
        return true;
      }

     /** \brief Return local basis for the given entity.
      */
      template<class EntityType>
      const typename Traits::FiniteElementType&
      find (const EntityType& e) const
      {
        auto idx = elementMapper_.index(e);

        if (not localView_[idx])
        {
          localView_[idx] = std::make_unique<typename DFBasis::LocalView>(basis_.localView());
          localView_[idx]->bind(e);
        }

        return localView_[idx]->tree().finiteElement();
      }

      std::size_t size(GeometryType gt) const
      {
        DUNE_THROW(NotImplemented, "!");
        return 0;
      }

      std::size_t maxLocalSize() const
      {
        DUNE_THROW(NotImplemented, "!");
        return 0;
      }

    private:
      const DFBasis& basis_;

      MultipleCodimMultipleGeomTypeMapper<typename DFBasis::GridView, MCMGElementLayout> elementMapper_;

      // The dune-functions local view for each grid element
      mutable std::vector<std::unique_ptr<typename DFBasis::LocalView> > localView_;
    };


    //=======================================
    // grid function space : single component case
    //=======================================

    //! collect types exported by a leaf grid function space
    /**
     * This is based on a global FiniteElementMap
     */
    template<typename G, typename L, typename C, typename B, typename O>
    struct DuneFunctionsGridFunctionSpaceTraits
    {
      //! True if this grid function space is composed of others.
      static const bool isComposite = false;

      //! the grid view where grid function is defined upon
      using GridView = impl::GridView<G>;

      //! the entity set of this function space.
      using EntitySet = impl::EntitySet<G>;

      using GridViewType = GridView;

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

    /** \brief A pdelab grid function space implemented by a dune-functions function space basis
     *
     *  \tparam DFBasis A dune-functions function space basis
     *  \tparam CE   Type for constraints assembler
     *  \tparam B    Backend type
     *  \tparam P    Parameter type. Possible types are
     * \link GridFunctionGeneralMapper \endlink (arbitrary number of unknowns per
     * entity) or \link GridFunctionRestrictedMapper \endlink (fixed number of unknowns per
     * entity) or \link GridFunctionStaticSize \endlink (number of unknowns per
     * entity, known at compile-time)
     */
    template<typename DFBasis, typename CE=NoConstraints,
             typename B=istl::VectorBackend<>, typename P=DefaultLeafOrderingTag>
    class DuneFunctionsGridFunctionSpace
      : public TypeTree::LeafNode
      , public GridFunctionSpaceBase<
                 DuneFunctionsGridFunctionSpace<DFBasis,CE,B,P>,
                 DuneFunctionsGridFunctionSpaceTraits<typename DFBasis::GridView,
                                                      DuneFunctionsFiniteElementMap<DFBasis>,CE,B,P>
                 >
      , public GridFunctionOutputParameters
      , public DataHandleProvider<DuneFunctionsGridFunctionSpace<DFBasis,CE,B,P> >
    {
      using GV = typename DFBasis::GridView;

      typedef TypeTree::TransformTree<DuneFunctionsGridFunctionSpace,gfs_to_ordering<DuneFunctionsGridFunctionSpace> > ordering_transformation;

      template<typename,typename>
      friend class GridFunctionSpaceBase;

      using FEM = DuneFunctionsFiniteElementMap<DFBasis>;

    public:
      //! export Traits class
      typedef DuneFunctionsGridFunctionSpaceTraits<GV,FEM,CE,B,P> Traits;

    private:

      typedef GridFunctionSpaceBase<DuneFunctionsGridFunctionSpace,Traits> BaseT;

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
        typedef typename std::conditional<
          std::is_same<
            CE,
            NoConstraints
            >::value,
          EmptyTransformation,
          ConstraintsTransformation<typename Ordering::Traits::DOFIndex,typename Ordering::Traits::ContainerIndex,E>
          >::type Type;

      private:
        ConstraintsContainer () {}
      };

      // ****************************************************************************************************
      // Construct from a dune-functions basis
      // ****************************************************************************************************

      //! constructor
      DuneFunctionsGridFunctionSpace (const DFBasis& dfBasis, const CE& ce, const B& backend = B(), const OrderingTag& ordering_tag = OrderingTag())
        : BaseT(backend,ordering_tag)
        , _es(dfBasis.gridView())
        , dfBasis_(dfBasis)
        , pfem(std::make_shared<FEM>(dfBasis))
        , _pce(stackobject_to_shared_ptr(ce))
      {
      }

      //! constructor
      DuneFunctionsGridFunctionSpace (const DFBasis& dfBasis, const B& backend = B(), const OrderingTag& ordering_tag = OrderingTag())
        : BaseT(backend,ordering_tag)
        , _es(dfBasis.gridView())
        , pfem(std::make_shared<FEM>(dfBasis))
        , _pce(std::make_shared<CE>())
      {}


      //! get grid view
      const typename Traits::GridView& gridView () const
      {
        return _es.gridView();
      }

      //! get EntitySet
      const typename Traits::EntitySet& entitySet () const
      {
        return _es;
      }

      //! get finite element map
      const FEM& finiteElementMap () const
      {
        return *pfem;
      }

      //! get finite element map
      std::shared_ptr<const FEM> finiteElementMapStorage () const
      {
        return pfem;
      }

      //! return constraints engine
      const typename Traits::ConstraintsType& constraints () const
      {
        return *_pce;
      }

      //! return storage of constraints engine
      std::shared_ptr<const CE> constraintsStorage () const
      {
        return _pce;
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
      std::shared_ptr<const Ordering> orderingStorage() const
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
      std::shared_ptr<Ordering> orderingStorage()
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
        _ordering = std::make_shared<Ordering>(ordering_transformation::transform(*this));
      }

      typename Traits::EntitySet _es;
      DFBasis dfBasis_;
      std::shared_ptr<FEM const> pfem;
      std::shared_ptr<CE const> _pce;

      mutable std::shared_ptr<Ordering> _ordering;
    };


  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_DUNEFUNCTIONSGRIDFUNCTIONSPACE_HH
