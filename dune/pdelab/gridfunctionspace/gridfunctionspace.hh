// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_GRIDFUNCTIONSPACE_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_GRIDFUNCTIONSPACE_HH

#include <cstddef>
#include <map>
#include <ostream>
#include <set>
#include <vector>
#include <memory>

#include <dune/common/deprecated.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/stdstreams.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/shared_ptr.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/interfaceswitch.hh>
#include <dune/localfunctions/common/localkey.hh>

#include <dune/common/version.hh>
#include <dune/typetree/typetree.hh>

#if DUNE_VERSION_GTE_REV(DUNE_COMMON,2,10,0)
#include <dune-pdelab-config.hh>
#endif

// This alias should be removed after a PDELab 2.7 release.
#if DUNE_VERSION_LT_REV(DUNE_TYPETREE,2,7,1)
namespace Dune {
  namespace TypeTree {
    template<std::size_t... i>
    using StaticTreePath = TreePath<i...>;
  }
}
#endif

#include <dune/pdelab/common/partitionviewentityset.hh>
#include <dune/pdelab/backend/interface.hh>
#include <dune/pdelab/backend/istl/descriptors.hh>
#include <dune/pdelab/constraints/noconstraints.hh>
#include <dune/pdelab/gridfunctionspace/compositegridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/datahandleprovider.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspacebase.hh>
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

#ifndef DOXYGEN

    namespace impl {

      // Helper structs to avoid compilation failures in the
      // backwards compatibility mode where users stuff a
      // GridView into a GridFunctionSpace.
      // In that case, we cannot extract the GridView type from
      // the GridView itself, so we use a std::conditional in the
      // Traits class to pick either one of the following structs
      // and then use the correct class to do the lookup.

      struct _lazy_identity
      {
        template<typename T>
        struct evaluate
        {
          using type = T;
        };
      };

      struct _lazy_extract_gridview
      {
        template<typename T>
        struct evaluate
        {
          using type = typename T::GridView;
        };
      };

      // Returns a GridView, regardless of whether GV_or_ES is a GridView or an EntitySet
      template<typename GV_or_ES>
      using GridView = typename std::conditional<
        isEntitySet<GV_or_ES>::value,
        impl::_lazy_extract_gridview,
        impl::_lazy_identity
        >::type::template evaluate<GV_or_ES>::type;


      // Returns an EntitySet, regardless of whether GV_or_ES is a GridView or an EntitySet
      template<typename GV_or_ES>
      using EntitySet = typename std::conditional<
        isEntitySet<GV_or_ES>::value,
        GV_or_ES,
        AllEntitySet<GV_or_ES>
        >::type;

    }

#endif // DOXYGEN

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

    /** \brief A grid function space.
     *  \details Representation of a function space in a grid. Here, the entity
     *    set.
     *
     *  \tparam ES   Entity Set type (See PartitionViewEntitySet). Respresents
     *               the sub set of entities where the finite element map has
     *               support.
     *  \tparam FEM  Type implementing FiniteElementMapInterface. A map from
     *               entity to local finite element.
     *  \tparam CE   Type for constraints assembler
     *  \tparam B    Backend type
     *  \tparam O    Ordering tag
     */
    template<typename ES, typename FEM, typename CE=NoConstraints,
             typename B=ISTL::VectorBackend<>, typename O=DefaultLeafOrderingTag>
    class GridFunctionSpace
      : public TypeTree::LeafNode
      , public GridFunctionSpaceBase<
                 GridFunctionSpace<ES,FEM,CE,B,O>,
                 GridFunctionSpaceTraits<ES,FEM,CE,B,O>
                 >
      , public GridFunctionOutputParameters
      , public DataHandleProvider<GridFunctionSpace<ES,FEM,CE,B,O> >
    {

      typedef TypeTree::TransformTree<GridFunctionSpace,gfs_to_ordering<GridFunctionSpace> > ordering_transformation;

      template<typename,typename>
      friend class GridFunctionSpaceBase;

    public:
      //! export Traits class
      typedef GridFunctionSpaceTraits<ES,FEM,CE,B,O> Traits;

    private:

      typedef GridFunctionSpaceBase<GridFunctionSpace,Traits> BaseT;

    public:

      typedef typename ES::Traits::template Codim<0>::Entity Element;
      typedef typename ES::Traits::template Codim<0>::Iterator ElementIterator;

      [[deprecated]]
      typedef O SizeTag;

      typedef O OrderingTag;

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
      // Construct from GridView
      // ****************************************************************************************************

      //! constructor
      GridFunctionSpace (const typename Traits::GridView& gridview, const FEM& fem, const CE& ce, const B& backend = B(), const OrderingTag& ordering_tag = OrderingTag())
#if DUNE_PDELAB_WARN_ON_GRIDVIEW_BASED_GFS
        [[deprecated("GridFunctionSpaces now internally use an EntitySet instead of a GridView, please replace the template parameter and the first constructor parameter by an EntitySet.")]]
#endif
        : BaseT(backend,ordering_tag)
        , pfem(stackobject_to_shared_ptr(fem))
        , _pce(stackobject_to_shared_ptr(ce))
      {
        this->setEntitySet(typename Traits::EntitySet{gridview});
      }

      //! constructor
      GridFunctionSpace (const typename Traits::GridView& gridview, const std::shared_ptr<const FEM>& fem, const std::shared_ptr<const CE>& ce, const B& backend = B(), const OrderingTag& ordering_tag = OrderingTag())
#if DUNE_PDELAB_WARN_ON_GRIDVIEW_BASED_GFS
        [[deprecated("GridFunctionSpaces now internally use an EntitySet instead of a GridView, please replace the template parameter and the first constructor parameter by an EntitySet.")]]
#endif
        : BaseT(backend,ordering_tag)
        , pfem(fem)
        , _pce(ce)
      {
        this->setEntitySet(typename Traits::EntitySet{gridview});
      }

      //! constructor
      GridFunctionSpace (const typename Traits::GridView& gridview, const FEM& fem, const B& backend = B(), const OrderingTag& ordering_tag = OrderingTag())
#if DUNE_PDELAB_WARN_ON_GRIDVIEW_BASED_GFS
        [[deprecated("GridFunctionSpaces now internally use an EntitySet instead of a GridView, please replace the template parameter and the first constructor parameter by an EntitySet.")]]
#endif
        : BaseT(backend,ordering_tag)
        , pfem(stackobject_to_shared_ptr(fem))
        , _pce(std::make_shared<CE>())
      {
        this->setEntitySet(typename Traits::EntitySet{gridview});
      }

      //! constructor
      GridFunctionSpace (const typename Traits::GridView& gridview, const std::shared_ptr<const FEM>& fem, const B& backend = B(), const OrderingTag& ordering_tag = OrderingTag())
#if DUNE_PDELAB_WARN_ON_GRIDVIEW_BASED_GFS
        [[deprecated("GridFunctionSpaces now internally use an EntitySet instead of a GridView, please replace the template parameter and the first constructor parameter by an EntitySet.")]]
#endif
        : BaseT(backend,ordering_tag)
        , pfem(fem)
        , _pce(std::make_shared<CE>())
      {
        this->setEntitySet(typename Traits::EntitySet{gridview});
      }


      // ****************************************************************************************************
      // Construct from EntitySet
      // ****************************************************************************************************


      /**
       * @brief Construct a new Grid Function Space object
       * @warning The entity set internals will be modified according to the
       *          finite element map used codimensions
       * @warning This version of the constructor takes a reference on the fem
       *          and the ce. Therefore, these objects shall live longer than
       *          the grid function space!
       *
       * @param entitySet     Copy of an entity set
       * @param fem           Finite Element Map
       * @param ce            Constraints Assembler
       * @param backend       Vector backend
       * @param ordering_tag  Ordering tag
       */
      GridFunctionSpace(typename Traits::EntitySet entitySet, const FEM &fem,
                        const CE &ce, const B &backend = B(),
                        const OrderingTag &ordering_tag = OrderingTag())
          : BaseT(backend, ordering_tag), pfem(stackobject_to_shared_ptr(fem)),
            _pce(stackobject_to_shared_ptr(ce)) {
        this->setEntitySet(std::move(entitySet));
      }

      /**
       * @copybrief GridFunctionSpace
       * @warning The entity set internals will be modified according to the
       *          finite element map used codimensions
       *
       * @param entitySet     Copy of an entity set
       * @param fem           Finite Element Map pointer
       * @param ce            Constraints Assembler pointer
       * @param backend       Vector backend
       * @param ordering_tag  Ordering tag
       */
      GridFunctionSpace(typename Traits::EntitySet entitySet,
                        const std::shared_ptr<const FEM> &fem,
                        const std::shared_ptr<const CE> &ce,
                        const B &backend = B(),
                        const OrderingTag &ordering_tag = OrderingTag())
          : BaseT(backend, ordering_tag)
          , pfem(fem)
          , _pce(ce)
      {
        this->setEntitySet(entitySet);
      }

      /**
       * @copybrief GridFunctionSpace
       * @warning The entity set internals will be modified according to the
       *          finite element map used codimensions
       * @warning This version of the constructor takes a reference on the fem.
       *          Therefore, these objects shall live longer than the grid
       *          function space!
       *
       * @param entitySet     Copy of an entity set
       * @param fem           Finite Element Map
       * @param backend       Vector backend
       * @param ordering_tag  Ordering tag
       */
      GridFunctionSpace(typename Traits::EntitySet entitySet,
                        const FEM &fem,
                        const B &backend = B(),
                        const OrderingTag &ordering_tag = OrderingTag())
          : BaseT(backend, ordering_tag)
          , pfem(stackobject_to_shared_ptr(fem))
          , _pce(std::make_shared<CE>())
      {
        this->setEntitySet(entitySet);
      }

      /**
       * @copybrief GridFunctionSpace
       * @warning The entity set internals will be modified according to the
       *          finite element map used codimensions
       *
       * @param entitySet     Copy of an entity set
       * @param fem           Finite Element Map pointer
       * @param backend       Vector backend
       * @param ordering_tag  Ordering tag
       */
      GridFunctionSpace(typename Traits::EntitySet entitySet,
                        const std::shared_ptr<const FEM> &fem,
                        const B &backend = B(),
                        const OrderingTag &ordering_tag = OrderingTag())
          : BaseT(backend, ordering_tag)
          , pfem(fem)
          , _pce(std::make_shared<CE>())
      {
        this->setEntitySet(entitySet);
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

      std::shared_ptr<FEM const> pfem;
      std::shared_ptr<CE const> _pce;

      mutable std::shared_ptr<Ordering> _ordering;
    };


  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_GRIDFUNCTIONSPACE_HH
