// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_LOCALFUNCTIONSPACE_HH
#define DUNE_PDELAB_LOCALFUNCTIONSPACE_HH

#include<vector>

#include <dune/common/stdstreams.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/localfunctions/common/interfaceswitch.hh>
#include <dune/localfunctions/common/localkey.hh>

#include <dune/pdelab/common/typetree.hh>
#include <dune/pdelab/common/multiindex.hh>
#include <dune/pdelab/gridfunctionspace/tags.hh>
#include <dune/pdelab/gridfunctionspace/localvector.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup GridFunctionSpace
    //! \ingroup PDELab
    //! \{

    //=======================================
    // local function space base: metaprograms
    //=======================================

    template<typename GFS, typename Tag>
    struct _build_dof_index_type
    {
      typedef Dune::PDELab::DOFIndex<std::size_t,TypeTree::TreeInfo<GFS>::depth,2> type;
    };

    template<typename GFS>
    struct _build_dof_index_type<GFS,SingleCodimMapper>
    {
      typedef SimpleDOFIndex<typename GFS::Traits::SizeType> type;
    };


    template<typename GFS>
    struct build_dof_index_type
    {
      typedef typename _build_dof_index_type<GFS,typename GFS::OrderingTag>::type type;
    };

    //! GridFunctionSpace to LocalFunctionSpace transformation.
    /**
     * gfs_to_lfs describes the transformation of a GridFunctionSpace tree to its corresponding
     * LocalFunctionSpace tree and holds any information that may be required for performing
     * the transformation.
     *
     * \warning The exact meaning of the template parameter is an implementation detail
     *          and may change at any time, as the only class that is supposed to instantiate
     *          the transformation is LocalFunctionSpace. Implementors of custom transformation
     *          descriptors should only use information exported by gfs_to_lfs. In particular,
     *          the registration declaration should not make any assumptions on GFS and just
     *          treat it as some kind of opaque parameter type.
     *
     * \tparam GFS  the \b root GridFunctionSpace that the resulting LocalFunctionSpace tree
     *              will be based on.
     */
    template<typename GFS>
    struct gfs_to_lfs {

      //! The MultiIndex type that will be used in the resulting LocalFunctionSpace tree.
      //typedef Dune::PDELab::MultiIndex<std::size_t,TypeTree::TreeInfo<GFS>::depth> MultiIndex;
      typedef typename build_dof_index_type<GFS>::type DOFIndex;

    };

    namespace {

      // the bogus template parameter is necessary to make GCC honor the friend declaration
      // in the LocalFunctionSpace (probably a GCC bug)
      template<typename = int>
      struct PropagateGlobalStorageVisitor
        : public TypeTree::TreeVisitor
        , public TypeTree::DynamicTraversal
      {

        template<typename LFS, typename Child, typename TreePath, typename ChildIndex>
        void beforeChild(const LFS& lfs, Child& child, TreePath treePath, ChildIndex childIndex) const
        {
          child._dof_indices = lfs._dof_indices;
        }
      };

      // This visitor is not used in standard PDELab code, but is necessary for MultiDomain
      // It is defined here due to the necessary friend declarations in the local function spaces.
      // for template parameter see above
      template<typename = int>
      struct ClearSizeVisitor
        : public TypeTree::TreeVisitor
        , public TypeTree::DynamicTraversal
      {

        template<typename Node, typename TreePath>
        void pre(Node& node, TreePath treePath)
        {
          leaf(node,treePath);
        }

        template<typename Node, typename TreePath>
        void leaf(Node& node, TreePath treePath)
        {
          node.offset = offset;
          node.n = 0;
        }

        ClearSizeVisitor(std::size_t offset_)
          : offset(offset_)
        {}

        const std::size_t offset;

      };


      template<typename Entity>
      struct ComputeSizeVisitor
        : public TypeTree::TreeVisitor
        , public TypeTree::DynamicTraversal
      {

        template<typename Node, typename TreePath>
        void pre(Node& node, TreePath treePath)
        {
          node.offset = offset;
        }

        template<typename Node, typename TreePath>
        void post(Node& node, TreePath treePath)
        {
          node.n = offset - node.offset;
        }

        template<typename Node, typename TreePath>
        void leaf(Node& node, TreePath treePath)
        {
          node.offset = offset;
          Node::FESwitch::setStore(node.pfe, node.pgfs->finiteElementMap().find(e));
          node.n = Node::FESwitch::basis(*node.pfe).size();
          offset += node.n;
        }

        ComputeSizeVisitor(const Entity& entity, std::size_t offset = 0)
          : e(entity)
          , offset(offset)
        {}

        const Entity& e;
        std::size_t offset;

      };


      template<typename Entity>
      struct FillIndicesVisitor
        : public TypeTree::TreeVisitor
        , public TypeTree::DynamicTraversal
      {

        template<typename Node, typename TreePath>
        void leaf(Node& node, TreePath treePath)
        {
          // setup DOFIndices for this finite element
          node.dofIndices(e,node._dof_indices->begin()+node.offset,node._dof_indices->begin()+node.offset+node.n);
        }

        template<typename Node, typename Child, typename TreePath, typename ChildIndex>
        void afterChild(const Node& node, const Child& child, TreePath treePath, ChildIndex childIndex)
        {
          for (std::size_t i = 0; i<child.n; ++i)
            {
              // update tree path for the DOFIndices of the child
              (*node._dof_indices)[child.offset+i].treeIndex().push_back(childIndex);
            }
        }

        FillIndicesVisitor(const Entity& entity)
          : e(entity)
        {}

        const Entity& e;
      };

    } // end empty namespace

    //=======================================
    // local function space base: base class
    //=======================================

    //! traits mapping global function space information to local function space
    template<typename GFS, typename DI>
    struct LocalFunctionSpaceBaseTraits
    {
      //! \brief Type of the underlying grid function space
      typedef GFS GridFunctionSpaceType;

      //! \brief Type of the underlying grid function space
      typedef GFS GridFunctionSpace;

      //! \brief Type of the grid view that the underlying grid function space is defined on.
      typedef typename GFS::Traits::GridViewType GridViewType;

      //! \brief Type of the grid view that the underlying grid function space is defined on.
      typedef typename GFS::Traits::GridViewType GridView;

      //! \brief Type of codim 0 entity in the grid
      typedef typename GridViewType::Traits::template Codim<0>::Entity Element;

      //! \brief Type to store indices from Backend
      typedef typename GFS::Traits::SizeType SizeType;

      //! \brief Type of container to store indices
      typedef typename std::vector<SizeType> IndexContainer;

      //! \brief Type of MultiIndex associated with this LocalFunctionSpace.
      typedef DI DOFIndex;

      //! \brief Type of container to store multiindices.
      typedef typename std::vector<DI> DOFIndexContainer;

    };

    template <typename GFS, typename DOFIndex>
    class LocalFunctionSpaceBaseNode
    {
      typedef typename GFS::Traits::Backend B;

      template<typename>
      friend struct PropagateGlobalStorageVisitor;

      template<typename>
      friend struct ComputeSizeVisitor;

      template<typename>
      friend struct FillIndicesVisitor;

      template<typename LFS, typename C, typename Tag>
      friend class LFSContainerIndexCacheBase;

    public:
      typedef LocalFunctionSpaceBaseTraits<GFS,DOFIndex> Traits;

      //! \brief construct from global function space
      LocalFunctionSpaceBaseNode (shared_ptr<const GFS> gfs)
        : pgfs(gfs)
        , _dof_index_storage()
        , _dof_indices(&_dof_index_storage)
        , n(0)
      {}

      //! \brief get current size
      typename Traits::IndexContainer::size_type size () const
      {
        return n;
      }

      //! \brief get maximum possible size (which is maxLocalSize from grid function space)
      typename Traits::IndexContainer::size_type maxSize () const
      {
        // _dof_indices is always as large as the max local size of the root GFS
        return _dof_indices->size();
      }

      //! \brief get size of an appropriate local vector object
      /**
         this is the number of dofs of the complete local function
         space tree, i.e. the size() of the root node. The local
         vector objects must always have this size and the localIndex
         method maps into the range [0,localVectorSize()[
       */
      typename Traits::IndexContainer::size_type localVectorSize () const
      {
        return _dof_indices->size();
      }

      //! \brief map index in this local function space to root local function space
      typename Traits::IndexContainer::size_type localIndex (typename Traits::IndexContainer::size_type index) const
      {
        return offset+index;
      }

      //! \brief Maps given index in this local function space to its corresponding global MultiIndex.
      /**
       * \param index  The local index value from the range 0,...,size()-1
       * \returns      A const reference to the associated, globally unique MultiIndex. Note that the returned
       *               object may (and must) be copied if it needs to be stored beyond the time of the next
       *               call to bind() on this LocalFunctionSpace (e.g. when the MultiIndex is used as a DOF
       *               identifier in a constraints container).
       */
      const typename Traits::DOFIndex& dofIndex(typename Traits::IndexContainer::size_type index) const
      {
        return (*_dof_indices)[offset + index];
      }

      //! \brief print debug information about this local function space
      void debug () const
      {
        std::cout << n << " indices = (";
        for (typename Traits::IndexContainer::size_type k=0; k<n; k++)
          std::cout << (*_dof_indices)[localIndex(k)] << " ";
        std::cout << ")" << std::endl;
      }

      //! Returns the GridFunctionSpace underlying this LocalFunctionSpace.
      const GFS& gridFunctionSpace() const
      {
        return *pgfs;
      }

    protected:
      //! \brief bind local function space to entity
      /**

         This is a generic implementation of the bind function. It is
         parametrized with the NodeType, which the type of the derived
         LocalFunctionSpaceNode. Handing the NodeType as a parammeter
         avoid the need for the CRTP construct, but all derived
         classes have to add a method bind, which forward to this
         method.

         \param node reference to the derived node, the address must be the same as this
         \param e entity to bind to
       */
      template<typename NodeType>
      void bind (NodeType& node, const typename Traits::Element& e);

      template<typename NodeType>
      void setup(NodeType& node)
      {
        _dof_index_storage.resize(gridFunctionSpace().ordering().maxLocalSize());
        TypeTree::applyToTree(node,PropagateGlobalStorageVisitor<>());
      }

      shared_ptr<GFS const> pgfs;
      typename Traits::DOFIndexContainer _dof_index_storage;
      typename Traits::DOFIndexContainer* _dof_indices;
      typename Traits::IndexContainer::size_type n;
      typename Traits::IndexContainer::size_type offset;
    };


    template <typename GFS, typename DOFIndex>
    template <typename NodeType>
    void LocalFunctionSpaceBaseNode<GFS,DOFIndex>::bind (NodeType& node,
         const typename LocalFunctionSpaceBaseNode<GFS,DOFIndex>::Traits::Element& e)
    {
      typedef typename LocalFunctionSpaceBaseNode<GFS,DOFIndex>::Traits::Element Element;
      assert(&node == this);

      // compute sizes
      ComputeSizeVisitor<Element> csv(e);
      TypeTree::applyToTree(node,csv);


      // initialize iterators and fill indices
      FillIndicesVisitor<Element> fiv(e);
      TypeTree::applyToTree(node,fiv);

    }

    //=======================================
    // local function space base: power implementation
    //=======================================

    //! traits for multi component local function space
    template<typename GFS, typename DOFIndex, typename N>
    struct PowerCompositeLocalFunctionSpaceTraits : public LocalFunctionSpaceBaseTraits<GFS,DOFIndex>
    {
      //! type of local function space node
      typedef N NodeType;
    };

    // local function space for a power grid function space
    template<typename GFS, typename DOFIndex, typename ChildLFS, std::size_t k>
    class PowerLocalFunctionSpaceNode :
      public LocalFunctionSpaceBaseNode<GFS,DOFIndex>,
      public TypeTree::PowerNode<ChildLFS,k>
    {
      typedef LocalFunctionSpaceBaseNode<GFS,DOFIndex> BaseT;
      typedef TypeTree::PowerNode<ChildLFS,k> TreeNode;

      template<typename>
      friend struct PropagateGlobalStorageVisitor;

      template<typename>
      friend struct ClearSizeVisitor;

      template<typename>
      friend struct ComputeSizeVisitor;

      template<typename>
      friend struct FillIndicesVisitor;

    public:
      typedef PowerCompositeLocalFunctionSpaceTraits<GFS,DOFIndex,PowerLocalFunctionSpaceNode> Traits;

      typedef PowerLocalFunctionSpaceTag ImplementationTag;

      //! \brief initialize with grid function space
      template<typename Transformation>
      PowerLocalFunctionSpaceNode (shared_ptr<const GFS> gfs,
                                   const Transformation& t,
                                   const array<shared_ptr<ChildLFS>,k>& children)
        : BaseT(gfs)
        , TreeNode(children)
      {}

      template<typename Transformation>
      PowerLocalFunctionSpaceNode (const GFS& gfs,
                                   const Transformation& t,
                                   const array<shared_ptr<ChildLFS>,k>& children)
        : BaseT(stackobject_to_shared_ptr(gfs))
        , TreeNode(children)
      {}

      //! \brief bind local function space to entity
      void bind (const typename Traits::Element& e)
      {
        // call method on base class, this avoid the barton neckman trick
        BaseT::bind(*this,e);
      }

    };


    // transformation template, we need a custom template in order to inject the DOFIndex type into the LocalFunctionSpace
    template<typename SourceNode, typename Transformation>
    struct power_gfs_to_lfs_template
    {
      template<typename TC>
      struct result
      {
        typedef PowerLocalFunctionSpaceNode<SourceNode,typename Transformation::DOFIndex,TC,SourceNode::CHILDREN> type;
      };
    };

    // register PowerGFS -> LocalFunctionSpace transformation
    template<typename PowerGridFunctionSpace, typename Params>
    Dune::PDELab::TypeTree::TemplatizedGenericPowerNodeTransformation<
      PowerGridFunctionSpace,
      gfs_to_lfs<Params>,
      power_gfs_to_lfs_template<PowerGridFunctionSpace,gfs_to_lfs<Params> >::template result
      >
    lookupNodeTransformation(PowerGridFunctionSpace* pgfs, gfs_to_lfs<Params>* t, PowerGridFunctionSpaceTag tag);


    //=======================================
    // local function space base: composite implementation
    //=======================================

    // local function space for a power grid function space
    template<typename GFS, typename DOFIndex, DUNE_TYPETREE_COMPOSITENODE_TEMPLATE_CHILDREN>
    class CompositeLocalFunctionSpaceNode
      : public LocalFunctionSpaceBaseNode<GFS,DOFIndex>
      , public DUNE_TYPETREE_COMPOSITENODE_BASETYPE
    {
      typedef LocalFunctionSpaceBaseNode<GFS,DOFIndex> BaseT;
      typedef DUNE_TYPETREE_COMPOSITENODE_BASETYPE NodeType;

      template<typename>
      friend struct PropagateGlobalStorageVisitor;

      template<typename>
      friend struct ClearSizeVisitor;

      template<typename>
      friend struct ComputeSizeVisitor;

      template<typename>
      friend struct FillIndicesVisitor;

    public:
      typedef PowerCompositeLocalFunctionSpaceTraits<GFS,DOFIndex,CompositeLocalFunctionSpaceNode> Traits;

      typedef CompositeLocalFunctionSpaceTag ImplementationTag;

      template<typename Transformation>
      CompositeLocalFunctionSpaceNode (shared_ptr<const GFS> gfs,
                                       const Transformation& t,
                                       DUNE_TYPETREE_COMPOSITENODE_STORAGE_CONSTRUCTOR_SIGNATURE)
        : BaseT(gfs)
        , NodeType(DUNE_TYPETREE_COMPOSITENODE_CHILDVARIABLES)
      {}

      template<typename Transformation>
      CompositeLocalFunctionSpaceNode (const GFS& gfs,
                                       const Transformation& t,
                                       DUNE_TYPETREE_COMPOSITENODE_STORAGE_CONSTRUCTOR_SIGNATURE)
        : BaseT(stackobject_to_shared_ptr(gfs))
        , NodeType(DUNE_TYPETREE_COMPOSITENODE_CHILDVARIABLES)
      {}

      //! \brief bind local function space to entity
      void bind (const typename Traits::Element& e)
      {
        // call method on base class, this avoid the barton neckman trick
        BaseT::bind(*this,e);
      }

    };

#if HAVE_VARIADIC_TEMPLATES

    // transformation template, we need a custom template in order to inject the MultiIndex type into the LocalFunctionSpace
    template<typename SourceNode, typename Transformation>
    struct variadic_composite_gfs_to_lfs_template
    {
      template<typename... TC>
      struct result
      {
        typedef CompositeLocalFunctionSpaceNode<SourceNode,typename Transformation::DOFIndex,TC...> type;
      };
    };

    // register CompositeGFS -> LocalFunctionSpace transformation (variadic version)
    template<typename CompositeGridFunctionSpace, typename Params>
    Dune::PDELab::TypeTree::TemplatizedGenericVariadicCompositeNodeTransformation<
      CompositeGridFunctionSpace,
      gfs_to_lfs<Params>,
      variadic_composite_gfs_to_lfs_template<CompositeGridFunctionSpace,gfs_to_lfs<Params> >::template result
      >
    lookupNodeTransformation(CompositeGridFunctionSpace* cgfs, gfs_to_lfs<Params>* t, CompositeGridFunctionSpaceTag tag);

#else

    // transformation template, we need a custom template in order to inject the MultiIndex type into the LocalFunctionSpace
    template<typename SourceNode, typename Transformation>
    struct composite_gfs_to_lfs_template
    {
      template<typename TC0,
               typename TC1,
               typename TC2,
               typename TC3,
               typename TC4,
               typename TC5,
               typename TC6,
               typename TC7,
               typename TC8,
               typename TC9>
      struct result
      {
        typedef CompositeLocalFunctionSpaceNode<SourceNode,typename Transformation::DOFIndex,TC0,TC1,TC2,TC3,TC4,TC5,TC6,TC7,TC8,TC9> type;
      };
    };

    // register CompositeGFS -> LocalFunctionSpace transformation (non-variadic version)
    template<typename CompositeGridFunctionSpace, typename Params>
    Dune::PDELab::TypeTree::TemplatizedGenericCompositeNodeTransformation<
      CompositeGridFunctionSpace,
      gfs_to_lfs<Params>,
      composite_gfs_to_lfs_template<CompositeGridFunctionSpace,gfs_to_lfs<Params> >::template result
      >
    lookupNodeTransformation(CompositeGridFunctionSpace* cgfs, gfs_to_lfs<Params>* t, CompositeGridFunctionSpaceTag tag);

#endif

    //=======================================
    // local function space base: single component implementation
    //=======================================

    //! traits for single component local function space
    template<typename GFS, typename DOFIndex, typename N>
    struct LeafLocalFunctionSpaceTraits : public PowerCompositeLocalFunctionSpaceTraits<GFS,DOFIndex,N>
    {
      //! Type of local finite element
      typedef typename GFS::Traits::FiniteElementType FiniteElementType;

      typedef typename GFS::Traits::FiniteElementType FiniteElement;

      //! \brief Type of constraints engine
      typedef typename GFS::Traits::ConstraintsType ConstraintsType;

      typedef typename GFS::Traits::ConstraintsType Constraints;

    };

    //! single component local function space
    template<typename GFS, typename DOFIndex>
    class LeafLocalFunctionSpaceNode
      : public LocalFunctionSpaceBaseNode<GFS,DOFIndex>
      , public TypeTree::LeafNode
    {
      typedef LocalFunctionSpaceBaseNode<GFS,DOFIndex> BaseT;

      template<typename>
      friend struct PropagateGlobalStorageVisitor;

      template<typename>
      friend struct ClearSizeVisitor;

      template<typename>
      friend struct ComputeSizeVisitor;

      template<typename>
      friend struct FillIndicesVisitor;

    public:
      typedef LeafLocalFunctionSpaceTraits<GFS,DOFIndex,LeafLocalFunctionSpaceNode> Traits;

      typedef LeafLocalFunctionSpaceTag ImplementationTag;

    private:
      typedef FiniteElementInterfaceSwitch<
      typename Traits::FiniteElementType
      > FESwitch;

    public:

      //! \brief initialize with grid function space
      template<typename Transformation>
      LeafLocalFunctionSpaceNode (shared_ptr<const GFS> gfs, const Transformation& t)
        : BaseT(gfs)
      {
      }

      template<typename Transformation>
      LeafLocalFunctionSpaceNode (const GFS& gfs, const Transformation& t)
        : BaseT(stackobject_to_shared_ptr(gfs))
      {
      }

      //! get finite element
      const typename Traits::FiniteElementType& finiteElement () const
      {
        return *pfe;
      }

      //! \brief get local finite element
      const typename Traits::FiniteElementType& localFiniteElement () const
        DUNE_DEPRECATED
      {
        return *pfe;
      }

      //! \brief get constraints engine
      const typename Traits::ConstraintsType& constraints () const
      {
        return this->pgfs->constraints();
      }

      //! Calculates the multiindices associated with the given entity.
      template<typename Entity, typename DOFIndexIterator>
      void dofIndices(const Entity& e, DOFIndexIterator it, DOFIndexIterator endit)
      {
        // get layout of entity
        const typename FESwitch::Coefficients &coeffs =
          FESwitch::coefficients(*pfe);

        typedef typename GFS::Traits::GridViewType GV;
        GV gv = this->gridFunctionSpace().gridView();

        const Dune::GenericReferenceElement<double,GV::Grid::dimension>& refEl =
          Dune::GenericReferenceElements<double,GV::Grid::dimension>::general(this->pfe->type());

        for (std::size_t i = 0; i < std::size_t(coeffs.size()); ++i, ++it)
          {
            // get geometry type of subentity
            Dune::GeometryType gt = refEl.type(coeffs.localKey(i).subEntity(),
                                              coeffs.localKey(i).codim());

            // evaluate consecutive index of subentity
            typename GV::IndexSet::IndexType index = gv.indexSet().subIndex(e,
                                                                            coeffs.localKey(i).subEntity(),
                                                                            coeffs.localKey(i).codim());

            // store data
            GFS::Ordering::Traits::DOFIndexAccessor::store(*it,gt,index,coeffs.localKey(i).index());

            // make sure we don't write past the end of the iterator range
            assert(it != endit);
          }
      }


      template<typename GC, typename LC>
      void insert_constraints (const LC& lc, GC& gc) const
      {
        // LC and GC are maps of maps
        typedef typename LC::const_iterator local_col_iterator;
        typedef typename LC::value_type::second_type::const_iterator local_row_iterator;
        typedef typename GC::iterator global_col_iterator;
        typedef typename GC::value_type::second_type global_row_type;

        for (local_col_iterator cit=lc.begin(); cit!=lc.end(); ++cit)
          {

            // look up entry in global map, if not found, insert an empty one.
            global_col_iterator gcit = gc.insert(std::make_pair(std::ref(this->dofIndex(cit->first)),global_row_type())).first;

            // copy row to global container with transformed indices
            for (local_row_iterator rit=(cit->second).begin(); rit!=(cit->second).end(); ++rit)
              gcit->second[this->dofIndex(rit->first)] = rit->second;
          }
      }

      //! \brief bind local function space to entity
      void bind (const typename Traits::Element& e)
      {
        // call method on base class, this avoid the barton neckman trick
        BaseT::bind(*this,e);
      }

    private:
      typename FESwitch::Store pfe;
    };

    // Register LeafGFS -> LocalFunctionSpace transformation
    template<typename GridFunctionSpace, typename Params>
    Dune::PDELab::TypeTree::GenericLeafNodeTransformation<
      GridFunctionSpace,
      gfs_to_lfs<Params>,
      LeafLocalFunctionSpaceNode<GridFunctionSpace,typename gfs_to_lfs<Params>::DOFIndex>
      >
    lookupNodeTransformation(GridFunctionSpace* gfs, gfs_to_lfs<Params>* t, LeafGridFunctionSpaceTag tag);

    //=======================================
    // local function facade
    //=======================================

    template <typename GFS, typename TAG=AnySpaceTag>
    class LocalFunctionSpace;

    /**
       \brief Create a local function space from a global function space

       The local function space can be tagged with on of the tags
       defined in localfunctionspacetags.hh. This allows to
       destinguish between trial and test space.

       If no TAG is specified the AnySpaceTag is used, which basicly
       states, that it is not clear, whether this is a trial of a test
       space.

       \extends LocalFunctionSpaceBaseNode
     */
    template <typename GFS, typename TAG>
    class LocalFunctionSpace :
      public Dune::PDELab::TypeTree::TransformTree<GFS,gfs_to_lfs<GFS> >::Type
    {
      typedef typename Dune::PDELab::TypeTree::TransformTree<GFS,gfs_to_lfs<GFS> >::Type BaseT;
      typedef typename BaseT::Traits::IndexContainer::size_type I;
      typedef typename BaseT::Traits::IndexContainer::size_type LocalIndex;

      template<typename>
      friend struct PropagateGlobalStorageVisitor;

      template<typename>
      friend struct ClearSizeVisitor;

      template<typename>
      friend struct ComputeSizeVisitor;

      template<typename>
      friend struct FillIndicesVisitor;

    public:
      typedef typename BaseT::Traits Traits;

      LocalFunctionSpace(const GFS & gfs)
        : BaseT(TypeTree::TransformTree<GFS,gfs_to_lfs<GFS> >::transform(gfs))
      {
        this->setup(*this);
      }

      LocalFunctionSpace(const LocalFunctionSpace & lfs)
        : BaseT(lfs)
      {
        // We need to reset the DOFIndex storage pointers in the new LFS tree,
        // as they are still pointing to the _dof_index_storage of the
        // old tree.
        this->_dof_indices = &(this->_dof_index_storage);
        this->setup(*this);
      }

      LocalIndex localIndex (typename Traits::IndexContainer::size_type index) const
      {
        return LocalIndex(BaseT::localIndex(index));
      }

    private:
      // we don't support getChild yet, so let's hide it!
      template<int i>
      void getChild () const;
      template<int i>
      void child () const;
    };

    // specialization for AnySpaceTag
    template <typename GFS>
    class LocalFunctionSpace<GFS, AnySpaceTag> :
      public Dune::PDELab::TypeTree::TransformTree<GFS,gfs_to_lfs<GFS> >::Type
    {
      typedef typename Dune::PDELab::TypeTree::TransformTree<GFS,gfs_to_lfs<GFS> >::Type BaseT;

      template<typename>
      friend struct PropagateGlobalStorageVisitor;

      template<typename>
      friend struct ClearSizeVisitor;

      template<typename>
      friend struct ComputeSizeVisitor;

      template<typename>
      friend struct FillIndicesVisitor;

    public:
      LocalFunctionSpace(const GFS & gfs)
        : BaseT(TypeTree::TransformTree<GFS,gfs_to_lfs<GFS> >::transform(gfs))
      {
        this->_dof_indices = &(this->_dof_index_storage);
        this->setup(*this);
      }

      LocalFunctionSpace(const LocalFunctionSpace & lfs)
        : BaseT(lfs)
      {
        // We need to reset the DOFIndex storage pointers in the new LFS tree,
        // as they are still pointing to the _dof_index_storage of the
        // old tree.
        this->_dof_indices = &(this->_dof_index_storage);
        this->setup(*this);
      }

    };

    //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif
