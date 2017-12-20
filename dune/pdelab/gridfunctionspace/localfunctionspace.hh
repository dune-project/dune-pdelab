// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_LOCALFUNCTIONSPACE_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_LOCALFUNCTIONSPACE_HH

#include<vector>

#include <dune/common/stdstreams.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/localfunctions/common/interfaceswitch.hh>
#include <dune/localfunctions/common/localkey.hh>

#include <dune/typetree/typetree.hh>

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


      template<typename Entity, bool fast>
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
          if (fast)
            {
              node.pfe = nullptr;
              node.n = node.pgfs->finiteElementMap().maxLocalSize();
              Node::FESwitch::setStore(node.pfe, node.pgfs->finiteElementMap().find(e));
            }
          else
            {
              Node::FESwitch::setStore(node.pfe, node.pgfs->finiteElementMap().find(e));
              node.n = Node::FESwitch::basis(*node.pfe).size();
            }
          offset += node.n;
        }

        ComputeSizeVisitor(const Entity& entity, std::size_t offset_ = 0)
          : e(entity)
          , offset(offset_)
        {}

        const Entity& e;
        std::size_t offset;

      };


      template<typename Entity, bool fast>
      struct FillIndicesVisitor
        : public TypeTree::TreeVisitor
        , public TypeTree::DynamicTraversal
      {

        template<typename Node, typename TreePath>
        void leaf(Node& node, TreePath treePath)
        {
          // setup DOFIndices for this finite element
          node.dofIndices(e,node._dof_indices->begin()+node.offset,node._dof_indices->begin()+node.offset+node.n,std::integral_constant<bool,fast>{});
        }

        template<typename Node, typename Child, typename TreePath, typename ChildIndex>
        void afterChild(const Node& node, const Child& child, TreePath treePath, ChildIndex childIndex)
        {
          // Just skip the entire function space structure handling in fast mode
          // This **really** breaks any attempt at using the DOFIndex for anything other
          // than as input to the FastDGGridOperator machine.
          // You have been warned!
          if (not fast)
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

      template<typename,bool>
      friend struct ComputeSizeVisitor;

      template<typename,bool>
      friend struct FillIndicesVisitor;

      template<typename LFS, typename C, typename Tag, bool>
      friend class LFSIndexCacheBase;

    public:
      typedef LocalFunctionSpaceBaseTraits<GFS,DOFIndex> Traits;

      //! \brief construct from global function space
      LocalFunctionSpaceBaseNode (std::shared_ptr<const GFS> gfs)
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

      std::size_t subSpaceDepth() const
      {
        return 0;
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

    public:
      template<typename NodeType>
      void setup(NodeType& node)
      {
        _dof_index_storage.resize(gridFunctionSpace().ordering().maxLocalSize());
        TypeTree::applyToTree(node,PropagateGlobalStorageVisitor<>());
      }

      std::shared_ptr<GFS const> pgfs;
      typename Traits::DOFIndexContainer _dof_index_storage;
      typename Traits::DOFIndexContainer* _dof_indices;
      typename Traits::IndexContainer::size_type n;
      typename Traits::IndexContainer::size_type offset;
    };

    //! traits for local function space on a gridview
    template<typename GFS, typename DOFIndex>
    struct GridViewLocalFunctionSpaceBaseTraits : public LocalFunctionSpaceBaseTraits<GFS,DOFIndex>
    {
      //! \brief Type of the grid view that the underlying grid function space is defined on.
      typedef typename GFS::Traits::GridViewType GridViewType;

      //! \brief Type of the grid view that the underlying grid function space is defined on.
      typedef typename GFS::Traits::GridViewType GridView;

      using EntitySet = typename GFS::Traits::EntitySet;

      //! \brief Type of codim 0 entity in the grid
      using Element = typename EntitySet::Element;
    };

    template <typename GFS, typename DOFIndex>
    class GridViewLocalFunctionSpaceBaseNode :
      public LocalFunctionSpaceBaseNode<GFS,DOFIndex>
    {
      typedef typename GFS::Traits::Backend B;
      typedef LocalFunctionSpaceBaseNode<GFS,DOFIndex> BaseT;

    public:
      typedef GridViewLocalFunctionSpaceBaseTraits<GFS,DOFIndex> Traits;

      //! \brief construct from global function space
      GridViewLocalFunctionSpaceBaseNode (std::shared_ptr<const GFS> gfs)
        : BaseT(gfs)
      {}

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
      template<typename NodeType, bool fast = false>
      void bind (NodeType& node, const typename Traits::Element& e, std::integral_constant<bool,fast> = std::integral_constant<bool,fast>{});
    };

    template <typename GFS, typename DOFIndex>
    template <typename NodeType, bool fast>
    void GridViewLocalFunctionSpaceBaseNode<GFS,DOFIndex>::bind (NodeType& node,
                                                                 const typename GridViewLocalFunctionSpaceBaseNode<GFS,DOFIndex>::Traits::Element& e,
                                                                 std::integral_constant<bool,fast>)
    {
      typedef typename GridViewLocalFunctionSpaceBaseNode<GFS,DOFIndex>::Traits::Element Element;
      assert(&node == this);

      // compute sizes
      ComputeSizeVisitor<Element,fast> csv(e);
      TypeTree::applyToTree(node,csv);


      // initialize iterators and fill indices
      FillIndicesVisitor<Element,fast> fiv(e);
      TypeTree::applyToTree(node,fiv);
    }

    //=======================================
    // local function space base: power implementation
    //=======================================

    //! traits for multi component local function space
    template<typename GFS, typename DOFIndex, typename N>
    struct PowerCompositeLocalFunctionSpaceTraits : public GridViewLocalFunctionSpaceBaseTraits<GFS,DOFIndex>
    {
      //! type of local function space node
      typedef N NodeType;
    };

    // local function space for a power grid function space
    template<typename GFS, typename DOFIndex, typename ChildLFS, std::size_t k>
    class PowerLocalFunctionSpaceNode :
      public GridViewLocalFunctionSpaceBaseNode<GFS,DOFIndex>,
      public TypeTree::PowerNode<ChildLFS,k>
    {
      typedef GridViewLocalFunctionSpaceBaseNode<GFS,DOFIndex> BaseT;
      typedef TypeTree::PowerNode<ChildLFS,k> TreeNode;

      template<typename>
      friend struct PropagateGlobalStorageVisitor;

      template<typename>
      friend struct ClearSizeVisitor;

      template<typename,bool>
      friend struct ComputeSizeVisitor;

      template<typename,bool>
      friend struct FillIndicesVisitor;

    public:
      typedef PowerCompositeLocalFunctionSpaceTraits<GFS,DOFIndex,PowerLocalFunctionSpaceNode> Traits;

      typedef PowerLocalFunctionSpaceTag ImplementationTag;

      //! \brief initialize with grid function space
      template<typename Transformation>
      PowerLocalFunctionSpaceNode (std::shared_ptr<const GFS> gfs,
                                   const Transformation& t,
                                   const std::array<std::shared_ptr<ChildLFS>,k>& children)
        : BaseT(gfs)
        , TreeNode(children)
      {}

      template<typename Transformation>
      PowerLocalFunctionSpaceNode (const GFS& gfs,
                                   const Transformation& t,
                                   const std::array<std::shared_ptr<ChildLFS>,k>& children)
        : BaseT(stackobject_to_shared_ptr(gfs))
        , TreeNode(children)
      {}

      //! \brief bind local function space to entity
      template<bool fast = false>
      void bind (const typename Traits::Element& e, std::integral_constant<bool,fast> fast_ = std::integral_constant<bool,fast>{})
      {
        // call method on base class, this avoid the barton neckman trick
        BaseT::bind(*this,e,fast_);
      }

    };


    // transformation template, we need a custom template in order to inject the DOFIndex type into the LocalFunctionSpace
    template<typename SourceNode, typename Transformation>
    struct power_gfs_to_lfs_template
    {
      template<typename TC>
      struct result
      {
        typedef PowerLocalFunctionSpaceNode<SourceNode,typename Transformation::DOFIndex,TC,TypeTree::StaticDegree<SourceNode>::value> type;
      };
    };

    // register PowerGFS -> LocalFunctionSpace transformation
    template<typename PowerGridFunctionSpace, typename Params>
    TypeTree::TemplatizedGenericPowerNodeTransformation<
      PowerGridFunctionSpace,
      gfs_to_lfs<Params>,
      power_gfs_to_lfs_template<PowerGridFunctionSpace,gfs_to_lfs<Params> >::template result
      >
    registerNodeTransformation(PowerGridFunctionSpace* pgfs, gfs_to_lfs<Params>* t, PowerGridFunctionSpaceTag* tag);


    //=======================================
    // local function space base: composite implementation
    //=======================================

    // local function space for a power grid function space
    template<typename GFS, typename DOFIndex, typename... Children>
    class CompositeLocalFunctionSpaceNode
      : public GridViewLocalFunctionSpaceBaseNode<GFS,DOFIndex>
      , public TypeTree::CompositeNode<Children...>
    {
      typedef GridViewLocalFunctionSpaceBaseNode<GFS,DOFIndex> BaseT;
      typedef TypeTree::CompositeNode<Children...> NodeType;

      template<typename>
      friend struct PropagateGlobalStorageVisitor;

      template<typename>
      friend struct ClearSizeVisitor;

      template<typename,bool>
      friend struct ComputeSizeVisitor;

      template<typename,bool>
      friend struct FillIndicesVisitor;

    public:
      typedef PowerCompositeLocalFunctionSpaceTraits<GFS,DOFIndex,CompositeLocalFunctionSpaceNode> Traits;

      typedef CompositeLocalFunctionSpaceTag ImplementationTag;

      template<typename Transformation>
      CompositeLocalFunctionSpaceNode (std::shared_ptr<const GFS> gfs,
                                       const Transformation& t,
                                       std::shared_ptr<Children>... children)
        : BaseT(gfs)
        , NodeType(children...)
      {}

      template<typename Transformation>
      CompositeLocalFunctionSpaceNode (const GFS& gfs,
                                       const Transformation& t,
                                       std::shared_ptr<Children>... children)
        : BaseT(stackobject_to_shared_ptr(gfs))
        , NodeType(children...)
      {}

      //! \brief bind local function space to entity
      template<bool fast = false>
      void bind (const typename Traits::Element& e, std::integral_constant<bool,fast> fast_ = std::integral_constant<bool,fast>{})
      {
        // call method on base class, this avoid the barton neckman trick
        BaseT::bind(*this,e,fast_);
      }

    };

    // transformation template, we need a custom template in order to inject the MultiIndex type into the LocalFunctionSpace
    template<typename SourceNode, typename Transformation>
    struct composite_gfs_to_lfs_template
    {
      template<typename... TC>
      struct result
      {
        typedef CompositeLocalFunctionSpaceNode<SourceNode,typename Transformation::DOFIndex,TC...> type;
      };
    };

    // register CompositeGFS -> LocalFunctionSpace transformation (variadic version)
    template<typename CompositeGridFunctionSpace, typename Params>
    TypeTree::TemplatizedGenericCompositeNodeTransformation<
      CompositeGridFunctionSpace,
      gfs_to_lfs<Params>,
      composite_gfs_to_lfs_template<CompositeGridFunctionSpace,gfs_to_lfs<Params> >::template result
      >
    registerNodeTransformation(CompositeGridFunctionSpace* cgfs, gfs_to_lfs<Params>* t, CompositeGridFunctionSpaceTag* tag);


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
      : public GridViewLocalFunctionSpaceBaseNode<GFS,DOFIndex>
      , public TypeTree::LeafNode
    {
      typedef GridViewLocalFunctionSpaceBaseNode<GFS,DOFIndex> BaseT;

      template<typename>
      friend struct PropagateGlobalStorageVisitor;

      template<typename>
      friend struct ClearSizeVisitor;

      template<typename,bool>
      friend struct ComputeSizeVisitor;

      template<typename,bool>
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
      LeafLocalFunctionSpaceNode (std::shared_ptr<const GFS> gfs, const Transformation& t)
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
        assert(pfe);
        return *pfe;
      }

      //! \brief get constraints engine
      const typename Traits::ConstraintsType& constraints () const
      {
        return this->pgfs->constraints();
      }

      //! Calculates the multiindices associated with the given entity.
      template<typename Entity, typename DOFIndexIterator, bool fast>
      void dofIndices(const Entity& e, DOFIndexIterator it, DOFIndexIterator endit, std::integral_constant<bool,fast>)
      {
        if (fast)
          {
            auto gt = e.type();
            auto index = this->gridFunctionSpace().entitySet().indexSet().index(e);
            GFS::Ordering::Traits::DOFIndexAccessor::store(*it,gt,index,0);
            ++it;
          }
        else
          {
            // get layout of entity
            const typename FESwitch::Coefficients &coeffs =
              FESwitch::coefficients(*pfe);

            using EntitySet = typename GFS::Traits::EntitySet;
            auto es = this->gridFunctionSpace().entitySet();

            auto refEl = Dune::ReferenceElements<double,EntitySet::dimension>::general(this->pfe->type());

            for (std::size_t i = 0; i < std::size_t(coeffs.size()); ++i, ++it)
              {
                // get geometry type of subentity
                auto gt = refEl.type(coeffs.localKey(i).subEntity(),
                  coeffs.localKey(i).codim());

                // evaluate consecutive index of subentity
                auto index = es.indexSet().subIndex(e,
                  coeffs.localKey(i).subEntity(),
                  coeffs.localKey(i).codim());

                // store data
                GFS::Ordering::Traits::DOFIndexAccessor::store(*it,gt,index,coeffs.localKey(i).index());

                // make sure we don't write past the end of the iterator range
                assert(it != endit);
              }
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
      template<bool fast = false>
      void bind (const typename Traits::Element& e, std::integral_constant<bool,fast> fast_ = std::integral_constant<bool,fast>{})
      {
        // call method on base class, this avoid the barton neckman trick
        BaseT::bind(*this,e,fast_);
      }

      //    private:
      typename FESwitch::Store pfe;
    };

    // Register LeafGFS -> LocalFunctionSpace transformation
    template<typename GridFunctionSpace, typename Params>
    TypeTree::GenericLeafNodeTransformation<
      GridFunctionSpace,
      gfs_to_lfs<Params>,
      LeafLocalFunctionSpaceNode<GridFunctionSpace,typename gfs_to_lfs<Params>::DOFIndex>
      >
    registerNodeTransformation(GridFunctionSpace* gfs, gfs_to_lfs<Params>* t, LeafGridFunctionSpaceTag* tag);

    //=======================================
    // local function facade
    //=======================================

    template <typename GFS, typename TAG=AnySpaceTag>
    class LocalFunctionSpace;

    /**
       \brief Create a local function space from a global function space

       The local function space can be tagged with one of the tags
       defined in localfunctionspacetags.hh. This allows to
       destinguish between trial and test space.

       If no TAG is specified the AnySpaceTag is used, which basically
       states, that it is not clear, whether this is a trial of a test
       space.

       \extends LocalFunctionSpaceBaseNode
     */
    template <typename GFS, typename TAG>
    class LocalFunctionSpace :
      public TypeTree::TransformTree<GFS,gfs_to_lfs<GFS> >::Type
    {
      typedef typename TypeTree::TransformTree<GFS,gfs_to_lfs<GFS> >::Type BaseT;
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
    // WARNING: If you modify this class, make sure to also fix the specialization in
    // subspacelocalfunctionspace.hh!
    template <typename GFS>
    class LocalFunctionSpace<GFS, AnySpaceTag> :
      public TypeTree::TransformTree<GFS,gfs_to_lfs<GFS> >::Type
    {
      typedef typename TypeTree::TransformTree<GFS,gfs_to_lfs<GFS> >::Type BaseT;

      template<typename>
      friend struct PropagateGlobalStorageVisitor;

      template<typename>
      friend struct ClearSizeVisitor;

      template<typename,bool>
      friend struct ComputeSizeVisitor;

      template<typename,bool>
      friend struct FillIndicesVisitor;

    public:

      LocalFunctionSpace(const GFS & gfs)
        : BaseT(TypeTree::TransformTree<GFS,gfs_to_lfs<GFS> >::transform(gfs))
      {
        this->_dof_indices = &(this->_dof_index_storage);
        this->setup(*this);
      }

      LocalFunctionSpace(std::shared_ptr<const GFS> pgfs)
        : BaseT(*TypeTree::TransformTree<GFS,gfs_to_lfs<GFS> >::transform_storage(pgfs))
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

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_LOCALFUNCTIONSPACE_HH
