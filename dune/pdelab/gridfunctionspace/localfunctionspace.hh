// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_LOCALFUNCTIONSPACE_HH
#define DUNE_PDELAB_LOCALFUNCTIONSPACE_HH

#include<vector>

#include <dune/common/stdstreams.hh>

#include <dune/localfunctions/common/interfaceswitch.hh>

#include "../common/typetree.hh"
#include "localindex.hh"
#include <dune/pdelab/gridfunctionspace/tags.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup GridFunctionSpace
    //! \ingroup PDELab
    //! \{

    //=======================================
    // local function space base: metaprograms
    //=======================================

    struct gfs_to_lfs {};

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
          child.global = lfs.global;
        }
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

        ComputeSizeVisitor(const Entity& entity)
          : e(entity)
          , offset(0)
        {}

        const Entity& e;
        std::size_t offset;

      };


      template<typename Entity, typename SizeType>
      struct FillIndicesVisitor
        : public TypeTree::TreeVisitor
        , public TypeTree::DynamicTraversal
      {

        template<typename Node, typename TreePath>
        void leaf(Node& node, TreePath treePath)
        {
          _global.resize(node.n);
          node.pgfs->globalIndices(*(node.pfe),e,_global); // get global indices for this finite element
          for (std::size_t i=0; i<node.n; i++) (*node.global)[node.offset+i]=_global[i];
        }

        template<typename Node, typename Child, typename TreePath, typename ChildIndex>
        void afterChild(const Node& node, const Child& child, TreePath treePath, ChildIndex childIndex)
        {
          for (std::size_t i = 0; i<child.n; ++i)
            (*node.global)[child.offset+i] = node.pgfs->subMap(childIndex,(*node.global)[child.offset+i]);
        }

        FillIndicesVisitor(const Entity& entity, std::size_t maxLocalSize)
          : e(entity)
        {
          _global.reserve(maxLocalSize);
        }

        const Entity& e;
        std::vector<SizeType> _global;

      };

    } // end empty namespace

    //=======================================
    // local function space base: base class
    //=======================================

    //! traits mapping global function space information to local function space
    template<typename GFS>
    struct LocalFunctionSpaceBaseTraits
    {
      //! \brief the grid view where grid function is defined upon
      typedef GFS GridFunctionSpaceType;

      //! \brief Type to store indices from Backend
      typedef typename GFS::Traits::GridViewType GridViewType;

      //! \brief Type of codim 0 entity in the grid
      typedef typename GridViewType::Traits::template Codim<0>::Entity Element;

      //! \brief Type to store indices from Backend
      typedef typename GFS::Traits::SizeType SizeType;

      //! \brief Type of container to store indices
      typedef typename std::vector<SizeType> IndexContainer;

    };

    template <typename GFS>
    class LocalFunctionSpaceBaseNode
    {
      typedef typename GFS::Traits::BackendType B;

      template<typename>
      friend struct PropagateGlobalStorageVisitor;

      template<typename>
      friend struct ComputeSizeVisitor;

      template<typename,typename>
      friend struct FillIndicesVisitor;

    public:
      typedef LocalFunctionSpaceBaseTraits<GFS> Traits;

      //! \brief construct from global function space
      LocalFunctionSpaceBaseNode (shared_ptr<const GFS> gfs) :
        pgfs(gfs), global_storage(gfs->maxLocalSize()), global(&global_storage), n(0)
      {}

      //! \brief get current size
      typename Traits::IndexContainer::size_type size () const
      {
        return n;
      }

      //! \brief get maximum possible size (which is maxLocalSize from grid function space)
      typename Traits::IndexContainer::size_type maxSize () const
      {
        return pgfs->maxLocalSize();
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
        return global->size();
      }

      //! \brief map index in this local function space to root local function space
      typename Traits::IndexContainer::size_type localIndex (typename Traits::IndexContainer::size_type index) const
      {
        return offset+index;
      }

      //! \brief map index in this local function space to global index space
      typename Traits::SizeType globalIndex (typename Traits::IndexContainer::size_type index) const
      {
        return (*global)[offset + index];
      }

      //! \brief extract coefficients for one element from container
      template<typename GC, typename LC>
      void vread (const GC& globalcontainer, LC& localcontainer) const
      {
        // assert(&global_storage == global); // make sure we call this method only on the root node!
        localcontainer.resize(n);
        for (typename Traits::IndexContainer::size_type k=0; k<n; ++k)
          localcontainer[typename LC::size_type(k)] = B::access(globalcontainer,(*global)[offset + k]);
      }

      //! \brief write back coefficients for one element to container
      template<typename GC, typename LC>
      void vwrite (const LC& localcontainer, GC& globalcontainer) const
      {
        // assert(&global_storage == global); // make sure we call this method only on the root node!
        for (typename Traits::IndexContainer::size_type k=0; k<n; ++k)
          B::access(globalcontainer,(*global)[offset + k]) = localcontainer[typename LC::size_type(k)];
      }

      //! \brief add coefficients for one element to container
      template<typename GC, typename LC>
      void vadd (const LC& localcontainer, GC& globalcontainer) const
      {
        // assert(&global_storage == global); // make sure we call this method only on the root node!
        for (typename Traits::IndexContainer::size_type k=0; k<n; ++k)
          B::access(globalcontainer,(*global)[offset + k]) += localcontainer[typename LC::size_type(k)];
      }

      //! \brief print debug information about this local function space
      void debug () const
      {
        std::cout << n << " indices = (";
        for (typename Traits::IndexContainer::size_type k=0; k<n; k++)
          std::cout << (*global)[offset + k] << " ";
        std::cout << ")" << std::endl;
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
        TypeTree::applyToTree(node,PropagateGlobalStorageVisitor<>());
      }

      shared_ptr<GFS const> pgfs;
      typename Traits::IndexContainer global_storage;
      typename Traits::IndexContainer* global;
      typename Traits::IndexContainer::size_type n;
      typename Traits::IndexContainer::size_type offset;
    };


    template <typename GFS>
    template <typename NodeType>
    void LocalFunctionSpaceBaseNode<GFS>::bind (NodeType& node,
      const typename LocalFunctionSpaceBaseNode<GFS>::Traits::Element& e)
    {
      typedef typename LocalFunctionSpaceBaseNode<GFS>::Traits::Element Element;
      assert(&node == this);

      // compute sizes
      ComputeSizeVisitor<Element> csv(e);
      TypeTree::applyToTree(node,csv);

      global_storage.resize(node.n);

      // initialize iterators and fill indices
      FillIndicesVisitor<Element,typename Traits::IndexContainer::size_type> fiv(e,node.maxSize());
      TypeTree::applyToTree(node,fiv);

      // apply upMap
      for (typename Traits::IndexContainer::size_type i=0; i<n; ++i)
        global_storage[i] = pgfs->upMap(global_storage[i]);
    }

    //=======================================
    // local function space base: power implementation
    //=======================================

    //! traits for multi component local function space
    template<typename GFS, typename N>
    struct PowerCompositeLocalFunctionSpaceTraits : public LocalFunctionSpaceBaseTraits<GFS>
    {
      //! type of local function space node
      typedef N NodeType;
    };

    // local function space for a power grid function space
    template<typename GFS, typename ChildLFS, std::size_t k>
    class PowerLocalFunctionSpaceNode :
      public LocalFunctionSpaceBaseNode<GFS>,
      public TypeTree::PowerNode<ChildLFS,k>
    {
      typedef LocalFunctionSpaceBaseNode<GFS> BaseT;
      typedef TypeTree::PowerNode<ChildLFS,k> TreeNode;

      template<typename>
      friend struct PropagateGlobalStorageVisitor;

      template<typename>
      friend struct ComputeSizeVisitor;

      template<typename,typename>
      friend struct FillIndicesVisitor;

    public:
      typedef PowerCompositeLocalFunctionSpaceTraits<GFS,PowerLocalFunctionSpaceNode> Traits;

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

    template<typename PowerGridFunctionSpace>
    Dune::PDELab::TypeTree::GenericPowerNodeTransformation<PowerGridFunctionSpace,gfs_to_lfs,PowerLocalFunctionSpaceNode>
    lookupNodeTransformation(PowerGridFunctionSpace* pgfs, gfs_to_lfs* t, PowerGridFunctionSpaceTag tag);


    //=======================================
    // local function space base: composite implementation
    //=======================================

    // local function space for a power grid function space
    template<typename GFS,DUNE_TYPETREE_COMPOSITENODE_TEMPLATE_CHILDREN>
    class CompositeLocalFunctionSpaceNode
      : public LocalFunctionSpaceBaseNode<GFS>
      , public DUNE_TYPETREE_COMPOSITENODE_BASETYPE
    {
      typedef LocalFunctionSpaceBaseNode<GFS> BaseT;
      typedef DUNE_TYPETREE_COMPOSITENODE_BASETYPE NodeType;

      template<typename>
      friend struct PropagateGlobalStorageVisitor;

      template<typename>
      friend struct ComputeSizeVisitor;

      template<typename,typename>
      friend struct FillIndicesVisitor;

    public:
      typedef PowerCompositeLocalFunctionSpaceTraits<GFS,CompositeLocalFunctionSpaceNode> Traits;

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
    template<typename CompositeGridFunctionSpace>
    Dune::PDELab::TypeTree::GenericVariadicCompositeNodeTransformation<CompositeGridFunctionSpace,gfs_to_lfs,CompositeLocalFunctionSpaceNode>
    lookupNodeTransformation(CompositeGridFunctionSpace* cgfs, gfs_to_lfs* t, CompositeGridFunctionSpaceTag tag);
#else
    template<typename CompositeGridFunctionSpace>
    Dune::PDELab::TypeTree::GenericCompositeNodeTransformation<CompositeGridFunctionSpace,gfs_to_lfs,CompositeLocalFunctionSpaceNode>
    lookupNodeTransformation(CompositeGridFunctionSpace* cgfs, gfs_to_lfs* t, CompositeGridFunctionSpaceTag tag);
#endif

    //=======================================
    // local function space base: single component implementation
    //=======================================

    //! traits for single component local function space
    template<typename GFS, typename N>
    struct LeafLocalFunctionSpaceTraits : public PowerCompositeLocalFunctionSpaceTraits<GFS,N>
    {
      //! Type of local finite element
      typedef typename GFS::Traits::FiniteElementType FiniteElementType;

      //! \brief Type of constraints engine
      typedef typename GFS::Traits::ConstraintsType ConstraintsType;
    };

    //! single component local function space
    template<typename GFS>
    class LeafLocalFunctionSpaceNode
      : public LocalFunctionSpaceBaseNode<GFS>
      , public TypeTree::LeafNode
    {
      typedef LocalFunctionSpaceBaseNode<GFS> BaseT;

      template<typename>
      friend struct PropagateGlobalStorageVisitor;

      template<typename>
      friend struct ComputeSizeVisitor;

      template<typename,typename>
      friend struct FillIndicesVisitor;

    public:
      typedef LeafLocalFunctionSpaceTraits<GFS,LeafLocalFunctionSpaceNode> Traits;

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

      /** \brief write back coefficients for one element to container */
      template<typename GC, typename LC>
      void mwrite (const LC& lc, GC& gc) const
      {
        // LC and GC are maps of maps
        typedef typename LC::const_iterator local_col_iterator;
        typedef typename LC::value_type::second_type local_row_type;
        typedef typename local_row_type::const_iterator local_row_iterator;
        typedef typename GC::iterator global_col_iterator;
        typedef typename GC::value_type::second_type global_row_type;

        for (local_col_iterator cit=lc.begin(); cit!=lc.end(); ++cit)
          {
            typename Traits::SizeType i = globalIndex(cit->first);
            // insert empty row in global container if necessary
            global_col_iterator gcit = gc.find(i);
            if (gcit==gc.end())
              gc[i] = global_row_type();

            // copy row to global container with transformed indices
            for (local_row_iterator rit=(cit->second).begin(); rit!=(cit->second).end(); ++rit)
              gc[i][i] = rit->second;
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

    template<typename GridFunctionSpace>
    Dune::PDELab::TypeTree::GenericLeafNodeTransformation<GridFunctionSpace,gfs_to_lfs,LeafLocalFunctionSpaceNode<GridFunctionSpace> >
    lookupNodeTransformation(GridFunctionSpace* gfs, gfs_to_lfs* t, LeafGridFunctionSpaceTag tag);

    //=======================================
    // local function facade
    //=======================================

    /**
       \brief Create a local function space from a global function space

       The local function space can be tagged with on of the tags
       defined in localfunctionspacetags.hh. This allows to
       destinguish between trial and test space.

       If no TAG is specified the AnySpaceTag is used, which basicly
       states, that it is not clear, whether this is a trial of a test
       space.
     */
    template <typename GFS, typename TAG=AnySpaceTag>
    class LocalFunctionSpace;

    // tagged version
    template <typename GFS, typename TAG>
    class LocalFunctionSpace :
      public GFS::LocalFunctionSpace
    {
      typedef typename GFS::LocalFunctionSpace BaseT;
      typedef typename BaseT::Traits::IndexContainer::size_type I;
      typedef typename LocalIndexTraits<I,TAG>::LocalIndex LocalIndex;

      template<typename>
      friend struct PropagateGlobalStorageVisitor;

      template<typename>
      friend struct ComputeSizeVisitor;

      template<typename,typename>
      friend struct FillIndicesVisitor;

    public:
      typedef typename BaseT::Traits Traits;

      LocalFunctionSpace(const GFS & gfs) : BaseT(TypeTree::TransformTree<GFS,gfs_to_lfs>::transform(gfs)) { this->setup(*this); }

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
      public GFS::LocalFunctionSpace
    {
      typedef typename GFS::LocalFunctionSpace BaseT;

      template<typename>
      friend struct PropagateGlobalStorageVisitor;

      template<typename>
      friend struct ComputeSizeVisitor;

      template<typename,typename>
      friend struct FillIndicesVisitor;

    public:
      LocalFunctionSpace(const GFS & gfs) : BaseT(TypeTree::TransformTree<GFS,gfs_to_lfs>::transform(gfs)) { this->setup(*this); }
    };

    //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif
