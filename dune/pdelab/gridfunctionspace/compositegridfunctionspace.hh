// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_COMPOSITEGRIDFUNCTIONSPACE_HH
#define DUNE_PDELAB_COMPOSITEGRIDFUNCTIONSPACE_HH

#include "gridfunctionspace.hh"
#include "../common/typetree.hh"

namespace Dune {
  namespace PDELab {

    //=======================================
    // composite grid function space
    //=======================================

    //! \addtogroup GridFunctionSpace grid function space
    //! \ingroup PDELab
    //! \{

    //! Mixin class providing DataHandle support for composite and power GridFunctionSpaces.
    template<typename GridFunctionSpace>
    class PowerCompositeDataHandleProvider
    {

#ifndef DOXYGEN

      //! Functor for dataHandleFixedSize()
      struct DataHandleFixedSize
      {

        template<typename Node, typename TreePath>
        bool operator()(const Node& node, TreePath treePath) const
        {
          return node.dataHandleFixedSize(dim,codim);
        }

        DataHandleFixedSize(int dimension, int codimension)
          : dim(dimension)
          , codim(codimension)
        {}

        const int dim;
        const int codim;

      };

      //! Functor for dataHandleContains()
      struct DataHandleContains
      {

        template<typename Node, typename TreePath>
        bool operator()(const Node& node, TreePath treePath) const
        {
          return node.dataHandleContains(dim,codim);
        }

        DataHandleContains(int dimension, int codimension)
          : dim(dimension)
          , codim(codimension)
        {}

        const int dim;
        const int codim;

      };

      //! Functor for dataHandleSize()
      template<typename Entity>
      struct DataHandleSize
      {

        template<typename Node, typename TreePath>
        bool operator()(const Node& node, TreePath treePath) const
        {
          return node.dataHandleSize(e);
        }

        DataHandleSize(const Entity& entity)
          : e(entity)
        {}

        const Entity& e;

      };

      //! Visitor for retrieving the global DOF indices of a given entity.
      template<typename Entity, typename Container>
      struct DataHandleGlobalIndicesVisitor
        : public TypeTree::TreeVisitor
      {

        template<typename Node, typename TreePath>
        void leaf(const Node& node, TreePath treePath)
        {
          pos += node.dataHandleGlobalIndices(e,g,pos,false);
        }

        template<typename CompositeGFS, typename Child, typename TreePath, typename ChildIndex>
        void beforeChild(CompositeGFS& cgfs, const Child& child, TreePath treePath, ChildIndex childIndex)
        {
          offsets.push_back(pos);
        }

        template<typename CompositeGFS, typename Child, typename TreePath, typename ChildIndex>
        void afterChild(CompositeGFS& cgfs, const Child& child, TreePath treePath, ChildIndex childIndex)
        {
          std::size_t offset = offsets.back();
          offsets.pop_back();
          for (std::size_t i = offset; i < pos; ++i)
            g[i] = cgfs.template subMap<ChildIndex::value>(g[i]);
        }

        DataHandleGlobalIndicesVisitor(const Entity& entity, Container& global)
          : e(entity)
          , g(global)
          , pos(0)
        {
          // a reasonable upper bound for the tree depth - this way, we avoid reallocations
          offsets.reserve(16);
        }

        const Entity& e;
        Container& g;
        std::size_t pos;
        std::vector<std::size_t> offsets;

      };

      const GridFunctionSpace& gfs() const
      {
        return static_cast<const GridFunctionSpace&>(*this);
      }

#endif // DOXYGEN

    public:

      //------------------------------
      // generic data handle interface
      //------------------------------

      //! returns true if data for this codim should be communicated
      bool dataHandleContains (int dim, int codim) const
      {
        return TypeTree::reduceOverLeafs(gfs(),DataHandleContains(dim,codim),std::logical_or<bool>(),false);
      }

      //! returns true if size per entity of given dim and codim is a constant
      bool dataHandleFixedSize (int dim, int codim) const
      {
        return TypeTree::reduceOverLeafs(gfs(),DataHandleFixedSize(dim,codim),std::logical_and<bool>(),true);
      }

      /*! how many objects of type DataType have to be sent for a given entity

        Note: Only the sender side needs to know this size.
      */
      template<typename EntityType>
      size_t dataHandleSize (const EntityType& e) const
      {
        return TypeTree::reduceOverLeafs(gfs(),DataHandleSize<EntityType>(e),std::plus<size_t>(),size_t(0));
      }

      //! return vector of global indices associated with the given entity
      template<typename EntityType, typename SizeType>
      void dataHandleGlobalIndices (const EntityType& e,
                                    std::vector<SizeType>& global) const
      {
        global.resize(dataHandleSize(e));
        DataHandleGlobalIndicesVisitor<EntityType,std::vector<SizeType> > visitor(e,global);
        TypeTree::applyToTree(gfs(),visitor);
      }

    };






    template<typename GridFunctionSpace, typename Node, typename Mapper>
    class PowerCompositeUpdateAndSetupProvider
    {

#ifndef DOXYGEN

      //! We put the actual method in a base class because we want to use it with different tree iteration patterns
      struct SetupVisitorBase
        : public TypeTree::DefaultVisitor
      {

        template<typename CompositeGFS, typename Child, typename TreePath, typename ChildIndex>
        void afterChild(CompositeGFS& cgfs, const Child& child, TreePath treePath, ChildIndex childIndex)
        {
          cgfs.childGlobalSize[ChildIndex::value] = child.globalSize();
          cgfs.childLocalSize[ChildIndex::value] = child.maxLocalSize();
        }

      };

      //! Visitor for setting up the GFS from pre-initialized children
      struct SetupVisitor
        : public SetupVisitorBase
        , public TypeTree::VisitDirectChildren
      {};

      //! Visitor for updating the complete GFS tree
      struct UpdateVisitor
        : public SetupVisitorBase
        , public TypeTree::VisitTree
      {

        template<typename LeafNode, typename TreePath>
        void leaf(const LeafNode& node, TreePath treePath)
        {
          node.update();
        }

        template<typename LeafNode, typename TreePath>
        void post(const LeafNode& node, TreePath treePath)
        {
          node.calculateSizes();
        }

      };

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
      typedef PowerCompositeGridFunctionSpaceTraits<typename Node::template Child<0>::Type::Traits::GridViewType,
                                                    typename Node::template Child<0>::Type::Traits::BackendType,
                                                    Mapper,
                                                    Node::CHILDREN> Traits;

      //! extract type of container storing Es
      template<typename E>
      struct VectorContainer
      {
        //! \brief define Type as the Type of a container of E's
        typedef typename Traits::BackendType::template VectorContainer<GridFunctionSpace,E> Type;
      private:
        VectorContainer ();
      };

      //! extract type for storing constraints
      template<typename E>
      struct ConstraintsContainer
      {
        //! \brief define Type as the Type of a container of E's
        typedef ConstraintsTransformation<typename Traits::SizeType,E> Type;
      private:
        ConstraintsContainer ();
      };


      //! recalculate sizes
      void update ()
      {
        TypeTree::applyToTree(gfs(),UpdateVisitor());
      }

      //! get dimension of root finite element space
      typename Traits::SizeType globalSize () const
      {
        return offset[Node::CHILDREN];
      }

      //! get dimension of this finite element space
      typename Traits::SizeType size () const
      {
        return offset[Node::CHILDREN];
      }

      //! get max dimension of shape function space
      typename Traits::SizeType maxLocalSize () const
      {
        // this is bullshit !
        return maxlocalsize;
      }

      //! get grid view
      const typename Traits::GridViewType& gridview () const
      {
        return gfs().template child<0>().gridview();
      }

    protected:

      typename Traits::SizeType childGlobalSize[Node::CHILDREN];
      typename Traits::SizeType childLocalSize[Node::CHILDREN];
      typename Traits::SizeType offset[Node::CHILDREN+1];
      typename Traits::SizeType maxlocalsize;

       void setup ()
      {
        TypeTree::applyToTree(gfs(),SetupVisitor());
        gfs().calculateSizes();
      }

    };


    template<typename BaseType, typename T>
    T& checkGridViewType(T& t)
    {
      dune_static_assert((is_same<typename BaseType::Traits::GridViewType,
                          typename T::Traits::GridViewType>::value),
                         "GridViewType must be equal in all components of composite grid function space");
      return t;
    }

    // this partial specialization is required for the non-variadic case
    template<typename BaseType>
    TypeTree::EmptyNode checkGridViewType(TypeTree::EmptyNode e)
    {
      return e;
    }

    /** \brief base class for tuples of grid function spaces
        base class that holds implementation of the methods
        this is the default version with lexicographic ordering
        \tparam Mapper is the ordering parameter. Use e.g.
        \link GridFunctionSpaceLexicographicMapper GridFunctionSpaceLexicographicMapper \endlink
        or \link  GridFunctionSpaceComponentBlockwiseMapper  GridFunctionSpaceComponentBlockwiseMapper \endlink
        or \link  GridFunctionSpaceBlockwiseMapper  GridFunctionSpaceBlockwiseMapper \endlink
        or \link  GridFunctionSpaceDynamicBlockwiseMapper  GridFunctionSpaceDynamicBlockwiseMapper \endlink
        \tparam Ti are all grid function spaces
    */
    template<typename Mapper, DUNE_TYPETREE_COMPOSITENODE_TEMPLATE_CHILDREN>
    class CompositeGridFunctionSpace
    {
    private:
      dune_static_assert(AlwaysFalse<Mapper>::value, "You seem to be using an unsupported Mapper");
    };

    template<DUNE_TYPETREE_COMPOSITENODE_TEMPLATE_CHILDREN_FOR_SPECIALIZATION>
    class CompositeGridFunctionSpace<GridFunctionSpaceLexicographicMapper,
                                     DUNE_TYPETREE_COMPOSITENODE_CHILDTYPES>
      : public DUNE_TYPETREE_COMPOSITENODE_BASETYPE
      , public PowerCompositeDataHandleProvider<CompositeGridFunctionSpace<
                                                  GridFunctionSpaceLexicographicMapper,
                                                  DUNE_TYPETREE_COMPOSITENODE_CHILDTYPES>
                                                >
      , public PowerCompositeUpdateAndSetupProvider<CompositeGridFunctionSpace<
                                                      GridFunctionSpaceLexicographicMapper,
                                                      DUNE_TYPETREE_COMPOSITENODE_CHILDTYPES>,
                                                    DUNE_TYPETREE_COMPOSITENODE_BASETYPE,
                                                    GridFunctionSpaceLexicographicMapper
                                                    >
    {
      typedef DUNE_TYPETREE_COMPOSITENODE_BASETYPE BaseT;

      typedef PowerCompositeUpdateAndSetupProvider<CompositeGridFunctionSpace,BaseT,GridFunctionSpaceLexicographicMapper> ImplementationBase;

      friend class PowerCompositeUpdateAndSetupProvider<CompositeGridFunctionSpace,BaseT,GridFunctionSpaceLexicographicMapper>;

    public:

      typedef typename ImplementationBase::Traits Traits;

      // define local function space parametrized by self
      typedef Dune::PDELab::CompositeLocalFunctionSpaceNode<CompositeGridFunctionSpace> LocalFunctionSpace;

      CompositeGridFunctionSpace(DUNE_TYPETREE_COMPOSITENODE_CONSTRUCTOR_SIGNATURE)
        : BaseT(DUNE_TYPETREE_COMPOSITENODE_CHILDVARIABLES_THROUGH_FUNCTION(checkGridViewType<typename BaseT::template Child<0>::Type>))
      {
        this->setup();
      }

      //! map index from our index set [0,size()-1] to root index set
      typename Traits::SizeType upMap (typename Traits::SizeType i) const
      {
        return i;
      }

      //! map index from child i's index set into our index set
      template<int i>
      typename Traits::SizeType subMap (typename Traits::SizeType j) const
      {
        return this->offset[i]+j;
      }

    private:

      void calculateSizes ()
      {
        Dune::dinfo << "CompositeGridFunctionSpace(lexicographic version):"
                    << std::endl;

        Dune::dinfo << "( ";
        this->offset[0] = 0;
        this->maxlocalsize = 0;
        for (std::size_t i=0; i<BaseT::CHILDREN; i++)
          {
            Dune::dinfo << this->childGlobalSize[i] << " ";
            this->offset[i+1] = this->offset[i]+this->childGlobalSize[i];
            this->maxlocalsize += this->childLocalSize[i];
          }
        Dune::dinfo << ") total size = " << this->offset[BaseT::CHILDREN]
                    << " max local size = " << this->maxlocalsize
                    << std::endl;
      }

    };

#if 0

    // tupel of grid function spaces
    // base class that holds implementation of the methods
    // specialization for blockwise ordering
    // P is the ordering parameter
    // Ti are all grid function spaces
    template<typename T0, typename T1, typename T2, typename T3,
             typename T4, typename T5, typename T6, typename T7, typename T8,
             int s0, int s1, int s2, int s3, int s4, int s5, int s6, int s7, int s8, int s9>
    class CompositeGridFunctionSpaceBase<GridFunctionSpaceComponentBlockwiseMapper<s0,s1,s2,s3,s4,s5,s6,s7,s8,s9>,
                                         T0,T1,T2,T3,T4,T5,T6,T7,T8>
      : public TypeTree::CompositeNode<T0,T1,T2,T3,T4,T5,T6,T7,T8>
    {
      typedef GridFunctionSpaceComponentBlockwiseMapper<s0,s1,s2,s3,s4,s5,s6,s7,s8,s9> BlockwiseMapper;
      typedef CompositeNode<CountingPointerStoragePolicy,T0,T1,T2,T3,T4,T5,T6,T7,T8> BaseT;
    public:
      //! export traits class
      typedef PowerCompositeGridFunctionSpaceTraits<typename Child<0>::Type::Traits::GridViewType,
                                                    typename Child<0>::Type::Traits::BackendType,
                                                    BlockwiseMapper,
                                                    BaseT::CHILDREN>
      Traits;

      //! extract type of container storing Es
      template<typename E>
      struct VectorContainer
      {
        //! \brief define Type as the Type of a container of E's
        typedef typename Traits::BackendType::template VectorContainer<CompositeGridFunctionSpaceBase,E> Type;
      private:
        VectorContainer () {}
      };

      //! extract type for storing constraints
      template<typename E>
      struct ConstraintsContainer
      {
        //! \brief define Type as the Type of a container of E's
        typedef ConstraintsTransformation<typename Traits::SizeType,E> Type;
      private:
        ConstraintsContainer () {}
      };

      // define local function space parametrized by self
      typedef Dune::PDELab::CompositeLocalFunctionSpaceNode<CompositeGridFunctionSpaceBase> LocalFunctionSpace;

      CompositeGridFunctionSpace(DUNE_TYPETREE_COMPOSITENODE_CONSTRUCTOR_SIGNATURE)
        : BaseT(DUNE_TYPETREE_COMPOSITENODE_CHILDVARIABLES_THROUGH_FUNCTION(checkGridViewType<typename Child<0>::Type>))
      {
        setup();
      }

      // get grid view
      const typename Traits::GridViewType& gridview () const
      {
        return this->template getChild<0>().gridview();
      }

      //! get dimension of root finite element space
      typename Traits::SizeType globalSize () const
      {
        return offset[BaseT::CHILDREN];
      }

      //! get dimension of this finite element space
      typename Traits::SizeType size () const
      {
        return offset[BaseT::CHILDREN];
      }

      // get max dimension of shape function space
      typename Traits::SizeType maxLocalSize () const
      {
        // this is bullshit !
        return maxlocalsize;
      }

      //! map index from our index set [0,size()-1] to root index set
      typename Traits::SizeType upMap (typename Traits::SizeType i) const
      {
        return i;
      }

      //! map index from child i's index set into our index set
      template<int i>
      typename Traits::SizeType subMap (typename Traits::SizeType j) const
      {
        // make the block sizes and offsets available in an array
        static const int blockSize[] = { s0, s1, s2, s3, s4, s5, s6, s7, s8, s9 };
        static const int blockOffset[] = { 0, s0, s0+s1, s0+s1+s2, s0+s1+s2+s3, s0+s1+s2+s3+s4,
                                           s0+s1+s2+s3+s4+s5, s0+s1+s2+s3+s4+s5+s6, s0+s1+s2+s3+s4+s5+s6+s7,
                                           s0+s1+s2+s3+s4+s5+s6+s7+s8, s0+s1+s2+s3+s4+s5+s6+s7+s8+s9 };
        return (j%BlockwiseMapper::size[i])
          +(j/BlockwiseMapper::size[i])*BlockwiseMapper::offset[BaseT::CHILDREN]
          +BlockwiseMapper::offset[i];
        return (j%blockSize[i])+(j/blockSize[i])*blockOffset[BaseT::CHILDREN]+blockOffset[i];
      }

      //------------------------------
      // generic data handle interface
      //------------------------------

      //! returns true if data for this codim should be communicated
      bool dataHandleContains (int dim, int codim) const
      {
        return CompositeGridFunctionSpaceBaseVisitChildMetaProgram<CompositeGridFunctionSpaceBase,BaseT::CHILDREN,0>::
          dataHandleContains(*this,dim,codim);
      }

      //! returns true if size per entity of given dim and codim is a constant
      bool dataHandleFixedSize (int dim, int codim) const
      {
        return CompositeGridFunctionSpaceBaseVisitChildMetaProgram<CompositeGridFunctionSpaceBase,BaseT::CHILDREN,0>::
          dataHandleFixedSize(*this,dim,codim);
      }

      /*! how many objects of type DataType have to be sent for a given entity

        Note: Only the sender side needs to know this size.
      */
      template<class EntityType>
      size_t dataHandleSize (const EntityType& e) const
      {
        return CompositeGridFunctionSpaceBaseVisitChildMetaProgram<CompositeGridFunctionSpaceBase,BaseT::CHILDREN,0>::
          dataHandleSize(*this,e);
      }

      //! return vector of global indices associated with the given entity
      template<class EntityType>
      void dataHandleGlobalIndices (const EntityType& e,
                                    std::vector<typename Traits::SizeType>& global) const
      {
        size_t n=dataHandleSize(e);
        global.resize(n);
        CompositeGridFunctionSpaceBaseVisitChildMetaProgram<CompositeGridFunctionSpaceBase,BaseT::CHILDREN,0>::
          dataHandleGlobalIndices(*this,e,global,0,childglobal);
      }

      //------------------------------

      // recalculate sizes
      void update ()
      {
        CompositeGridFunctionSpaceBaseVisitChildMetaProgram<CompositeGridFunctionSpaceBase,BaseT::CHILDREN,0>::
          update(*this);
        setup();
      }

    protected:
      void setup ()
      {
        Dune::dinfo << "CompositeGridFunctionSpace(blockwise version):"
                    << std::endl;

        CompositeGridFunctionSpaceBaseVisitChildMetaProgram<CompositeGridFunctionSpaceBase,BaseT::CHILDREN,0>::
          setup(*this,childGlobalSize,childLocalSize);

        // make the block sizes available in an array
        static const int blockSize[] = { s0, s1, s2, s3, s4, s5, s6, s7, s8, s9 };
        // check for compatible sizes
        for (int i=1; i<BaseT::CHILDREN; i++)
          {
            if (childLocalSize[i]%blockSize[i]!=0)
              DUNE_THROW(Exception,
                         "number of DOFs (" << childLocalSize[i] << ") per component "
                         "must be a multiple of the BlockSize (" << blockSize[i] << ")");
            if (childGlobalSize[i]/blockSize[i]!=childGlobalSize[0]/blockSize[0])
              DUNE_THROW(Exception, "components must be of equal size");
          }

        Dune::dinfo << "( ";
        offset[0] = 0;
        maxlocalsize = 0;
        for (int i=0; i<BaseT::CHILDREN; i++)
          {
            Dune::dinfo << childGlobalSize[i] << " ";
            offset[i+1] = offset[i]+childGlobalSize[i];
            maxlocalsize += childLocalSize[i];
          }
        Dune::dinfo << ") total size = " << offset[BaseT::CHILDREN]
                    << " max local size = " << maxlocalsize
                    << std::endl;
        childglobal.resize(maxlocalsize);
      }

    private:
      typename Traits::SizeType childGlobalSize[BaseT::CHILDREN];
      typename Traits::SizeType childLocalSize[BaseT::CHILDREN];
      typename Traits::SizeType offset[BaseT::CHILDREN+1];
      typename Traits::SizeType maxlocalsize;
      mutable std::vector<typename Traits::SizeType> childglobal;
    };

    template<typename T0, typename T1, typename T2, typename T3,
             typename T4, typename T5, typename T6, typename T7, typename T8>
    class CompositeGridFunctionSpaceBase<GridFunctionSpaceBlockwiseMapper,
                                       T0,T1,T2,T3,T4,T5,T6,T7,T8>
      : public CompositeGridFunctionSpaceBase<GridFunctionSpaceComponentBlockwiseMapper<1>,
                                            T0,T1,T2,T3,T4,T5,T6,T7,T8>
    {
    protected:
      using CompositeGridFunctionSpaceBase<GridFunctionSpaceComponentBlockwiseMapper<1>,
                                           T0,T1,T2,T3,T4,T5,T6,T7,T8>::setup;
    };

    /**
        \brief Tupel of grid function spaces base class that holds
        implementation of the methods specialization for dynamic
        blockwise ordering
    */
    template<DUNE_TYPETREE_COMPOSITENODE_TEMPLATE_CHILDREN>
    class CompositeGridFunctionSpaceBase<GridFunctionSpaceDynamicBlockwiseMapper,DUNE_TYPETREE_COMPOSITENODE_CHILDTYPES>
      : public DUNE_TYPETREE_COMPOSITENODE_BASETYPE
    {
      typedef GridFunctionSpaceDynamicBlockwiseMapper BlockwiseMapper;
      typedef DUNE_TYPETREE_COMPOSITENODE_BASETYPE BaseT;
    public:
      //! export traits class
      typedef PowerCompositeGridFunctionSpaceTraits<typename Child<0>::Type::Traits::GridViewType,
                                                    typename Child<0>::Type::Traits::BackendType,
                                                    BlockwiseMapper,
                                                    BaseT::CHILDREN>
      Traits;

      //! extract type of container storing Es
      template<typename E>
      struct VectorContainer
      {
        //! \brief define Type as the Type of a container of E's
        typedef typename Traits::BackendType::template VectorContainer<CompositeGridFunctionSpaceBase,E> Type;
      private:
        VectorContainer () {}
      };

      //! extract type for storing constraints
      template<typename E>
      struct ConstraintsContainer
      {
        //! \brief define Type as the Type of a container of E's
        typedef ConstraintsTransformation<typename Traits::SizeType,E> Type;
      private:
        ConstraintsContainer () {}
      };

      // define local function space parametrized by self
      typedef Dune::PDELab::CompositeLocalFunctionSpaceNode<CompositeGridFunctionSpaceBase> LocalFunctionSpace;

      CompositeGridFunctionSpace(DUNE_TYPETREE_COMPOSITENODE_CONSTRUCTOR_SIGNATURE)
        : BaseT(DUNE_TYPETREE_COMPOSITENODE_CHILDVARIABLES_THROUGH_FUNCTION(checkGridViewType<typename Child<0>::Type>))
      {
        setup();
      }

      // get grid view
      const typename Traits::GridViewType& gridview () const
      {
        return this->template getChild<0>().gridview();
      }

      //! get dimension of root finite element space
      typename Traits::SizeType globalSize () const
      {
        return offset[BaseT::CHILDREN];
      }

      //! get dimension of this finite element space
      typename Traits::SizeType size () const
      {
        return offset[BaseT::CHILDREN];
      }

      // get max dimension of shape function space
      typename Traits::SizeType maxLocalSize () const
      {
        // this is bullshit !
        return maxlocalsize;
      }

      //! map index from our index set [0,size()-1] to root index set
      typename Traits::SizeType upMap (typename Traits::SizeType i) const
      {
        return i;
      }

      //! map index from child i's index set into our index set
      template<int i>
      typename Traits::SizeType subMap (typename Traits::SizeType j) const
      {
        BlockIndexRangeIterator & it = blockIndexIterators[i];
        while(it->first < j)++it;
        while(it->first > j)--it;
        const typename Traits::SizeType global_index = it->second + j - it->first;
        return global_index;
      }

      //------------------------------
      // generic data handle interface
      //------------------------------

      //! returns true if data for this codim should be communicated
      bool dataHandleContains (int dim, int codim) const
      {
        return CompositeGridFunctionSpaceBaseVisitChildMetaProgram
          <CompositeGridFunctionSpaceBase,BaseT::CHILDREN,0>::
          dataHandleContains(*this,dim,codim);
      }

      //! returns true if size per entity of given dim and codim is a constant
      bool dataHandleFixedSize (int dim, int codim) const
      {
        return CompositeGridFunctionSpaceBaseVisitChildMetaProgram
          <CompositeGridFunctionSpaceBase,BaseT::CHILDREN,0>::
          dataHandleFixedSize(*this,dim,codim);
      }

      /*! how many objects of type DataType have to be sent for a given entity

        Note: Only the sender side needs to know this size.
      */
      template<class EntityType>
      size_t dataHandleSize (const EntityType& e) const
      {
        return CompositeGridFunctionSpaceBaseVisitChildMetaProgram
          <CompositeGridFunctionSpaceBase,BaseT::CHILDREN,0>::
          dataHandleSize(*this,e);
      }

      //! return vector of global indices associated with the given entity
      template<class EntityType>
      void dataHandleGlobalIndices (const EntityType& e,
                                    std::vector<typename Traits::SizeType>& global) const
      {
        size_t n=dataHandleSize(e);
        global.resize(n);
        CompositeGridFunctionSpaceBaseVisitChildMetaProgram<CompositeGridFunctionSpaceBase,BaseT::CHILDREN,0>::
          dataHandleGlobalIndices(*this,e,global,0,childglobal);
      }

      //------------------------------

      // recalculate sizes
      void update ()
      {
        CompositeGridFunctionSpaceBaseVisitChildMetaProgram<CompositeGridFunctionSpaceBase,BaseT::CHILDREN,0>::
          update(*this);
        setup();
      }

    protected:
      void setup ()
      {
        Dune::dinfo << "CompositeGridFunctionSpace(blockwise version):"
                    << std::endl;

        CompositeGridFunctionSpaceBaseVisitChildMetaProgram<CompositeGridFunctionSpaceBase,BaseT::CHILDREN,0>::
          setup(*this,childGlobalSize,childLocalSize);


        typedef typename Traits::GridViewType GridView;
        const GridView & gv = gridview();

        // Initialize offset array for each child

        blockIndices.clear();
        blockIndices.resize(BaseT::CHILDREN);

        std::vector<std::vector<typename Traits::SizeType> > childOffsets;
        childOffsets.resize(BaseT::CHILDREN);
        for(int i=0; i<BaseT::CHILDREN; ++i)
          childOffsets[i].resize(gv.size(0)+1);


        // Iterate grid and determine the offsets for each child
        typedef typename GridView::template Codim<0>::Iterator Iterator;
        Iterator it = gv.template begin<0>();
        const Iterator eit = gv.template end<0>();
        typename Traits::SizeType running_index(0);
        for(; it!=eit; ++it){
          typename Traits::SizeType e_index = gv.indexSet().index(*it);

          // Loop over children (realized by meta-program)
          DynamicBlockwiseMapperImp::GetChildOffsetsMetaProgram<CompositeGridFunctionSpaceBase,BaseT::CHILDREN,0>::
            getChildOffsets(*this,*it,childOffsets);

          for(int i=0; i<BaseT::CHILDREN; ++i){
            if(childOffsets[i][e_index+1] != childOffsets[i][e_index]){
              // Add new block index range element to list
              blockIndices[i].push_back(SizeTypePair(childOffsets[i][e_index],running_index));

              // Update running index
              running_index += childOffsets[i][e_index+1] - childOffsets[i][e_index];
            }
          }
        }

        // Insert a "stop entry" at end of every list
        for(int i=0; i<BaseT::CHILDREN; ++i)
          blockIndices[i].push_back(SizeTypePair(childOffsets[i][gv.size(0)],running_index));


        blockIndexIterators.clear();
        for(int i=0; i<BaseT::CHILDREN; ++i)
          blockIndexIterators.push_back(blockIndices[i].begin());

        // Gather the global information
        Dune::dinfo << "( ";
        offset[0] = 0;
        maxlocalsize = 0;
        for (int i=0; i<BaseT::CHILDREN; i++)
          {
            Dune::dinfo << childGlobalSize[i] << " ";
            offset[i+1] = offset[i]+childGlobalSize[i];
            maxlocalsize += childLocalSize[i];
          }
        Dune::dinfo << ") total size = " << offset[BaseT::CHILDREN]
                    << " max local size = " << maxlocalsize
                    << std::endl;
        childglobal.resize(maxlocalsize);
      }

    private:
      typename Traits::SizeType childGlobalSize[BaseT::CHILDREN];
      typename Traits::SizeType childLocalSize[BaseT::CHILDREN];
      typename Traits::SizeType offset[BaseT::CHILDREN+1];
      typename Traits::SizeType maxlocalsize;
      mutable std::vector<typename Traits::SizeType> childglobal;
      typedef std::pair<typename Traits::SizeType,typename Traits::SizeType> SizeTypePair;
      typedef std::list<SizeTypePair> BlockIndexRangeList;
      std::vector<BlockIndexRangeList> blockIndices;
      typedef typename BlockIndexRangeList::const_iterator BlockIndexRangeIterator;
      mutable std::vector<BlockIndexRangeIterator> blockIndexIterators;
    };

#endif

    //! \addtogroup GridFunctionSpace
    //! \{

    /** \brief grid function space composed of other grid function spaces

        Composes a tuple of arbitray grid function spaces into a grid function space.
        The ordering of the resulting unknowns can be done lexicographically or block-wise.
        \tparam Mapper is the ordering parameter. Use e.g.
        \link GridFunctionSpaceLexicographicMapper GridFunctionSpaceLexicographicMapper \endlink
        or \link  GridFunctionSpaceComponentBlockwiseMapper  GridFunctionSpaceComponentBlockwiseMapper \endlink
        or \link  GridFunctionSpaceBlockwiseMapper  GridFunctionSpaceBlockwiseMapper \endlink
        or \link  GridFunctionSpaceDynamicBlockwiseMapper  GridFunctionSpaceDynamicBlockwiseMapper \endlink
    */

  }

}
//! \}
#endif
