// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_COMPOSITEGRIDFUNCTIONSPACE_HH
#define DUNE_PDELAB_COMPOSITEGRIDFUNCTIONSPACE_HH

#include "gridfunctionspace.hh"

namespace Dune {
  namespace PDELab {

    //=======================================
    // composite grid function space
    //=======================================

    //! \addtogroup GridFunctionSpace grid function space
    //! \ingroup PDELab
    //! \{

    template<typename T, int n, int i>
    struct CompositeGridFunctionSpaceBaseVisitChildMetaProgram // visit child of inner node
    {
      template<typename Int>
      static void setup (T& t, Int childGlobalSize[], Int childLocalSize[])
      {
        childGlobalSize[i] = t.template getChild<i>().globalSize();
        childLocalSize[i] = t.template getChild<i>().maxLocalSize();
        CompositeGridFunctionSpaceBaseVisitChildMetaProgram<T,n,i+1>::
          setup(t,childGlobalSize,childLocalSize);
      }
      static void update (T& t)
      {
        t.template getChild<i>().update();
        CompositeGridFunctionSpaceBaseVisitChildMetaProgram<T,n,i+1>::update(t);
      }
      static bool dataHandleContains (const T& t, int dim, int codim)
      {
        return t.template getChild<i>().dataHandleContains(dim,codim) || 
          CompositeGridFunctionSpaceBaseVisitChildMetaProgram<T,n,i+1>::dataHandleContains(t,dim,codim);
      }
      static bool dataHandleFixedSize (const T& t, int dim, int codim)
      {
        return t.template getChild<i>().dataHandleFixedSize(dim,codim) && 
          CompositeGridFunctionSpaceBaseVisitChildMetaProgram<T,n,i+1>::dataHandleFixedSize(t,dim,codim);
      }
      template<class EntityType>
      static size_t dataHandleSize (const T& t, const EntityType& e)
      {
        return t.template getChild<i>().dataHandleSize(e) + 
          CompositeGridFunctionSpaceBaseVisitChildMetaProgram<T,n,i+1>::dataHandleSize(t,e);
      }
      template<class EntityType, class C>
      static void dataHandleGlobalIndices (const T& t, const EntityType& e, C& global, size_t ng, C& childglobal)
      {
        size_t nc=t.template getChild<i>().dataHandleSize(e);
        childglobal.resize(nc);
        t.template getChild<i>().dataHandleGlobalIndices(e,childglobal);
        for (size_t j=0; j<childglobal.size(); j++)
          global[ng+j] = t.template subMap<i>(childglobal[j]);
        CompositeGridFunctionSpaceBaseVisitChildMetaProgram<T,n,i+1>::dataHandleGlobalIndices(t,e,global,ng+nc,childglobal);
      }
    };

    template<typename T, int n>
    struct CompositeGridFunctionSpaceBaseVisitChildMetaProgram<T,n,n> // end of child recursion
    {
      template<typename Int>
      static void setup (T& t, Int childGlobalSize[], Int childLocalSize[])
      {
      }
      static void update (T& t)
      {
      }
      static bool dataHandleContains (const T& t, int dim, int codim)
      {
        return false;
      }
      static bool dataHandleFixedSize (const T& t, int dim, int codim)
      {
        return true;
      }
      template<class EntityType>
      static size_t dataHandleSize (const T& t, const EntityType& e)
      {
        return 0;
      }
      template<class EntityType, class C>
      static void dataHandleGlobalIndices (const T& t, const EntityType& e, C& global, size_t n_, C& childglobal)
      {
      }
    };

    namespace
    {
      // TMP to count the number of non empty entries
      template< typename T0, typename T1, typename T2, typename T3,
                typename T4, typename T5, typename T6, typename T7, typename T8>
      struct NonEmptyChilds
      {
        enum{ value = 9 };
      };
      
      template< typename T0, typename T1, typename T2, typename T3,
                typename T4, typename T5, typename T6, typename T7>
      struct NonEmptyChilds<T0,T1,T2,T3,T4,T5,T6,T7,EmptyChild>
      {
        enum{ value = 8 };
      };

      template< typename T0, typename T1, typename T2, typename T3,
                typename T4, typename T5, typename T6>
      struct NonEmptyChilds<T0,T1,T2,T3,T4,T5,T6,EmptyChild,EmptyChild>
      {
        enum{ value = 7 };
      };
      
      template< typename T0, typename T1, typename T2, typename T3,
                typename T4, typename T5>
      struct NonEmptyChilds<T0,T1,T2,T3,T4,T5,EmptyChild,EmptyChild,EmptyChild>
      {
        enum{ value = 6 };
      };

      template< typename T0, typename T1, typename T2, typename T3,
                typename T4>
      struct NonEmptyChilds<T0,T1,T2,T3,T4,EmptyChild,EmptyChild,EmptyChild,EmptyChild>
      {
        enum{ value = 5 };
      };

      template< typename T0, typename T1, typename T2, typename T3>
      struct NonEmptyChilds<T0,T1,T2,T3,EmptyChild,EmptyChild,EmptyChild,EmptyChild,EmptyChild>
      {
        enum{ value = 4 };
      };
      
      template< typename T0, typename T1, typename T2>
      struct NonEmptyChilds<T0,T1,T2,EmptyChild,EmptyChild,EmptyChild,EmptyChild,EmptyChild,EmptyChild>
      {
        enum{ value = 3 };
      };
      
      template< typename T0, typename T1>
      struct NonEmptyChilds<T0,T1,EmptyChild,EmptyChild,EmptyChild,EmptyChild,EmptyChild,EmptyChild,EmptyChild>
      {
        enum{ value = 2 };
      };
           
      template< typename T0>
      struct NonEmptyChilds<T0,EmptyChild,EmptyChild,EmptyChild,EmptyChild,EmptyChild,EmptyChild,EmptyChild,EmptyChild>
      {
        enum{ value = 1 };
      };
     
      template<>
      struct NonEmptyChilds<EmptyChild,EmptyChild,EmptyChild,EmptyChild,EmptyChild,EmptyChild,EmptyChild,EmptyChild,EmptyChild>
      {
        enum{ value = 1 };
      };
    }
    

    template<typename Mapper, typename T0, typename T1, typename T2, typename T3,
             typename T4, typename T5, typename T6, typename T7, typename T8>
    class CompositeGridFunctionSpace;

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
    template<typename Mapper, typename T0, typename T1, typename T2=EmptyChild, typename T3=EmptyChild,
             typename T4=EmptyChild, typename T5=EmptyChild, typename T6=EmptyChild,
             typename T7=EmptyChild, typename T8=EmptyChild>
    class CompositeGridFunctionSpaceBase;

    template<typename T0, typename T1, typename T2, typename T3,
             typename T4, typename T5, typename T6, typename T7, typename T8>
    class CompositeGridFunctionSpaceBase<GridFunctionSpaceLexicographicMapper,
                                         T0,T1,T2,T3,T4,T5,T6,T7,T8>
      : public CompositeNode<CountingPointerStoragePolicy,T0,T1,T2,T3,T4,T5,T6,T7,T8>,
        public Countable
    {
      typedef CompositeNode<CountingPointerStoragePolicy,T0,T1,T2,T3,T4,T5,T6,T7,T8> BaseT;
    public:
      //! export traits class
      typedef PowerCompositeGridFunctionSpaceTraits<typename T0::Traits::GridViewType, 
                                                    typename T0::Traits::BackendType,
                                                    GridFunctionSpaceLexicographicMapper,
                                                    NonEmptyChilds<T0,T1,T2,T3,T4,T5,
                                                                   T6,T7,T8>::value>
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


      // it is part of the trick to have a constructor without arguments
      // setting of the children is then done by the constructors
      // of the specialized derived classes
      CompositeGridFunctionSpaceBase ()
      {
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
        return offset[i]+j;
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
        Dune::dinfo << "CompositeGridFunctionSpace(lexicographic version):"
                    << std::endl;

        CompositeGridFunctionSpaceBaseVisitChildMetaProgram<CompositeGridFunctionSpaceBase,BaseT::CHILDREN,0>::
          setup(*this,childGlobalSize,childLocalSize);

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
      : public CompositeNode<CountingPointerStoragePolicy,T0,T1,T2,T3,T4,T5,T6,T7,T8>,
        public Countable
    {
      typedef GridFunctionSpaceComponentBlockwiseMapper<s0,s1,s2,s3,s4,s5,s6,s7,s8,s9> BlockwiseMapper;
      typedef CompositeNode<CountingPointerStoragePolicy,T0,T1,T2,T3,T4,T5,T6,T7,T8> BaseT;
    public:
      //! export traits class
      typedef PowerCompositeGridFunctionSpaceTraits<typename T0::Traits::GridViewType, 
                                                    typename T0::Traits::BackendType,
                                                    BlockwiseMapper,
                                                    NonEmptyChilds<T0,T1,T2,T3,T4,T5,
                                                                   T6,T7,T8>::value>
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

      // it is part of the trick to have a constructor without arguments
      // setting of the children is then done by the constructors
      // of the specialized derived classes
      CompositeGridFunctionSpaceBase ()
      {
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
    template<typename T0, typename T1, typename T2, typename T3,
             typename T4, typename T5, typename T6, typename T7, typename T8>
    class CompositeGridFunctionSpaceBase<GridFunctionSpaceDynamicBlockwiseMapper,
                                         T0,T1,T2,T3,T4,T5,T6,T7,T8>
      : public CompositeNode<CountingPointerStoragePolicy,T0,T1,T2,T3,T4,T5,T6,T7,T8>,
        public Countable
    {
      typedef GridFunctionSpaceDynamicBlockwiseMapper BlockwiseMapper;
      typedef CompositeNode<CountingPointerStoragePolicy,T0,T1,T2,T3,T4,T5,T6,T7,T8> BaseT;
    public:
      //! export traits class
      typedef PowerCompositeGridFunctionSpaceTraits<typename T0::Traits::GridViewType, 
                                                    typename T0::Traits::BackendType,
                                                    BlockwiseMapper,
                                                    NonEmptyChilds<T0,T1,T2,T3,T4,T5,
                                                                   T6,T7,T8>::value>
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

      // it is part of the trick to have a constructor without arguments
      // setting of the children is then done by the constructors
      // of the specialized derived classes
      CompositeGridFunctionSpaceBase ()
      {}

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
    template<typename Mapper, typename T0, typename T1, typename T2=EmptyChild, typename T3=EmptyChild,
             typename T4=EmptyChild, typename T5=EmptyChild, typename T6=EmptyChild,
             typename T7=EmptyChild, typename T8=EmptyChild>
    class CompositeGridFunctionSpace : 
      public CompositeGridFunctionSpaceBase<Mapper,T0,T1,T2,T3,T4,T5,T6,T7,T8>
    {
      typedef CompositeGridFunctionSpaceBase<Mapper,T0,T1,T2,T3,T4,T5,T6,T7,T8> BaseT;

    public:
      //! export traits class
      typedef typename BaseT::Traits Traits;

      //! Construct a PowerGridFunction with k clones of the function t
      CompositeGridFunctionSpace (T0& t0, T1& t1, T2& t2, T3& t3, T4& t4, T5& t5, T6& t6, T7& t7, T8& t8)
      {
        dune_static_assert((is_same<typename T0::Traits::GridViewType,
                            typename T1::Traits::GridViewType>::value),  
                           "GridViewType must be equal in all components of composite grid function space");
        dune_static_assert((is_same<typename T0::Traits::BackendType,
                            typename T1::Traits::BackendType>::value),  
                           "BackendType must be equal in all components of composite grid function space");

        dune_static_assert((is_same<typename T0::Traits::GridViewType,
                            typename T2::Traits::GridViewType>::value),  
                           "GridViewType must be equal in all components of composite grid function space");
        dune_static_assert((is_same<typename T0::Traits::BackendType,
                            typename T2::Traits::BackendType>::value),  
                           "BackendType must be equal in all components of composite grid function space");

        dune_static_assert((is_same<typename T0::Traits::GridViewType,
                            typename T3::Traits::GridViewType>::value),  
                           "GridViewType must be equal in all components of composite grid function space");
        dune_static_assert((is_same<typename T0::Traits::BackendType,
                            typename T3::Traits::BackendType>::value),  
                           "BackendType must be equal in all components of composite grid function space");

        dune_static_assert((is_same<typename T0::Traits::GridViewType,
                            typename T4::Traits::GridViewType>::value),  
                           "GridViewType must be equal in all components of composite grid function space");
        dune_static_assert((is_same<typename T0::Traits::BackendType,
                            typename T4::Traits::BackendType>::value),  
                           "BackendType must be equal in all components of composite grid function space");

        dune_static_assert((is_same<typename T0::Traits::GridViewType,
                            typename T5::Traits::GridViewType>::value),  
                           "GridViewType must be equal in all components of composite grid function space");
        dune_static_assert((is_same<typename T0::Traits::BackendType,
                            typename T5::Traits::BackendType>::value),  
                           "BackendType must be equal in all components of composite grid function space");

        dune_static_assert((is_same<typename T0::Traits::GridViewType,
                            typename T6::Traits::GridViewType>::value),  
                           "GridViewType must be equal in all components of composite grid function space");
        dune_static_assert((is_same<typename T0::Traits::BackendType,
                            typename T6::Traits::BackendType>::value),  
                           "BackendType must be equal in all components of composite grid function space");

        dune_static_assert((is_same<typename T0::Traits::GridViewType,
                            typename T7::Traits::GridViewType>::value),  
                           "GridViewType must be equal in all components of composite grid function space");
        dune_static_assert((is_same<typename T0::Traits::BackendType,
                            typename T7::Traits::BackendType>::value),  
                           "BackendType must be equal in all components of composite grid function space");

        dune_static_assert((is_same<typename T0::Traits::GridViewType,
                            typename T8::Traits::GridViewType>::value),  
                           "GridViewType must be equal in all components of composite grid function space");
        dune_static_assert((is_same<typename T0::Traits::BackendType,
                            typename T8::Traits::BackendType>::value),  
                           "BackendType must be equal in all components of composite grid function space");


        this->template setChild<0>(t0);
        this->template setChild<1>(t1);
        this->template setChild<2>(t2);
        this->template setChild<3>(t3);
        this->template setChild<4>(t4);
        this->template setChild<5>(t5);
        this->template setChild<6>(t6);
        this->template setChild<7>(t7);
        this->template setChild<8>(t8);

        BaseT::setup();
      }
    };

    //! \}

    template<typename P, typename T0, typename T1>
    class CompositeGridFunctionSpace<P,T0,T1> 
      : public CompositeGridFunctionSpaceBase<P,T0,T1,EmptyChild,EmptyChild,EmptyChild,
                                              EmptyChild,EmptyChild,EmptyChild,EmptyChild>
    {
      typedef CompositeGridFunctionSpaceBase<P,T0,T1,EmptyChild,EmptyChild,EmptyChild,
                                             EmptyChild,EmptyChild,EmptyChild,EmptyChild> BaseT;
    public:
      //! export traits class
      typedef typename BaseT::Traits Traits;

      //! Construct a PowerGridFunction with k clones of the function t
      CompositeGridFunctionSpace (T0& t0, T1& t1)
      {
        dune_static_assert((is_same<typename T0::Traits::GridViewType,
                            typename T1::Traits::GridViewType>::value),  
                           "GridViewType must be equal in all components of composite grid function space");
        dune_static_assert((is_same<typename T0::Traits::BackendType,
                            typename T1::Traits::BackendType>::value),  
                           "BackendType must be equal in all components of composite grid function space");
        this->template setChild<0>(t0);
        this->template setChild<1>(t1);
        BaseT::setup();
      }
    };

    template<typename P, typename T0, typename T1, typename T2>
    class CompositeGridFunctionSpace<P,T0,T1,T2> 
      : public CompositeGridFunctionSpaceBase<P,T0,T1,T2,EmptyChild,EmptyChild,
                                              EmptyChild,EmptyChild,EmptyChild,EmptyChild>
    {
      typedef CompositeGridFunctionSpaceBase<P,T0,T1,T2,EmptyChild,EmptyChild,
                                             EmptyChild,EmptyChild,EmptyChild,EmptyChild> BaseT;
    public:
      //! export traits class
      typedef typename BaseT::Traits Traits;

      //! Construct a PowerGridFunction with k clones of the function t
      CompositeGridFunctionSpace (T0& t0, T1& t1, T2& t2)
      {
        dune_static_assert((is_same<typename T0::Traits::GridViewType,
                            typename T1::Traits::GridViewType>::value),  
                           "GridViewType must be equal in all components of composite grid function space");
        dune_static_assert((is_same<typename T0::Traits::BackendType,
                            typename T1::Traits::BackendType>::value),  
                           "BackendType must be equal in all components of composite grid function space");
        dune_static_assert((is_same<typename T0::Traits::GridViewType,
                            typename T2::Traits::GridViewType>::value),  
                           "GridViewType must be equal in all components of composite grid function space");
        dune_static_assert((is_same<typename T0::Traits::BackendType,
                            typename T2::Traits::BackendType>::value),  
                           "BackendType must be equal in all components of composite grid function space");
        this->template setChild<0>(t0);
        this->template setChild<1>(t1);
        this->template setChild<2>(t2);
        BaseT::setup();
      }
    };

    template<typename P, typename T0, typename T1, typename T2, typename T3>
    class CompositeGridFunctionSpace<P,T0,T1,T2,T3> 
      : public CompositeGridFunctionSpaceBase<P,T0,T1,T2,T3,EmptyChild,
                                              EmptyChild,EmptyChild,EmptyChild,EmptyChild>
    {
      typedef CompositeGridFunctionSpaceBase<P,T0,T1,T2,T3,EmptyChild,
                                             EmptyChild,EmptyChild,EmptyChild,EmptyChild> BaseT;
    public:
      //! export traits class
      typedef typename BaseT::Traits Traits;

      //! Construct a PowerGridFunction with k clones of the function t
      CompositeGridFunctionSpace (T0& t0, T1& t1, T2& t2, T3& t3)
      {
        dune_static_assert((is_same<typename T0::Traits::GridViewType,
                            typename T1::Traits::GridViewType>::value),  
                           "GridViewType must be equal in all components of composite grid function space");
        dune_static_assert((is_same<typename T0::Traits::BackendType,
                            typename T1::Traits::BackendType>::value),  
                           "BackendType must be equal in all components of composite grid function space");

        dune_static_assert((is_same<typename T0::Traits::GridViewType,
                            typename T2::Traits::GridViewType>::value),  
                           "GridViewType must be equal in all components of composite grid function space");
        dune_static_assert((is_same<typename T0::Traits::BackendType,
                            typename T2::Traits::BackendType>::value),  
                           "BackendType must be equal in all components of composite grid function space");

        dune_static_assert((is_same<typename T0::Traits::GridViewType,
                            typename T3::Traits::GridViewType>::value),  
                           "GridViewType must be equal in all components of composite grid function space");
        dune_static_assert((is_same<typename T0::Traits::BackendType,
                            typename T3::Traits::BackendType>::value),  
                           "BackendType must be equal in all components of composite grid function space");

        this->template setChild<0>(t0);
        this->template setChild<1>(t1);
        this->template setChild<2>(t2);
        this->template setChild<3>(t3);

        BaseT::setup();
      }
    };

    template<typename P, typename T0, typename T1, typename T2, typename T3, typename T4>
    class CompositeGridFunctionSpace<P,T0,T1,T2,T3,T4> 
      : public CompositeGridFunctionSpaceBase<P,T0,T1,T2,T3,T4,
                                              EmptyChild,EmptyChild,EmptyChild,EmptyChild>
    {
      typedef CompositeGridFunctionSpaceBase<P,T0,T1,T2,T3,T4,
                                             EmptyChild,EmptyChild,EmptyChild,EmptyChild> BaseT;
    public:
      //! export traits class
      typedef typename BaseT::Traits Traits;

      //! Construct a PowerGridFunction with k clones of the function t
      CompositeGridFunctionSpace (T0& t0, T1& t1, T2& t2, T3& t3, T4& t4)
      {
        dune_static_assert((is_same<typename T0::Traits::GridViewType,
                            typename T1::Traits::GridViewType>::value),  
                           "GridViewType must be equal in all components of composite grid function space");
        dune_static_assert((is_same<typename T0::Traits::BackendType,
                            typename T1::Traits::BackendType>::value),  
                           "BackendType must be equal in all components of composite grid function space");

        dune_static_assert((is_same<typename T0::Traits::GridViewType,
                            typename T2::Traits::GridViewType>::value),  
                           "GridViewType must be equal in all components of composite grid function space");
        dune_static_assert((is_same<typename T0::Traits::BackendType,
                            typename T2::Traits::BackendType>::value),  
                           "BackendType must be equal in all components of composite grid function space");

        dune_static_assert((is_same<typename T0::Traits::GridViewType,
                            typename T3::Traits::GridViewType>::value),  
                           "GridViewType must be equal in all components of composite grid function space");
        dune_static_assert((is_same<typename T0::Traits::BackendType,
                            typename T3::Traits::BackendType>::value),  
                           "BackendType must be equal in all components of composite grid function space");

        dune_static_assert((is_same<typename T0::Traits::GridViewType,
                            typename T4::Traits::GridViewType>::value),  
                           "GridViewType must be equal in all components of composite grid function space");
        dune_static_assert((is_same<typename T0::Traits::BackendType,
                            typename T4::Traits::BackendType>::value),  
                           "BackendType must be equal in all components of composite grid function space");

        this->template setChild<0>(t0);
        this->template setChild<1>(t1);
        this->template setChild<2>(t2);
        this->template setChild<3>(t3);
        this->template setChild<4>(t4);

        BaseT::setup();
      }
    };

    template<typename P, typename T0, typename T1, typename T2, typename T3, typename T4,
             typename T5>
    class CompositeGridFunctionSpace<P,T0,T1,T2,T3,T4,T5> 
      : public CompositeGridFunctionSpaceBase<P,T0,T1,T2,T3,T4,
                                              T5,EmptyChild,EmptyChild,EmptyChild>
    {
      typedef CompositeGridFunctionSpaceBase<P,T0,T1,T2,T3,T4,
                                             T5,EmptyChild,EmptyChild,EmptyChild> BaseT;
    public:
      //! export traits class
      typedef typename BaseT::Traits Traits;

      //! Construct a PowerGridFunction with k clones of the function t
      CompositeGridFunctionSpace (T0& t0, T1& t1, T2& t2, T3& t3, T4& t4, T5& t5)
      {
        dune_static_assert((is_same<typename T0::Traits::GridViewType,
                            typename T1::Traits::GridViewType>::value),  
                           "GridViewType must be equal in all components of composite grid function space");
        dune_static_assert((is_same<typename T0::Traits::BackendType,
                            typename T1::Traits::BackendType>::value),  
                           "BackendType must be equal in all components of composite grid function space");

        dune_static_assert((is_same<typename T0::Traits::GridViewType,
                            typename T2::Traits::GridViewType>::value),  
                           "GridViewType must be equal in all components of composite grid function space");
        dune_static_assert((is_same<typename T0::Traits::BackendType,
                            typename T2::Traits::BackendType>::value),  
                           "BackendType must be equal in all components of composite grid function space");

        dune_static_assert((is_same<typename T0::Traits::GridViewType,
                            typename T3::Traits::GridViewType>::value),  
                           "GridViewType must be equal in all components of composite grid function space");
        dune_static_assert((is_same<typename T0::Traits::BackendType,
                            typename T3::Traits::BackendType>::value),  
                           "BackendType must be equal in all components of composite grid function space");

        dune_static_assert((is_same<typename T0::Traits::GridViewType,
                            typename T4::Traits::GridViewType>::value),  
                           "GridViewType must be equal in all components of composite grid function space");
        dune_static_assert((is_same<typename T0::Traits::BackendType,
                            typename T4::Traits::BackendType>::value),  
                           "BackendType must be equal in all components of composite grid function space");

        dune_static_assert((is_same<typename T0::Traits::GridViewType,
                            typename T5::Traits::GridViewType>::value),  
                           "GridViewType must be equal in all components of composite grid function space");
        dune_static_assert((is_same<typename T0::Traits::BackendType,
                            typename T5::Traits::BackendType>::value),  
                           "BackendType must be equal in all components of composite grid function space");

        this->template setChild<0>(t0);
        this->template setChild<1>(t1);
        this->template setChild<2>(t2);
        this->template setChild<3>(t3);
        this->template setChild<4>(t4);
        this->template setChild<5>(t5);

        BaseT::setup();
      }
    };

    template<typename P, typename T0, typename T1, typename T2, typename T3, typename T4,
             typename T5, typename T6>
    class CompositeGridFunctionSpace<P,T0,T1,T2,T3,T4,T5,T6> 
      : public CompositeGridFunctionSpaceBase<P,T0,T1,T2,T3,T4,
                                              T5,T6,EmptyChild,EmptyChild>
    {
      typedef CompositeGridFunctionSpaceBase<P,T0,T1,T2,T3,T4,
                                             T5,T6,EmptyChild,EmptyChild> BaseT;
    public:
      //! export traits class
      typedef typename BaseT::Traits Traits;

      //! Construct a PowerGridFunction with k clones of the function t
      CompositeGridFunctionSpace (T0& t0, T1& t1, T2& t2, T3& t3, T4& t4, T5& t5, T6& t6)
      {
        dune_static_assert((is_same<typename T0::Traits::GridViewType,
                            typename T1::Traits::GridViewType>::value),  
                           "GridViewType must be equal in all components of composite grid function space");
        dune_static_assert((is_same<typename T0::Traits::BackendType,
                            typename T1::Traits::BackendType>::value),  
                           "BackendType must be equal in all components of composite grid function space");

        dune_static_assert((is_same<typename T0::Traits::GridViewType,
                            typename T2::Traits::GridViewType>::value),  
                           "GridViewType must be equal in all components of composite grid function space");
        dune_static_assert((is_same<typename T0::Traits::BackendType,
                            typename T2::Traits::BackendType>::value),  
                           "BackendType must be equal in all components of composite grid function space");

        dune_static_assert((is_same<typename T0::Traits::GridViewType,
                            typename T3::Traits::GridViewType>::value),  
                           "GridViewType must be equal in all components of composite grid function space");
        dune_static_assert((is_same<typename T0::Traits::BackendType,
                            typename T3::Traits::BackendType>::value),  
                           "BackendType must be equal in all components of composite grid function space");

        dune_static_assert((is_same<typename T0::Traits::GridViewType,
                            typename T4::Traits::GridViewType>::value),  
                           "GridViewType must be equal in all components of composite grid function space");
        dune_static_assert((is_same<typename T0::Traits::BackendType,
                            typename T4::Traits::BackendType>::value),  
                           "BackendType must be equal in all components of composite grid function space");

        dune_static_assert((is_same<typename T0::Traits::GridViewType,
                            typename T5::Traits::GridViewType>::value),  
                           "GridViewType must be equal in all components of composite grid function space");
        dune_static_assert((is_same<typename T0::Traits::BackendType,
                            typename T5::Traits::BackendType>::value),  
                           "BackendType must be equal in all components of composite grid function space");

        dune_static_assert((is_same<typename T0::Traits::GridViewType,
                            typename T6::Traits::GridViewType>::value),  
                           "GridViewType must be equal in all components of composite grid function space");
        dune_static_assert((is_same<typename T0::Traits::BackendType,
                            typename T6::Traits::BackendType>::value),  
                           "BackendType must be equal in all components of composite grid function space");

        this->template setChild<0>(t0);
        this->template setChild<1>(t1);
        this->template setChild<2>(t2);
        this->template setChild<3>(t3);
        this->template setChild<4>(t4);
        this->template setChild<5>(t5);
        this->template setChild<6>(t6);

        BaseT::setup();
      }
    };

    template<typename P, typename T0, typename T1, typename T2, typename T3, typename T4,
             typename T5, typename T6, typename T7>
    class CompositeGridFunctionSpace<P,T0,T1,T2,T3,T4,T5,T6,T7> 
      : public CompositeGridFunctionSpaceBase<P,T0,T1,T2,T3,T4,
                                              T5,T6,T7,EmptyChild>
    {
      typedef CompositeGridFunctionSpaceBase<P,T0,T1,T2,T3,T4,
                                             T5,T6,T7,EmptyChild> BaseT;
    public:
      //! export traits class
      typedef typename BaseT::Traits Traits;

      //! Construct a PowerGridFunction with k clones of the function t
      CompositeGridFunctionSpace (T0& t0, T1& t1, T2& t2, T3& t3, T4& t4, T5& t5, T6& t6, T7& t7)
      {
        dune_static_assert((is_same<typename T0::Traits::GridViewType,
                            typename T1::Traits::GridViewType>::value),  
                           "GridViewType must be equal in all components of composite grid function space");
        dune_static_assert((is_same<typename T0::Traits::BackendType,
                            typename T1::Traits::BackendType>::value),  
                           "BackendType must be equal in all components of composite grid function space");

        dune_static_assert((is_same<typename T0::Traits::GridViewType,
                            typename T2::Traits::GridViewType>::value),  
                           "GridViewType must be equal in all components of composite grid function space");
        dune_static_assert((is_same<typename T0::Traits::BackendType,
                            typename T2::Traits::BackendType>::value),  
                           "BackendType must be equal in all components of composite grid function space");

        dune_static_assert((is_same<typename T0::Traits::GridViewType,
                            typename T3::Traits::GridViewType>::value),  
                           "GridViewType must be equal in all components of composite grid function space");
        dune_static_assert((is_same<typename T0::Traits::BackendType,
                            typename T3::Traits::BackendType>::value),  
                           "BackendType must be equal in all components of composite grid function space");

        dune_static_assert((is_same<typename T0::Traits::GridViewType,
                            typename T4::Traits::GridViewType>::value),  
                           "GridViewType must be equal in all components of composite grid function space");
        dune_static_assert((is_same<typename T0::Traits::BackendType,
                            typename T4::Traits::BackendType>::value),  
                           "BackendType must be equal in all components of composite grid function space");

        dune_static_assert((is_same<typename T0::Traits::GridViewType,
                            typename T5::Traits::GridViewType>::value),  
                           "GridViewType must be equal in all components of composite grid function space");
        dune_static_assert((is_same<typename T0::Traits::BackendType,
                            typename T5::Traits::BackendType>::value),  
                           "BackendType must be equal in all components of composite grid function space");

        dune_static_assert((is_same<typename T0::Traits::GridViewType,
                            typename T6::Traits::GridViewType>::value),  
                           "GridViewType must be equal in all components of composite grid function space");
        dune_static_assert((is_same<typename T0::Traits::BackendType,
                            typename T6::Traits::BackendType>::value),  
                           "BackendType must be equal in all components of composite grid function space");

        dune_static_assert((is_same<typename T0::Traits::GridViewType,
                            typename T7::Traits::GridViewType>::value),  
                           "GridViewType must be equal in all components of composite grid function space");
        dune_static_assert((is_same<typename T0::Traits::BackendType,
                            typename T7::Traits::BackendType>::value),  
                           "BackendType must be equal in all components of composite grid function space");

        this->template setChild<0>(t0);
        this->template setChild<1>(t1);
        this->template setChild<2>(t2);
        this->template setChild<3>(t3);
        this->template setChild<4>(t4);
        this->template setChild<5>(t5);
        this->template setChild<6>(t6);
        this->template setChild<7>(t7);

        BaseT::setup();
      }
    };

  }

}
//! \}
#endif
