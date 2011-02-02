// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_COMPOSITEGRIDFUNCTIONSPACE_HH
#define DUNE_PDELAB_COMPOSITEGRIDFUNCTIONSPACE_HH

#include "powercompositegridfunctionspacemixins.hh"

namespace Dune {
  namespace PDELab {

    //=======================================
    // composite grid function space
    //=======================================

    //! \addtogroup GridFunctionSpace grid function space
    //! \ingroup PDELab
    //! \{

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
      , public PowerCompositeGridFunctionSpaceBase<CompositeGridFunctionSpace<
                                                     GridFunctionSpaceLexicographicMapper,
                                                     DUNE_TYPETREE_COMPOSITENODE_CHILDTYPES>,
                                                   DUNE_TYPETREE_COMPOSITENODE_BASETYPE,
                                                   GridFunctionSpaceLexicographicMapper
                                                   >
    {
      typedef DUNE_TYPETREE_COMPOSITENODE_BASETYPE BaseT;

      typedef PowerCompositeGridFunctionSpaceBase<CompositeGridFunctionSpace,BaseT,GridFunctionSpaceLexicographicMapper> ImplementationBase;

      friend class PowerCompositeGridFunctionSpaceBase<CompositeGridFunctionSpace,BaseT,GridFunctionSpaceLexicographicMapper>;

    public:

      typedef CompositeGridFunctionSpaceTag ImplementationTag;

      typedef typename ImplementationBase::Traits Traits;

      // define local function space parametrized by self
      typedef typename Dune::PDELab::TypeTree::TransformTree<CompositeGridFunctionSpace,gfs_to_lfs>::Type LocalFunctionSpace;

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
        return subMap(i,j);
      }

      typename Traits::SizeType subMap (typename Traits::SizeType i, typename Traits::SizeType j) const
      {
        return this->offset[i]+j;
      }

    private:

      using ImplementationBase::childLocalSize;
      using ImplementationBase::childGlobalSize;
      using ImplementationBase::maxlocalsize;
      using ImplementationBase::offset;

      void calculateSizes ()
      {
        Dune::dinfo << "CompositeGridFunctionSpace(lexicographic version):"
                    << std::endl;

        Dune::dinfo << "( ";
        offset[0] = 0;
        maxlocalsize = 0;
        for (std::size_t i=0; i<BaseT::CHILDREN; i++)
          {
            Dune::dinfo << childGlobalSize[i] << " ";
            offset[i+1] = offset[i]+childGlobalSize[i];
            maxlocalsize += childLocalSize[i];
          }
        Dune::dinfo << ") total size = " << offset[BaseT::CHILDREN]
                    << " max local size = " << maxlocalsize
                    << std::endl;
      }

    };

    // tupel of grid function spaces
    // base class that holds implementation of the methods
    // specialization for blockwise ordering
    // P is the ordering parameter
    // Ti are all grid function spaces
    template<DUNE_TYPETREE_COMPOSITENODE_TEMPLATE_CHILDREN_FOR_SPECIALIZATION,
             int s0, int s1, int s2, int s3, int s4, int s5, int s6, int s7, int s8, int s9>
    class CompositeGridFunctionSpace<GridFunctionSpaceComponentBlockwiseMapper<s0,s1,s2,s3,s4,s5,s6,s7,s8,s9>,
                                     DUNE_TYPETREE_COMPOSITENODE_CHILDTYPES>
      : public DUNE_TYPETREE_COMPOSITENODE_BASETYPE
      , public PowerCompositeGridFunctionSpaceBase<CompositeGridFunctionSpace<
                                                     GridFunctionSpaceComponentBlockwiseMapper<s0,s1,s2,s3,s4,s5,s6,s7,s8,s9>,
                                                     DUNE_TYPETREE_COMPOSITENODE_CHILDTYPES>,
                                                   DUNE_TYPETREE_COMPOSITENODE_BASETYPE,
                                                   GridFunctionSpaceComponentBlockwiseMapper<s0,s1,s2,s3,s4,s5,s6,s7,s8,s9>
                                                   >
    {
      typedef GridFunctionSpaceComponentBlockwiseMapper<s0,s1,s2,s3,s4,s5,s6,s7,s8,s9> BlockwiseMapper;
      typedef DUNE_TYPETREE_COMPOSITENODE_BASETYPE BaseT;

      typedef PowerCompositeGridFunctionSpaceBase<CompositeGridFunctionSpace,BaseT,BlockwiseMapper> ImplementationBase;

      friend class PowerCompositeGridFunctionSpaceBase<CompositeGridFunctionSpace,BaseT,BlockwiseMapper>;

    public:

      typedef CompositeGridFunctionSpaceTag ImplementationTag;

      //! export traits class
      typedef typename ImplementationBase::Traits Traits;

      // define local function space parametrized by self
      typedef typename Dune::PDELab::TypeTree::TransformTree<CompositeGridFunctionSpace,gfs_to_lfs>::Type LocalFunctionSpace;

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
        return subMap(i,j);
      }

      typename Traits::SizeType subMap (typename Traits::SizeType i, typename Traits::SizeType j) const
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


    private:

      using ImplementationBase::childLocalSize;
      using ImplementationBase::childGlobalSize;
      using ImplementationBase::maxlocalsize;
      using ImplementationBase::offset;

      void calculateSizes ()
      {
        Dune::dinfo << "CompositeGridFunctionSpace(blockwise version):"
                    << std::endl;

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
      }

    };

    template<DUNE_TYPETREE_COMPOSITENODE_TEMPLATE_CHILDREN_FOR_SPECIALIZATION>
    class CompositeGridFunctionSpace<GridFunctionSpaceBlockwiseMapper,
                                     DUNE_TYPETREE_COMPOSITENODE_CHILDTYPES>
      : public CompositeGridFunctionSpace<GridFunctionSpaceComponentBlockwiseMapper<1>,
                                          DUNE_TYPETREE_COMPOSITENODE_CHILDTYPES>
    {

      typedef CompositeGridFunctionSpace<GridFunctionSpaceComponentBlockwiseMapper<1>,
                                         DUNE_TYPETREE_COMPOSITENODE_CHILDTYPES> BaseT;

    public:

      CompositeGridFunctionSpace(DUNE_TYPETREE_COMPOSITENODE_CONSTRUCTOR_SIGNATURE)
        : BaseT(DUNE_TYPETREE_COMPOSITENODE_CHILDVARIABLES)
      {}

    };


    /**
        \brief Tupel of grid function spaces base class that holds
        implementation of the methods specialization for dynamic
        blockwise ordering
    */
    template<DUNE_TYPETREE_COMPOSITENODE_TEMPLATE_CHILDREN_FOR_SPECIALIZATION>
    class CompositeGridFunctionSpace<GridFunctionSpaceDynamicBlockwiseMapper,DUNE_TYPETREE_COMPOSITENODE_CHILDTYPES>
      : public DUNE_TYPETREE_COMPOSITENODE_BASETYPE
      , public PowerCompositeGridFunctionSpaceBase<CompositeGridFunctionSpace<
                                                     GridFunctionSpaceDynamicBlockwiseMapper,
                                                     DUNE_TYPETREE_COMPOSITENODE_CHILDTYPES>,
                                                   DUNE_TYPETREE_COMPOSITENODE_BASETYPE,
                                                   GridFunctionSpaceDynamicBlockwiseMapper
                                                   >
    {
      typedef GridFunctionSpaceDynamicBlockwiseMapper BlockwiseMapper;
      typedef DUNE_TYPETREE_COMPOSITENODE_BASETYPE BaseT;

      typedef PowerCompositeGridFunctionSpaceBase<CompositeGridFunctionSpace,BaseT,BlockwiseMapper> ImplementationBase;

      friend class PowerCompositeGridFunctionSpaceBase<CompositeGridFunctionSpace,BaseT,BlockwiseMapper>;

    public:

      typedef CompositeGridFunctionSpaceTag ImplementationTag;

      //! export traits class
      typedef typename ImplementationBase::Traits Traits;

      // define local function space parametrized by self
      typedef typename Dune::PDELab::TypeTree::TransformTree<CompositeGridFunctionSpace,gfs_to_lfs>::Type LocalFunctionSpace;

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
        return subMap(i,j);
      }

      typename Traits::SizeType subMap (typename Traits::SizeType i, typename Traits::SizeType j) const
      {
        BlockIndexRangeIterator & it = blockIndexIterators[i];
        while(it->first < j)++it;
        while(it->first > j)--it;
        const typename Traits::SizeType global_index = it->second + j - it->first;
        return global_index;
      }


    private:

      using ImplementationBase::childLocalSize;
      using ImplementationBase::childGlobalSize;
      using ImplementationBase::maxlocalsize;
      using ImplementationBase::offset;

      void calculateSizes ()
      {
        Dune::dinfo << "CompositeGridFunctionSpace(blockwise version):"
                    << std::endl;

        typedef typename Traits::GridViewType GridView;
        const GridView & gv = this->gridview();

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
          DynamicBlockwiseMapperImp::GetChildOffsetsMetaProgram<CompositeGridFunctionSpace,BaseT::CHILDREN,0>::
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
      }

    private:
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

  }

}
//! \}
#endif
