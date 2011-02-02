// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_POWERGRIDFUNCTIONSPACE_HH
#define DUNE_PDELAB_POWERGRIDFUNCTIONSPACE_HH

#include "powercompositegridfunctionspacebase.hh"

namespace Dune {
  namespace PDELab {

    //=======================================
    // power grid function space
    //=======================================

    /** \brief base class for tuples of grid function spaces
        product of identical grid function spaces
        base class that holds implementation of the methods

        PGFS(T,k) = {T}^k

        \tparam T the underlying are all grid function spaces
        \tparam k power factor
        \tparam Mapper is the ordering parameter. Use e.g.
        \link GridFunctionSpaceLexicographicMapper GridFunctionSpaceLexicographicMapper \endlink
        or \link  GridFunctionSpaceComponentBlockwiseMapper  GridFunctionSpaceComponentBlockwiseMapper \endlink
        or \link  GridFunctionSpaceBlockwiseMapper  GridFunctionSpaceBlockwiseMapper \endlink
        or \link  GridFunctionSpaceDynamicBlockwiseMapper  GridFunctionSpaceDynamicBlockwiseMapper \endlink
    */
    template<typename T, std::size_t k, typename Mapper = GridFunctionSpaceLexicographicMapper>
    class PowerGridFunctionSpace
    {
    private:
      dune_static_assert(AlwaysFalse<Mapper>::value, "You seem to be using an unsupported Mapper");
    };

    // Specialization for GridFunctionSpaceLexicographicMapper
    template<typename T, std::size_t k>
    class PowerGridFunctionSpace<T,k,GridFunctionSpaceLexicographicMapper>
      : public TypeTree::PowerNode<T,k>
      , public PowerCompositeGridFunctionSpaceBase<PowerGridFunctionSpace<
                                                     T,
                                                     k,
                                                     GridFunctionSpaceLexicographicMapper>,
                                                   typename T::Traits::GridViewType,
                                                   typename T::Traits::BackendType,
                                                   GridFunctionSpaceLexicographicMapper,
                                                   k
                                                   >
    {

      typedef TypeTree::PowerNode<T,k> BaseT;

      typedef PowerCompositeGridFunctionSpaceBase<PowerGridFunctionSpace,
                                                  typename T::Traits::GridViewType,
                                                  typename T::Traits::BackendType,
                                                  GridFunctionSpaceLexicographicMapper,
                                                  k> ImplementationBase;

      friend class PowerCompositeGridFunctionSpaceBase<PowerGridFunctionSpace,
                                                       typename T::Traits::GridViewType,
                                                       typename T::Traits::BackendType,
                                                       GridFunctionSpaceLexicographicMapper,
                                                       k>;

    public:

      typedef PowerGridFunctionSpaceTag ImplementationTag;

      //! export traits class
      typedef typename ImplementationBase::Traits Traits;

      // define local function space parametrized by self
      typedef typename Dune::PDELab::TypeTree::TransformTree<PowerGridFunctionSpace,gfs_to_lfs>::Type LocalFunctionSpace;

      PowerGridFunctionSpace(T& c)
        : BaseT(c)
      {
        this->setup();
      }

      template<std::size_t K = 2>
      PowerGridFunctionSpace (typename enable_if<K == BaseT::CHILDREN,T>::type& c0,
                              T& c1)
        : BaseT(c0,c1)
      {
        this->setup();
      }

      template<std::size_t K = 3>
      PowerGridFunctionSpace (typename enable_if<K == BaseT::CHILDREN,T>::type& c0,
                              T& c1,
                              T& c2)
        : BaseT(c0,c1,c2)
      {
        this->setup();
      }

      template<std::size_t K = 4>
      PowerGridFunctionSpace (typename enable_if<K == BaseT::CHILDREN,T>::type& c0,
                              T& c1,
                              T& c2,
                              T& c3)
        : BaseT(c0,c1,c2,c3)
      {
        this->setup();
      }

      template<std::size_t K = 5>
      PowerGridFunctionSpace (typename enable_if<K == BaseT::CHILDREN,T>::type& c0,
                              T& c1,
                              T& c2,
                              T& c3,
                              T& c4)
        : BaseT(c0,c1,c2,c3,c4)
      {
        this->setup();
      }

      template<std::size_t K = 6>
      PowerGridFunctionSpace (typename enable_if<K == BaseT::CHILDREN,T>::type& c0,
                              T& c1,
                              T& c2,
                              T& c3,
                              T& c4,
                              T& c5)
        : BaseT(c0,c1,c2,c3,c4,c5)
      {
        this->setup();
      }

      template<std::size_t K = 7>
      PowerGridFunctionSpace (typename enable_if<K == BaseT::CHILDREN,T>::type& c0,
                              T& c1,
                              T& c2,
                              T& c3,
                              T& c4,
                              T& c5,
                              T& c6)
        : BaseT(c0,c1,c2,c3,c4,c5,c6)
      {
        this->setup();
      }

      template<std::size_t K = 8>
      PowerGridFunctionSpace (typename enable_if<K == BaseT::CHILDREN,T>::type& c0,
                              T& c1,
                              T& c2,
                              T& c3,
                              T& c4,
                              T& c5,
                              T& c6,
                              T& c7)
        : BaseT(c0,c1,c2,c3,c4,c5,c6,c7)
      {
        this->setup();
      }

      template<std::size_t K = 9>
      PowerGridFunctionSpace (typename enable_if<K == BaseT::CHILDREN,T>::type& c0,
                              T& c1,
                              T& c2,
                              T& c3,
                              T& c4,
                              T& c5,
                              T& c6,
                              T& c7,
                              T& c8)
        : BaseT(c0,c1,c2,c3,c4,c5,c6,c7,c8)
      {
        this->setup();
      }

      template<std::size_t K = 10>
      PowerGridFunctionSpace (typename enable_if<K == BaseT::CHILDREN,T>::type& c0,
                              T& c1,
                              T& c2,
                              T& c3,
                              T& c4,
                              T& c5,
                              T& c6,
                              T& c7,
                              T& c8,
                              T& c9)
        : BaseT(c0,c1,c2,c3,c4,c5,c6,c7,c8,c9)
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
        Dune::dinfo << "PowerGridFunctionSpace(lexicographic version):"
                    << std::endl;
        Dune::dinfo << "( ";
        offset[0] = 0;
        maxlocalsize = 0;
        for (std::size_t i=0; i<k; i++)
          {
            Dune::dinfo << childGlobalSize[i] << " ";
            offset[i+1] = offset[i]+childGlobalSize[i];
            maxlocalsize += childLocalSize[i];
          }
        Dune::dinfo << ") total size = " << offset[k]
                    << " max local size = " << maxlocalsize
                    << std::endl;
      }

    };


    // product of identical grid function spaces
    // base class that holds implementation of the methods
    // specialization for blockwise ordering
    template<typename T, std::size_t k, int s>
    class PowerGridFunctionSpace<T,k,GridFunctionSpaceComponentBlockwiseMapper<s> >
      : public TypeTree::PowerNode<T,k>
      , public PowerCompositeGridFunctionSpaceBase<PowerGridFunctionSpace<
                                                     T,
                                                     k,
                                                     GridFunctionSpaceComponentBlockwiseMapper<s> >,
                                                   typename T::Traits::GridViewType,
                                                   typename T::Traits::BackendType,
                                                   GridFunctionSpaceComponentBlockwiseMapper<s>,
                                                   k
                                                   >
    {

      typedef TypeTree::PowerNode<T,k> BaseT;

      typedef GridFunctionSpaceComponentBlockwiseMapper<s> BlockwiseMapper;

      typedef PowerCompositeGridFunctionSpaceBase<PowerGridFunctionSpace,
                                                  typename T::Traits::GridViewType,
                                                  typename T::Traits::BackendType,
                                                  BlockwiseMapper,
                                                  k> ImplementationBase;

      friend class PowerCompositeGridFunctionSpaceBase<PowerGridFunctionSpace,
                                                       typename T::Traits::GridViewType,
                                                       typename T::Traits::BackendType,
                                                       BlockwiseMapper,
                                                       k>;

    public:

      typedef PowerGridFunctionSpaceTag ImplementationTag;

      //! export traits class
      typedef typename ImplementationBase::Traits Traits;

      // define local function space parametrized by self
      typedef typename Dune::PDELab::TypeTree::TransformTree<PowerGridFunctionSpace,gfs_to_lfs>::Type LocalFunctionSpace;

      PowerGridFunctionSpace(T& c)
        : BaseT(c)
      {
        this->setup();
      }

      template<std::size_t K = 2>
      PowerGridFunctionSpace (typename enable_if<K == BaseT::CHILDREN,T>::type& c0,
                              T& c1)
        : BaseT(c0,c1)
      {
        this->setup();
      }

      template<std::size_t K = 3>
      PowerGridFunctionSpace (typename enable_if<K == BaseT::CHILDREN,T>::type& c0,
                              T& c1,
                              T& c2)
        : BaseT(c0,c1,c2)
      {
        this->setup();
      }

      template<std::size_t K = 4>
      PowerGridFunctionSpace (typename enable_if<K == BaseT::CHILDREN,T>::type& c0,
                              T& c1,
                              T& c2,
                              T& c3)
        : BaseT(c0,c1,c2,c3)
      {
        this->setup();
      }

      template<std::size_t K = 5>
      PowerGridFunctionSpace (typename enable_if<K == BaseT::CHILDREN,T>::type& c0,
                              T& c1,
                              T& c2,
                              T& c3,
                              T& c4)
        : BaseT(c0,c1,c2,c3,c4)
      {
        this->setup();
      }

      template<std::size_t K = 6>
      PowerGridFunctionSpace (typename enable_if<K == BaseT::CHILDREN,T>::type& c0,
                              T& c1,
                              T& c2,
                              T& c3,
                              T& c4,
                              T& c5)
        : BaseT(c0,c1,c2,c3,c4,c5)
      {
        this->setup();
      }

      template<std::size_t K = 7>
      PowerGridFunctionSpace (typename enable_if<K == BaseT::CHILDREN,T>::type& c0,
                              T& c1,
                              T& c2,
                              T& c3,
                              T& c4,
                              T& c5,
                              T& c6)
        : BaseT(c0,c1,c2,c3,c4,c5,c6)
      {
        this->setup();
      }

      template<std::size_t K = 8>
      PowerGridFunctionSpace (typename enable_if<K == BaseT::CHILDREN,T>::type& c0,
                              T& c1,
                              T& c2,
                              T& c3,
                              T& c4,
                              T& c5,
                              T& c6,
                              T& c7)
        : BaseT(c0,c1,c2,c3,c4,c5,c6,c7)
      {
        this->setup();
      }

      template<std::size_t K = 9>
      PowerGridFunctionSpace (typename enable_if<K == BaseT::CHILDREN,T>::type& c0,
                              T& c1,
                              T& c2,
                              T& c3,
                              T& c4,
                              T& c5,
                              T& c6,
                              T& c7,
                              T& c8)
        : BaseT(c0,c1,c2,c3,c4,c5,c6,c7,c8)
      {
        this->setup();
      }

      template<std::size_t K = 10>
      PowerGridFunctionSpace (typename enable_if<K == BaseT::CHILDREN,T>::type& c0,
                              T& c1,
                              T& c2,
                              T& c3,
                              T& c4,
                              T& c5,
                              T& c6,
                              T& c7,
                              T& c8,
                              T& c9)
        : BaseT(c0,c1,c2,c3,c4,c5,c6,c7,c8,c9)
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
        return (j%s)+(j/s)*k*s+i*s;
      }


    private:

      using ImplementationBase::childLocalSize;
      using ImplementationBase::childGlobalSize;
      using ImplementationBase::maxlocalsize;
      using ImplementationBase::offset;

      void calculateSizes ()
      {
        Dune::dinfo << "PowerGridFunctionSpace(blockwise version):"
                    << std::endl;
        Dune::dinfo << "( ";
        offset[0] = 0;
        maxlocalsize = 0;
        for (std::size_t i=0; i<k; i++)
          {
            offset[i+1] = offset[i]+childGlobalSize[i];
            Dune::dinfo << childGlobalSize[i] << "[" << offset[i] << "] ";
            maxlocalsize += childLocalSize[i];
          }
        Dune::dinfo << ") total size = " << offset[k]
                    << " max local size = " << maxlocalsize
                    << std::endl;
        /* check the local block size */
        if (childLocalSize[0]%s != 0)
          DUNE_THROW(Exception,
                     "number of DOFs (" << childLocalSize[0] << ") per component "
                     "must be a multiple of the BlockSize (" << s << ")");
        for (std::size_t i=1; i<k; i++)
          if (childLocalSize[i]!=childLocalSize[0])
            DUNE_THROW(Exception, "components must be of equal size");
      }

    };

    template<typename T, std::size_t k>
    class PowerGridFunctionSpace<T,k,GridFunctionSpaceBlockwiseMapper >
      : public PowerGridFunctionSpace<T,k,GridFunctionSpaceComponentBlockwiseMapper<1> >
    {

      typedef PowerGridFunctionSpace<T,k,GridFunctionSpaceComponentBlockwiseMapper<1> > BaseT;

    public:

      typedef typename Dune::PDELab::TypeTree::TransformTree<PowerGridFunctionSpace,gfs_to_lfs>::Type LocalFunctionSpace;

      PowerGridFunctionSpace(T& c)
        : BaseT(c)
      {
      }

      template<std::size_t K = 2>
      PowerGridFunctionSpace (typename enable_if<K == BaseT::CHILDREN,T>::type& c0,
                              T& c1)
        : BaseT(c0,c1)
      {
      }

      template<std::size_t K = 3>
      PowerGridFunctionSpace (typename enable_if<K == BaseT::CHILDREN,T>::type& c0,
                              T& c1,
                              T& c2)
        : BaseT(c0,c1,c2)
      {
      }

      template<std::size_t K = 4>
      PowerGridFunctionSpace (typename enable_if<K == BaseT::CHILDREN,T>::type& c0,
                              T& c1,
                              T& c2,
                              T& c3)
        : BaseT(c0,c1,c2,c3)
      {
      }

      template<std::size_t K = 5>
      PowerGridFunctionSpace (typename enable_if<K == BaseT::CHILDREN,T>::type& c0,
                              T& c1,
                              T& c2,
                              T& c3,
                              T& c4)
        : BaseT(c0,c1,c2,c3,c4)
      {
      }

      template<std::size_t K = 6>
      PowerGridFunctionSpace (typename enable_if<K == BaseT::CHILDREN,T>::type& c0,
                              T& c1,
                              T& c2,
                              T& c3,
                              T& c4,
                              T& c5)
        : BaseT(c0,c1,c2,c3,c4,c5)
      {
      }

      template<std::size_t K = 7>
      PowerGridFunctionSpace (typename enable_if<K == BaseT::CHILDREN,T>::type& c0,
                              T& c1,
                              T& c2,
                              T& c3,
                              T& c4,
                              T& c5,
                              T& c6)
        : BaseT(c0,c1,c2,c3,c4,c5,c6)
      {
      }

      template<std::size_t K = 8>
      PowerGridFunctionSpace (typename enable_if<K == BaseT::CHILDREN,T>::type& c0,
                              T& c1,
                              T& c2,
                              T& c3,
                              T& c4,
                              T& c5,
                              T& c6,
                              T& c7)
        : BaseT(c0,c1,c2,c3,c4,c5,c6,c7)
      {
      }

      template<std::size_t K = 9>
      PowerGridFunctionSpace (typename enable_if<K == BaseT::CHILDREN,T>::type& c0,
                              T& c1,
                              T& c2,
                              T& c3,
                              T& c4,
                              T& c5,
                              T& c6,
                              T& c7,
                              T& c8)
        : BaseT(c0,c1,c2,c3,c4,c5,c6,c7,c8)
      {
      }

      template<std::size_t K = 10>
      PowerGridFunctionSpace (typename enable_if<K == BaseT::CHILDREN,T>::type& c0,
                              T& c1,
                              T& c2,
                              T& c3,
                              T& c4,
                              T& c5,
                              T& c6,
                              T& c7,
                              T& c8,
                              T& c9)
        : BaseT(c0,c1,c2,c3,c4,c5,c6,c7,c8,c9)
      {
      }

    };

    /**
        \brief Tupel of grid function spaces base class that holds
        implementation of the methods specialization for dynamic
        blockwise ordering
    */
    template<typename T, std::size_t k>
    class PowerGridFunctionSpace<T,k,GridFunctionSpaceDynamicBlockwiseMapper >
      : public TypeTree::PowerNode<T,k>
      , public PowerCompositeGridFunctionSpaceBase<PowerGridFunctionSpace<
                                                     T,
                                                     k,
                                                     GridFunctionSpaceDynamicBlockwiseMapper>,
                                                   typename T::Traits::GridViewType,
                                                   typename T::Traits::BackendType,
                                                   GridFunctionSpaceDynamicBlockwiseMapper,
                                                   k
                                                   >
    {

      typedef TypeTree::PowerNode<T,k> BaseT;

      typedef PowerCompositeGridFunctionSpaceBase<PowerGridFunctionSpace,
                                                  typename T::Traits::GridViewType,
                                                  typename T::Traits::BackendType,
                                                  GridFunctionSpaceDynamicBlockwiseMapper,
                                                  k> ImplementationBase;

      friend class PowerCompositeGridFunctionSpaceBase<PowerGridFunctionSpace,
                                                       typename T::Traits::GridViewType,
                                                       typename T::Traits::BackendType,
                                                       GridFunctionSpaceDynamicBlockwiseMapper,
                                                       k>;

    public:

      typedef PowerGridFunctionSpaceTag ImplementationTag;

      //! export traits class
      typedef typename ImplementationBase::Traits Traits;

      // define local function space parametrized by self
      typedef typename Dune::PDELab::TypeTree::TransformTree<PowerGridFunctionSpace,gfs_to_lfs>::Type LocalFunctionSpace;

      PowerGridFunctionSpace(T& c)
        : BaseT(c)
      {
        this->setup();
      }

      template<std::size_t K = 2>
      PowerGridFunctionSpace (typename enable_if<K == BaseT::CHILDREN,T>::type& c0,
                              T& c1)
        : BaseT(c0,c1)
      {
        this->setup();
      }

      template<std::size_t K = 3>
      PowerGridFunctionSpace (typename enable_if<K == BaseT::CHILDREN,T>::type& c0,
                              T& c1,
                              T& c2)
        : BaseT(c0,c1,c2)
      {
        this->setup();
      }

      template<std::size_t K = 4>
      PowerGridFunctionSpace (typename enable_if<K == BaseT::CHILDREN,T>::type& c0,
                              T& c1,
                              T& c2,
                              T& c3)
        : BaseT(c0,c1,c2,c3)
      {
        this->setup();
      }

      template<std::size_t K = 5>
      PowerGridFunctionSpace (typename enable_if<K == BaseT::CHILDREN,T>::type& c0,
                              T& c1,
                              T& c2,
                              T& c3,
                              T& c4)
        : BaseT(c0,c1,c2,c3,c4)
      {
        this->setup();
      }

      template<std::size_t K = 6>
      PowerGridFunctionSpace (typename enable_if<K == BaseT::CHILDREN,T>::type& c0,
                              T& c1,
                              T& c2,
                              T& c3,
                              T& c4,
                              T& c5)
        : BaseT(c0,c1,c2,c3,c4,c5)
      {
        this->setup();
      }

      template<std::size_t K = 7>
      PowerGridFunctionSpace (typename enable_if<K == BaseT::CHILDREN,T>::type& c0,
                              T& c1,
                              T& c2,
                              T& c3,
                              T& c4,
                              T& c5,
                              T& c6)
        : BaseT(c0,c1,c2,c3,c4,c5,c6)
      {
        this->setup();
      }

      template<std::size_t K = 8>
      PowerGridFunctionSpace (typename enable_if<K == BaseT::CHILDREN,T>::type& c0,
                              T& c1,
                              T& c2,
                              T& c3,
                              T& c4,
                              T& c5,
                              T& c6,
                              T& c7)
        : BaseT(c0,c1,c2,c3,c4,c5,c6,c7)
      {
        this->setup();
      }

      template<std::size_t K = 9>
      PowerGridFunctionSpace (typename enable_if<K == BaseT::CHILDREN,T>::type& c0,
                              T& c1,
                              T& c2,
                              T& c3,
                              T& c4,
                              T& c5,
                              T& c6,
                              T& c7,
                              T& c8)
        : BaseT(c0,c1,c2,c3,c4,c5,c6,c7,c8)
      {
        this->setup();
      }

      template<std::size_t K = 10>
      PowerGridFunctionSpace (typename enable_if<K == BaseT::CHILDREN,T>::type& c0,
                              T& c1,
                              T& c2,
                              T& c3,
                              T& c4,
                              T& c5,
                              T& c6,
                              T& c7,
                              T& c8,
                              T& c9)
        : BaseT(c0,c1,c2,c3,c4,c5,c6,c7,c8,c9)
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
        return this->subMap(i,j);
      }

      //! map index from child i's index set into our index set
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

        typedef typename Traits::GridViewType GridView;
        const GridView & gv = this->gridview();

        // Initialize offset array for each child
        blockIndices.clear();
        blockIndices.resize(BaseT::CHILDREN);

        std::vector<std::vector<typename Traits::SizeType> > childOffsets;
        childOffsets.resize(BaseT::CHILDREN);
        for(std::size_t i=0; i<BaseT::CHILDREN; ++i)
          childOffsets[i].resize(gv.size(0)+1);


        // Iterate grid and determine the offsets for each child
        typedef typename GridView::template Codim<0>::Iterator Iterator;
        Iterator it = gv.template begin<0>();
        const Iterator eit = gv.template end<0>();
        typename Traits::SizeType running_index(0);
        for(; it!=eit; ++it){
          typename Traits::SizeType e_index = gv.indexSet().index(*it);

          // Loop over children (realized by meta-program)
          DynamicBlockwiseMapperImp::GetChildOffsetsMetaProgram<PowerGridFunctionSpace,BaseT::CHILDREN,0>::
            getChildOffsets(*this,*it,childOffsets);

          for(std::size_t i=0; i<BaseT::CHILDREN; ++i){
            if(childOffsets[i][e_index+1] != childOffsets[i][e_index]){
              // Add new block index range element to list
              blockIndices[i].push_back(SizeTypePair(childOffsets[i][e_index],running_index));

              // Update running index
              running_index += childOffsets[i][e_index+1] - childOffsets[i][e_index];
            }
          }
        }

        // Insert a "stop entry" at end of every list
        for(std::size_t i=0; i<BaseT::CHILDREN; ++i)
          blockIndices[i].push_back(SizeTypePair(childOffsets[i][gv.size(0)],running_index));


        blockIndexIterators.clear();
        for(std::size_t i=0; i<BaseT::CHILDREN; ++i)
          blockIndexIterators.push_back(blockIndices[i].begin());


        Dune::dinfo << "PowerGridFunctionSpace(blockwise version):"
                    << std::endl;
        Dune::dinfo << "( ";
        offset[0] = 0;
        maxlocalsize = 0;
        for (std::size_t i=0; i<k; i++)
          {
            offset[i+1] = offset[i]+childGlobalSize[i];
            Dune::dinfo << childGlobalSize[i] << "[" << offset[i] << "] ";
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

    /// \brief product of identical grid function spaces
    ///
    /// the specializations of this class just set the members
    /// all the methods are generic in the implementation
    /// \tparam T The type of the underlying grid function space
    /// \tparam k how many identical function space to use in the composition
    /// \tparam P The type of the mapper used. The mapper maps each degree of freedom
    /// of each function space to a unique index. Use e.g.
    /// \link GridFunctionSpaceLexicographicMapper GridFunctionSpaceLexicographicMapper \endlink
    /// or \link  GridFunctionSpaceComponentBlockwiseMapper  GridFunctionSpaceComponentBlockwiseMapper \endlink
    /// or \link  GridFunctionSpaceBlockwiseMapper  GridFunctionSpaceBlockwiseMapper \endlink



  }
}
#endif
