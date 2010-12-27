// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_LOCALFUNCTIONSPACE_HH
#define DUNE_PDELAB_LOCALFUNCTIONSPACE_HH

#include<vector>

#include <dune/common/stdstreams.hh>

#include <dune/localfunctions/common/interfaceswitch.hh>

#include "../common/multitypetree.hh"
#include "../common/cpstoragepolicy.hh"

#include "localindex.hh"

namespace Dune {
  namespace PDELab {

    //! \addtogroup GridFunctionSpace
    //! \ingroup PDELab
    //! \{

    //=======================================
    // local function space base: metaprograms
    //=======================================

    namespace {
      
      template<typename T, bool isleaf, typename E, typename GC, typename Int = typename GC::size_type>
      struct LocalFunctionSpaceBaseVisitNodeMetaProgram;

      template<typename T, typename E, typename GC, typename Int, int n, int i>
      struct LocalFunctionSpaceBaseVisitChildMetaProgram // visit i'th child of inner node
      {
        typedef typename T::template Child<i>::Type C;
        LocalFunctionSpaceBaseVisitNodeMetaProgram<C,C::isLeaf,E,GC,Int> childTMP;
        LocalFunctionSpaceBaseVisitChildMetaProgram<T,E,GC,Int,n,i+1> siblingTMP;
        void fill_indices (T& t, const E& e, Int& offset, GC * const global)
        {
          // vist children of node t in order
          Int initial_offset = offset; // remember initial offset to compute size later
          childTMP.fill_indices(t.template getChild<i>(),e,offset,global);
          for (Int j=initial_offset; j<offset; j++)
            (*global)[j] = t.pgfs->template subMap<i>((*global)[j]);
          // visit siblings
          siblingTMP.fill_indices(t,e,offset,global);
        }
        void compute_size (T& t, const E& e, Int& size)
        {
          if (i == 0) // braucht man dieses if?
            t.offset = size;
          // now we are at a multi component local function space
          // vist children of node t in order
          Int initial_size = size; // remember initial offset to compute size later
          childTMP.compute_size(t.template getChild<i>(),e,size);
          t.n = size-initial_size;
          // visit siblings
          siblingTMP.compute_size(t,e,size);
        }
      };

      template<typename T, typename E, typename GC, typename Int, int n>
      struct LocalFunctionSpaceBaseVisitChildMetaProgram<T,E,GC,Int,n,n> // end of child recursion
      {
        void fill_indices (T& t, const E& e, Int& offset, GC * const global)
        {
          return;
        }
        void compute_size (T& t, const E& e, Int& offset)
        {
          return;
        }
      };

      template<typename T, bool isleaf, typename E, typename GC, typename Int> 
      struct LocalFunctionSpaceBaseVisitNodeMetaProgram // visit inner node
      {
        LocalFunctionSpaceBaseVisitChildMetaProgram<T,E,GC,Int,T::CHILDREN,0> childTMP;
        void fill_indices (T& t, const E& e, Int& offset, GC * const global)
        {
          // std::cout << "OFFSET: " << t.offset << " SIZE: " << t.n << std::endl;
          // now we are at a multi component local function space
          t.global = global;
          assert(t.offset == offset);
          childTMP.fill_indices(t,e,offset,global);
        }
        void compute_size (T& t, const E& e, Int& size)
        {
          t.offset = size;
          // now we are at a multi component local function space
          Int initial_size = size; // remember initial offset to compute size later
          childTMP.compute_size(t,e,size);
          t.n = size-initial_size;
        }
      };

      template<typename T, typename E, typename GC, typename Int> 
      struct LocalFunctionSpaceBaseVisitNodeMetaProgram<T,true,E,GC,Int> // visit leaf node 
      {
        std::vector<typename T::Traits::GridFunctionSpaceType::Traits::SizeType> _global;
        void fill_indices (T& t, const E& e, Int& offset, GC * const global)
        {
          // now we are at a single component local function space
          // which is part of a multi component local function space
          t.global = global;
          assert(t.offset == offset);
          _global.resize(t.n);
          t.pgfs->globalIndices(*(t.pfe),e,_global); // get global indices for this finite element
          for (Int i=0; i<t.n; i++) (*global)[offset+i]=_global[i]; 
          offset += t.n; // append this chunk
        }
        void compute_size (T& t, const E& e, Int& size)
        {
          t.offset = size;
          // now we are at a single component local function space
          // which is part of a multi component local function space
          T::FESwitch::setStore(t.pfe, t.pgfs->finiteElementMap().find(e));
          // determine size of this chunk
          t.n = T::FESwitch::basis(*t.pfe).size();
          size += t.n; // append this chunk
        }
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
    public:
      typedef LocalFunctionSpaceBaseTraits<GFS> Traits;

      //! \brief empty constructor (needed for CopyStoragePolicy)
      LocalFunctionSpaceBaseNode () {}
      
      //! \brief construct from global function space
      LocalFunctionSpaceBaseNode (const GFS& gfs) : 
        pgfs(&gfs), global_storage(gfs.maxLocalSize()), global(0), n(0)
      {}
      
      //! \brief initialize with grid function space
      void setup (const GFS& gfs)
      {
        pgfs = &gfs;
      }

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

      CountingPointer<GFS const> pgfs;
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
      static
        LocalFunctionSpaceBaseVisitNodeMetaProgram<NodeType,NodeType::isLeaf,
                                                    typename Traits::Element,
                                                    typename Traits::IndexContainer>
        TMP;
      // we should only call bind on out selfs
      assert(&node == this);
      
      // make offset
      typename Traits::IndexContainer::size_type size=0;
      
      // compute sizes
      TMP.compute_size(node,e,size);
      assert(size == n);
      
      // initialize iterators and fill indices
      size = 0;
      TMP.fill_indices(node,e,size,&global_storage);
      assert(global == &global_storage);
      assert(size <= global_storage.size());
      assert(offset == 0);
      
      // apply upMap
      for (typename Traits::IndexContainer::size_type i=0; i<size; i++)
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

    template<typename GFS>
    class PowerLocalFunctionSpaceNode;

    // local function space for a power grid function space
    template<typename GFS>
    class PowerLocalFunctionSpaceNode :
      public LocalFunctionSpaceBaseNode<GFS>,
      public PowerNode<typename GFS::template Child<0>::Type::LocalFunctionSpace::Traits::NodeType,
                       GFS::CHILDREN,CopyStoragePolicy>
    {
      typedef LocalFunctionSpaceBaseNode<GFS> BaseT;

      // friend decl for bind meta program
      template<typename T, bool b, typename E, typename GC, typename Int> 
      friend struct LocalFunctionSpaceBaseVisitNodeMetaProgram;
      template<typename T, typename E, typename GC, typename Int, int n, int i>
      friend struct LocalFunctionSpaceBaseVisitChildMetaProgram;

    public:
      typedef PowerCompositeLocalFunctionSpaceTraits<GFS,PowerLocalFunctionSpaceNode> Traits;

      //! \brief empty constructor (needed for CopyStoragePolicy)
      PowerLocalFunctionSpaceNode ()
      {
      }

      //! \brief initialize with grid function space
      PowerLocalFunctionSpaceNode (const GFS& gfs)  : 
        BaseT(gfs)
      {
        setup(gfs);
      }

      //! \brief initialize with grid function space
      void setup (const GFS& gfs)
      {
        BaseT::setup(gfs);
        for (int i=0; i<GFS::CHILDREN; i++)
          {
            dinfo << "setting up child " << i << " of " << GFS::CHILDREN << std::endl;
            this->getChild(i).setup(this->pgfs->getChild(i));
          }
      }

      //! \brief bind local function space to entity
      void bind (const typename Traits::Element& e)
      {
        // call method on base class, this avoid the barton neckman trick
        BaseT::bind(*this,e);
      }

    };

    //=======================================
    // local function space base: composite implementation
    //=======================================

#ifndef DOXYGEN
    namespace {
      template<typename T, int n, int i>
      struct CompositeLocalFunctionSpaceNodeVisitChildMetaProgram // visit child of inner node
      {
        template<typename GFS>
        static void setup (T& t, const GFS& gfs)
        {
          dinfo << "setting up child " << i << " of " << n << std::endl;
          t.template getChild<i>().setup(gfs.template getChild<i>());
          CompositeLocalFunctionSpaceNodeVisitChildMetaProgram<T,n,i+1>::
            setup(t,gfs);
        }
      };
      
      template<typename T, int n>
      struct CompositeLocalFunctionSpaceNodeVisitChildMetaProgram<T,n,n> // end of child recursion
      {
        template<typename GFS>
        static void setup (T& t, const GFS& gfs)
        {
        }
      };
    }
#endif

    template<typename GFS, int k>
    struct CompositeLocalFunctionSpaceNodeBaseTypeTraits
    {
#ifdef DOXYGEN
      typedef NodeType Type;
#endif
    };

#ifndef DOXYGEN
    template<typename GFS>
    struct CompositeLocalFunctionSpaceNodeBaseTypeTraits<GFS,2>
    {
      typedef CompositeNode<CopyStoragePolicy, 
                            typename GFS::template Child<0>::Type::LocalFunctionSpace::Traits::NodeType,
                            typename GFS::template Child<1>::Type::LocalFunctionSpace::Traits::NodeType> Type;
    };

    template<typename GFS>
    struct CompositeLocalFunctionSpaceNodeBaseTypeTraits<GFS,3>
    {
      typedef CompositeNode<CopyStoragePolicy, 
                            typename GFS::template Child<0>::Type::LocalFunctionSpace::Traits::NodeType,
                            typename GFS::template Child<1>::Type::LocalFunctionSpace::Traits::NodeType,
                            typename GFS::template Child<2>::Type::LocalFunctionSpace::Traits::NodeType> Type;
    };

    template<typename GFS>
    struct CompositeLocalFunctionSpaceNodeBaseTypeTraits<GFS,4>
    {
      typedef CompositeNode<CopyStoragePolicy, 
                            typename GFS::template Child<0>::Type::LocalFunctionSpace::Traits::NodeType,
                            typename GFS::template Child<1>::Type::LocalFunctionSpace::Traits::NodeType,
                            typename GFS::template Child<2>::Type::LocalFunctionSpace::Traits::NodeType,
                            typename GFS::template Child<3>::Type::LocalFunctionSpace::Traits::NodeType> Type;
    };

    template<typename GFS>
    struct CompositeLocalFunctionSpaceNodeBaseTypeTraits<GFS,5>
    {
      typedef CompositeNode<CopyStoragePolicy, 
                            typename GFS::template Child<0>::Type::LocalFunctionSpace::Traits::NodeType,
                            typename GFS::template Child<1>::Type::LocalFunctionSpace::Traits::NodeType,
                            typename GFS::template Child<2>::Type::LocalFunctionSpace::Traits::NodeType,
                            typename GFS::template Child<3>::Type::LocalFunctionSpace::Traits::NodeType,
                            typename GFS::template Child<4>::Type::LocalFunctionSpace::Traits::NodeType> Type;
    };

    template<typename GFS>
    struct CompositeLocalFunctionSpaceNodeBaseTypeTraits<GFS,6>
    {
      typedef CompositeNode<CopyStoragePolicy, 
                            typename GFS::template Child<0>::Type::LocalFunctionSpace::Traits::NodeType,
                            typename GFS::template Child<1>::Type::LocalFunctionSpace::Traits::NodeType,
                            typename GFS::template Child<2>::Type::LocalFunctionSpace::Traits::NodeType,
                            typename GFS::template Child<3>::Type::LocalFunctionSpace::Traits::NodeType,
                            typename GFS::template Child<4>::Type::LocalFunctionSpace::Traits::NodeType,
                            typename GFS::template Child<5>::Type::LocalFunctionSpace::Traits::NodeType> Type;
    };

    template<typename GFS>
    struct CompositeLocalFunctionSpaceNodeBaseTypeTraits<GFS,7>
    {
      typedef CompositeNode<CopyStoragePolicy, 
                            typename GFS::template Child<0>::Type::LocalFunctionSpace::Traits::NodeType,
                            typename GFS::template Child<1>::Type::LocalFunctionSpace::Traits::NodeType,
                            typename GFS::template Child<2>::Type::LocalFunctionSpace::Traits::NodeType,
                            typename GFS::template Child<3>::Type::LocalFunctionSpace::Traits::NodeType,
                            typename GFS::template Child<4>::Type::LocalFunctionSpace::Traits::NodeType,
                            typename GFS::template Child<5>::Type::LocalFunctionSpace::Traits::NodeType,
                            typename GFS::template Child<6>::Type::LocalFunctionSpace::Traits::NodeType> Type;
    };

    template<typename GFS>
    struct CompositeLocalFunctionSpaceNodeBaseTypeTraits<GFS,8>
    {
      typedef CompositeNode<CopyStoragePolicy, 
                            typename GFS::template Child<0>::Type::LocalFunctionSpace::Traits::NodeType,
                            typename GFS::template Child<1>::Type::LocalFunctionSpace::Traits::NodeType,
                            typename GFS::template Child<2>::Type::LocalFunctionSpace::Traits::NodeType,
                            typename GFS::template Child<3>::Type::LocalFunctionSpace::Traits::NodeType,
                            typename GFS::template Child<4>::Type::LocalFunctionSpace::Traits::NodeType,
                            typename GFS::template Child<5>::Type::LocalFunctionSpace::Traits::NodeType,
                            typename GFS::template Child<6>::Type::LocalFunctionSpace::Traits::NodeType,
                            typename GFS::template Child<7>::Type::LocalFunctionSpace::Traits::NodeType> Type;
    };

    template<typename GFS>
    struct CompositeLocalFunctionSpaceNodeBaseTypeTraits<GFS,9>
    {
      typedef CompositeNode<CopyStoragePolicy, 
                            typename GFS::template Child<0>::Type::LocalFunctionSpace::Traits::NodeType,
                            typename GFS::template Child<1>::Type::LocalFunctionSpace::Traits::NodeType,
                            typename GFS::template Child<2>::Type::LocalFunctionSpace::Traits::NodeType,
                            typename GFS::template Child<3>::Type::LocalFunctionSpace::Traits::NodeType,
                            typename GFS::template Child<4>::Type::LocalFunctionSpace::Traits::NodeType,
                            typename GFS::template Child<5>::Type::LocalFunctionSpace::Traits::NodeType,
                            typename GFS::template Child<6>::Type::LocalFunctionSpace::Traits::NodeType,
                            typename GFS::template Child<7>::Type::LocalFunctionSpace::Traits::NodeType,
                            typename GFS::template Child<8>::Type::LocalFunctionSpace::Traits::NodeType> Type;
     
    };
#endif

    // local function space for a power grid function space
    template<typename GFS> 
    class CompositeLocalFunctionSpaceNode :
      public LocalFunctionSpaceBaseNode<GFS>,
      public CompositeLocalFunctionSpaceNodeBaseTypeTraits<GFS,GFS::CHILDREN>::Type
    {
      typedef LocalFunctionSpaceBaseNode<GFS> BaseT;

      // friend decl for bind meta program
      template<typename T, bool b, typename E, typename GC, typename Int> 
      friend struct LocalFunctionSpaceBaseVisitNodeMetaProgram;
      template<typename T, typename E, typename GC, typename Int, int n, int i>
      friend struct LocalFunctionSpaceBaseVisitChildMetaProgram;

    public:
      typedef PowerCompositeLocalFunctionSpaceTraits<GFS,CompositeLocalFunctionSpaceNode> Traits;

      //! \brief empty constructor (needed for CopyStoragePolicy)
      CompositeLocalFunctionSpaceNode ()
      {
      }

      //! \brief initialize with grid function space
      CompositeLocalFunctionSpaceNode (const GFS& gfs)  : BaseT(gfs)
      {
        setup(gfs);
      }

      //! \brief initialize with grid function space
      void setup (const GFS& gfs)
      {
        BaseT::setup(gfs);
        CompositeLocalFunctionSpaceNodeVisitChildMetaProgram<CompositeLocalFunctionSpaceNode,GFS::CHILDREN,0>::
          setup(*this,gfs);
      }

      //! \brief bind local function space to entity
      void bind (const typename Traits::Element& e)
      {
        // call method on base class, this avoid the barton neckman trick
        BaseT::bind(*this,e);
      }

    };

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
    class LeafLocalFunctionSpaceNode :
      public LocalFunctionSpaceBaseNode<GFS>,
      public LeafNode
    {
      typedef LocalFunctionSpaceBaseNode<GFS> BaseT;

      // friend decl for bind meta program
      template<typename T, bool b, typename E, typename GC, typename Int> 
      friend struct LocalFunctionSpaceBaseVisitNodeMetaProgram;
      template<typename T, typename E, typename GC, typename Int, int n, int i>
      friend struct LocalFunctionSpaceBaseVisitChildMetaProgram;

    public:
      typedef LeafLocalFunctionSpaceTraits<GFS,LeafLocalFunctionSpaceNode> Traits;

    private:
      typedef FiniteElementInterfaceSwitch<
      typename Traits::FiniteElementType
      > FESwitch;

    public:
      //! \brief empty constructor (needed for CopyStoragePolicy)
      LeafLocalFunctionSpaceNode ()
      {
      }

      //! \brief initialize with grid function space
      LeafLocalFunctionSpaceNode (const GFS& gfs) : BaseT(gfs)
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
    public:
      typedef typename BaseT::Traits Traits;

      LocalFunctionSpace(const GFS & gfs) : BaseT(gfs) {}

      LocalIndex localIndex (typename Traits::IndexContainer::size_type index) const
      {
        return LocalIndex(BaseT::localIndex(index));
      }

    private:
      // we don't support getChild yet, so let's hide it!
      template<int i>
      void getChild () const;
    };

    // specialization for AnySpaceTag
    template <typename GFS>
    class LocalFunctionSpace<GFS, AnySpaceTag> : 
      public GFS::LocalFunctionSpace
    {
      typedef typename GFS::LocalFunctionSpace BaseT;
    public:
      LocalFunctionSpace(const GFS & gfs) : BaseT(gfs) {}
    };

    //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif
