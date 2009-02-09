// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_LOCALFUNCTIONSPACE_HH
#define DUNE_PDELAB_LOCALFUNCTIONSPACE_HH

#include<vector>

#include"../common/multitypetree.hh"
#include"../common/cpstoragepolicy.hh"

namespace Dune {
  namespace PDELab {

    //! \addtogroup GridFunctionSpace
    //! \ingroup PDELab
    //! \{

    //=======================================
    // local function space base: metaprograms
    //=======================================

	template<typename T, bool isleaf, typename E, typename It, typename Int>
	struct LocalFunctionSpaceBaseVisitNodeMetaProgram;

	template<typename T, typename E, typename It, typename Int, int n, int i>
	struct LocalFunctionSpaceBaseVisitChildMetaProgram // visit i'th child of inner node
	{
	  static void bind_localfunctionspace_to_element (T& t, const E& e, It begin, Int& offset)
	  {
        // vist children of node t in order
		typedef typename T::template Child<i>::Type C;
        Int initial_offset = offset; // remember initial offset to compute size later
        LocalFunctionSpaceBaseVisitNodeMetaProgram<C,C::isLeaf,E,It,Int>::
          bind_localfunctionspace_to_element(t.template getChild<i>(),e,begin,offset);
        for (int j=initial_offset; j<offset; j++)
          begin[j] = t.pgfs->template subMap<i>(begin[j]);
        LocalFunctionSpaceBaseVisitChildMetaProgram<T,E,It,Int,n,i+1>::
          bind_localfunctionspace_to_element(t,e,begin,offset);
	  }
	};

	template<typename T, typename E, typename It, typename Int, int n>
	struct LocalFunctionSpaceBaseVisitChildMetaProgram<T,E,It,Int,n,n> // end of child recursion
	{
	  static void bind_localfunctionspace_to_element (T& t, const E& e, It begin, Int& offset)
	  {
        return;
	  }
	};

	template<typename T, bool isleaf, typename E, typename It, typename Int> 
	struct LocalFunctionSpaceBaseVisitNodeMetaProgram // visit inner node
	{
	  static void bind_localfunctionspace_to_element (T& t, const E& e, It begin, Int& offset)
	  {
        // now we are at a multi component local function space
        Int initial_offset = offset; // remember initial offset to compute size later
        t.i = begin+initial_offset; // begin is always the first entry in the vector
		LocalFunctionSpaceBaseVisitChildMetaProgram<T,E,It,Int,T::CHILDREN,0>::
          bind_localfunctionspace_to_element(t,e,begin,offset);
        t.n = offset-initial_offset;
	  }
	};

	template<typename T, typename E, typename It, typename Int> 
	struct LocalFunctionSpaceBaseVisitNodeMetaProgram<T,true,E,It,Int> // visit leaf node 
	{
	  static void bind_localfunctionspace_to_element (T& t, const E& e, It begin, Int& offset)
	  {
        // now we are at a single component local function space
        // which is part of a multi component local function space
        t.i = begin+offset; // begin is always the first entry in the vector
        t.plfem = &(((t.pgfs)->localFiniteElementMap()).find(e));
        t.n = t.plfem->localBasis().size(); // determine size of this chunk
        std::vector<typename T::Traits::GridFunctionSpaceType::Traits::SizeType> global(t.n);
        t.pgfs->globalIndices(*(t.plfem),e,global); // get global indices for this finite element
        for (Int i=0; i<t.n; i++) t.i[i]=t.pgfs->upMap(global[i]); 
        offset += t.n; // append this chunk
	  }
	};

    //=======================================
    // local function space base: power implementation
    //=======================================

	//! traits for multi component local function space
	template<typename GFS, typename N>	
    struct PowerCompositeLocalFunctionSpaceTraits
    {
	  //! \brief the grid view where grid function is defined upon
	  typedef GFS GridFunctionSpaceType;

      //! type of local function space node
      typedef N NodeType;

	  //! \brief Type to store indices from Backend
      typedef typename GFS::Traits::GridViewType GridViewType;

      //! \brief Type of codim 0 entity in the grid
      typedef typename GridViewType::Traits::template Codim<0>::Entity Element;

	  //! \brief Type to store indices from Backend
      typedef typename GFS::Traits::SizeType SizeType;

	  //! \brief Type of container to store indices
      typedef typename std::vector<SizeType> IndexContainer;
    };


    // local function space for a power grid function space
    template<typename GFS> 
    class PowerLocalFunctionSpaceNode 
      : public PowerNode<typename GFS::template Child<0>::Type::LocalFunctionSpace::Traits::NodeType,
                         GFS::CHILDREN,CopyStoragePolicy>
    {
      template<typename T, bool b, typename E, typename It, typename Int> 
      friend struct LocalFunctionSpaceBaseVisitNodeMetaProgram;
      template<typename T, typename E, typename It, typename Int, int n, int i>
      friend struct LocalFunctionSpaceBaseVisitChildMetaProgram;

      typedef typename GFS::Traits::BackendType B;
 	  typedef typename GFS::Traits::GridViewType::Traits::template Codim<0>::Entity Element;

      typedef typename GFS::template Child<0>::Type::LocalFunctionSpace::Traits::NodeType NodeType;
      typedef PowerNode<NodeType,GFS::CHILDREN,CopyStoragePolicy> BaseT;

   public:
      typedef PowerCompositeLocalFunctionSpaceTraits<GFS,PowerLocalFunctionSpaceNode> Traits;

      //! \brief empty constructor
      PowerLocalFunctionSpaceNode ()
      {
      }

      //! \brief initialize with grid function space
      PowerLocalFunctionSpaceNode (const GFS& gfs)  : pgfs(&gfs)
      {
        setup(gfs);
      }

      //! \brief initialize with grid function space
      void setup (const GFS& gfs) const
      {
        for (int i=0; i<GFS::CHILDREN; i++)
          this->getChild(i).setup(pgfs->getChild(i));
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

      /** \brief extract coefficients for one element from container */  
      template<typename GC, typename LC>
      void vread (const GC& globalcontainer, LC& localcontainer) const
      {
        localcontainer.resize(n);
        for (typename Traits::IndexContainer::size_type k=0; k<n; ++k)
          localcontainer[k] = B::const_access(globalcontainer,i[k]);
      }

      /** \brief write back coefficients for one element to container */  
      template<typename GC, typename LC>
      void vwrite (const LC& localcontainer, GC& globalcontainer) const
      {
        for (typename Traits::IndexContainer::size_type k=0; k<n; ++k)
          B::access(globalcontainer,i[k]) = localcontainer[k];
      }

      /** \brief add coefficients for one element to container */  
      template<typename GC, typename LC>
      void vadd (const LC& localcontainer, GC& globalcontainer) const
      {
        for (typename Traits::IndexContainer::size_type k=0; k<n; ++k)
          B::access(globalcontainer,i[k]) += localcontainer[k];
      }

      void debug () const
      {
        std::cout << n << " indices = (";
        for (typename Traits::IndexContainer::size_type k=0; k<n; k++)
          std::cout << i[k] << " ";
        std::cout << ")" << std::endl;
      }

    private:
	  CP<GFS const> pgfs;
      typename Traits::IndexContainer::iterator i;
      typename Traits::IndexContainer::size_type n;
     };


    // local function space description that can be bound to an element
    // depends on a grid function space
    template<typename GFS>
    class PowerLocalFunctionSpace : public PowerLocalFunctionSpaceNode<GFS>
    {
      typedef PowerLocalFunctionSpaceNode<GFS> BaseT;

    public:
      typedef typename BaseT::Traits Traits;

      PowerLocalFunctionSpace (const GFS& gfs) 
        : BaseT(gfs), global(gfs.maxLocalSize())
      {}

      //! \brief bind local function space to entity
      void bind (const typename Traits::Element& e)
      {
        typename Traits::IndexContainer::size_type offset=0;

        // is implemented as a template metaprogram over BaseT
        LocalFunctionSpaceBaseVisitNodeMetaProgram<BaseT,BaseT::isLeaf,
          typename Traits::Element,
          typename Traits::IndexContainer::iterator,
          typename Traits::IndexContainer::size_type>::
          bind_localfunctionspace_to_element(*this,e,global.begin(),offset);
        global.resize(offset); // now the size is known
      }

    private:
      typename BaseT::Traits::IndexContainer global;
    };


    //=======================================
    // local function space base: single component implementation
    //=======================================

	//! traits for single component local function space
	template<typename GFS, typename N>	
	struct LocalFunctionSpaceTraits : public PowerCompositeLocalFunctionSpaceTraits<GFS,N>
	{
	  //! \brief Type of local finite element
      typedef typename GFS::Traits::LocalFiniteElementType LocalFiniteElementType;
	};

    template<typename GFS> 
    class LocalFunctionSpaceNode : public LeafNode
    {
      template<typename T, bool b, typename E, typename It, typename Int> 
      friend struct LocalFunctionSpaceBaseVisitNodeMetaProgram;
      template<typename T, typename E, typename It, typename Int, int n, int i>
      friend struct LocalFunctionSpaceBaseVisitChildMetaProgram;

      typedef typename GFS::Traits::BackendType B;
 	  typedef typename GFS::Traits::GridViewType::Traits::template Codim<0>::Entity Element;

   public:
      typedef LocalFunctionSpaceTraits<GFS,LocalFunctionSpaceNode> Traits;

      //! \brief empty constructor
      LocalFunctionSpaceNode ()
      {
      }

      //! \brief initialize with grid function space
      LocalFunctionSpaceNode (const GFS& gfs) : pgfs(&gfs)
      {
      }

      //! \brief initialize with grid function space
      void setup (const GFS& gfs) const
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

      //! \brief get local finite element
      const typename Traits::LocalFiniteElementType& localFiniteElement () const
      {
        return *plfem;
      }

      /** \brief extract coefficients for one element from container */  
      template<typename GC, typename LC>
      void vread (const GC& globalcontainer, LC& localcontainer) const
      {
        localcontainer.resize(n);
        for (typename Traits::IndexContainer::size_type k=0; k<n; ++k)
          localcontainer[k] = B::const_access(globalcontainer,i[k]);
      }

      /** \brief write back coefficients for one element to container */  
      template<typename GC, typename LC>
      void vwrite (const LC& localcontainer, GC& globalcontainer) const
      {
        for (typename Traits::IndexContainer::size_type k=0; k<n; ++k)
          B::access(globalcontainer,i[k]) = localcontainer[k];
      }

      /** \brief add coefficients for one element to container */  
      template<typename GC, typename LC>
      void vadd (const LC& localcontainer, GC& globalcontainer) const
      {
        for (typename Traits::IndexContainer::size_type k=0; k<n; ++k)
          B::access(globalcontainer,i[k]) += localcontainer[k];
      }

      void debug () const
      {
        std::cout << n << " indices = (";
        for (typename Traits::IndexContainer::size_type k=0; k<n; k++)
          std::cout << i[k] << " ";
        std::cout << ")" << std::endl;
      }

    private:
	  mutable CP<GFS const> pgfs;
      typename Traits::IndexContainer::iterator i;
      typename Traits::IndexContainer::size_type n;
      const typename Traits::LocalFiniteElementType* plfem;
    };


    // local function space description that can be bound to an element
    // depends on a grid function space
    template<typename GFS>
    class LocalFunctionSpace : public LocalFunctionSpaceNode<GFS>
    {
      typedef LocalFunctionSpaceNode<GFS> BaseT;

    public:
      typedef typename BaseT::Traits Traits;

      LocalFunctionSpace (const GFS& gfs) 
        : BaseT(gfs), global(gfs.maxLocalSize())
      {}

      //! \brief bind local function space to entity
      void bind (const typename Traits::Element& e)
      {
        typename Traits::IndexContainer::size_type offset=0;

        // is implemented as a template metaprogram over BaseT
        LocalFunctionSpaceBaseVisitNodeMetaProgram<BaseT,BaseT::isLeaf,
          typename Traits::Element,
          typename Traits::IndexContainer::iterator,
          typename Traits::IndexContainer::size_type>::
          bind_localfunctionspace_to_element(*this,e,global.begin(),offset);
        global.resize(offset); // now the size is known
      }

    private:
      typename BaseT::Traits::IndexContainer global;
    };

    //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif
