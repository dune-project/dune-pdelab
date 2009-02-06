// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_LEAFGRIDFUNCTIONSPACE_HH
#define DUNE_PDELAB_LEAFGRIDFUNCTIONSPACE_HH

#include<vector>
#include<set>
#include<map>

#include<dune/common/exceptions.hh>
#include<dune/grid/common/referenceelements.hh>

#include"../common/countingptr.hh"
#include"../common/multitypetree.hh"
#include"../common/cpstoragepolicy.hh"

namespace Dune {
  namespace PDELab {

    //=======================================
    // grid function space
    //=======================================

	//! collect types exported by a leaf grid function space
	template<typename G, class L, typename B>
	struct GridFunctionSpaceTraits
	{
	  //! \brief the grid view where grid function is defined upon
	  typedef G GridViewType;

	  //! \brief vector backend
	  typedef B BackendType;

	  //! \brief short cut for size type exported by Backend
	  typedef typename B::size_type SizeType;

	  //! \brief local finite element
	  typedef L LocalFiniteElementMapType;

	  //! \brief local finite element
	  typedef typename L::Traits::LocalFiniteElementType LocalFiniteElementType;

	  //! \brief type for constraints coefficients taken as range type from finite element
	  typedef typename L::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType CoefficientType;

	  //! \brief type to represent one row of constraints
	  typedef std::map<SizeType,CoefficientType> RowType;

	  //! \brief type to store constraints
	  typedef std::map<SizeType,RowType> TransformationType;
	};


	//! Simple Backend for std::vector
	class StdVectorBackend
	{
	public:
	  //! container construction
	  template<typename T, typename E>
	  class VectorContainer : public std::vector<E>
	  {
		typedef std::vector<E> BaseT;
	  public:
		typedef E ElementType;

		VectorContainer (const T& t) : BaseT(t.globalSize()) {}
		VectorContainer (const T& t, const E& e) : BaseT(t.globalSize(),e) {}
		VectorContainer& operator= (const E& e) // set all elements to same value
		{
		  for (typename BaseT::size_type i=0; i<BaseT::size(); i++)
			(*this)[i] = e;
		  return *this;
		}
	  };

	  //! extract type of container element 
	  template<class C>
	  struct Value
	  {
        //! type of a container element
		typedef typename C::value_type Type;
	  };

	  //! The size type
	  typedef std::vector<int>::size_type size_type;

	  // get const_reference to container element
	  // we can assume C to be std::vector<T>
	  template<typename C>
	  static const typename C::value_type& const_access (const C& c, size_type i)
	  {
		return c.operator[](i);
	  }

	  // get non const_reference to container element 
	  // note: this method does not depend on T!
	  template<typename C>
	  static typename C::value_type& access (C& c, size_type i)
	  {
		return c.operator[](i);
	  }
	};


    // forward declaration of local function space
    template<typename GFS> class LocalFunctionSpace;


	// template metaprogram to evaluate subindex without compile-time codim parameter
	// C : codimension
	template<int C, typename IS, typename E>
	class EvalSubIndexMetaProgram
	{
	public:
	  static void subIndex (const IS& is, const E& e, unsigned int i, 
							unsigned int c, unsigned int& index)
	  {
		if (c==C)
		  {
			index = is.template subIndex<C>(e,i);
			return;
		  }
		EvalSubIndexMetaProgram<C-1,IS,E>::subIndex(is,e,i,c,index);	
	  }
	};

	template<typename IS, typename E>
	class EvalSubIndexMetaProgram<0,IS,E>
	{
	public:
	  static void subIndex (const IS& is, const E& e, unsigned int i, 
							unsigned int c, unsigned int& index)
	  {
		if (c==0)
		  index = is.index(e); // there is only one codim 0 entity
	  }
	};

	template<int DIM, typename IS, typename E>
	unsigned int eval_subindex (const IS& is, const E& e, unsigned int i, unsigned int c)
	{
	  unsigned int index=0;
	  EvalSubIndexMetaProgram<DIM,IS,E>::subIndex(is,e,i,c,index);
	  return index;
	}

	// this class may be used to pass compile-time
	// parameters to the implementation 
	struct GridFunctionGeneralMapper
	{
	  enum {dummy=0} ;
	};


	// mapper for layouts with arbitrary number of entries per entity
	// GV : Type implementing GridView
	// FEM  : Type implementing LocalFiniteElementMapInterface
	// B : Backend type
	// P : Parameter type
	template<typename GV, typename LFEM, typename B=StdVectorBackend, 
			 typename P=GridFunctionGeneralMapper>
	class GridFunctionSpace : public Countable, public LeafNode
	{
	public:
	  typedef GridFunctionSpaceTraits<GV,LFEM,B> Traits;
	  typedef typename GV::Traits::template Codim<0>::Entity Element;
	  typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;

	  // extract type of container storing Es
	  template<typename T>
	  struct VectorContainer
	  {
		//! \brief define Type as the Type of a container of E's
		typedef typename B::template VectorContainer<GridFunctionSpace,T> Type;	
	  };

      // define local function space parametrized by self 
      typedef LocalFunctionSpace<GridFunctionSpace> LocalFunctionSpace;

	  // constructor
	  GridFunctionSpace (const GV& gridview, const LFEM& lfem) 
		: gv(gridview), plfem(&lfem)
	  {
		update();
	  }

	  // get grid view
	  const GV& gridview () const
	  {
		return gv;
	  }

	  // get finite element map, I think we dont need it
	  const LFEM& localFiniteElementMap () const
	  {
		return *plfem;
	  }

	  // get dimension of finite element space
	  typename Traits::SizeType globalSize () const
	  {
		return nglobal;
	  }

	  // get max dimension of shape function space
	  typename Traits::SizeType maxLocalSize () const
	  {
		return nlocal;
	  }

	  // map index [0,globalSize-1] to root index set
	  typename Traits::SizeType upMap (typename Traits::SizeType i) const
	  {
		return i;
	  }

	  // compute global indices for one element
	  void globalIndices (const typename Traits::LocalFiniteElementType& lfe, 
                          const Element& e, 
						  std::vector<typename Traits::SizeType>& global) const
	  {
		// get layout of entity
		const typename Traits::LocalFiniteElementType::Traits::LocalCoefficientsType&
		  lc = lfe.localCoefficients();
		global.resize(lc.size());

		for (unsigned int i=0; i<lc.size(); ++i)
		  {
			// get geometry type of subentity 
			Dune::GeometryType gt=Dune::ReferenceElements<double,GV::Grid::dimension>
			  ::general(lfe.type()).type(lc.localIndex(i).subentity(),lc.localIndex(i).codim());

			// evaluate consecutive index of subentity
			int index = eval_subindex<GV::Grid::dimension>(gv.indexSet(),e,
														   lc.localIndex(i).subentity(),
														   lc.localIndex(i).codim());
		
			// now compute 
			global[i] = offset[(gtoffset.find(gt)->second)+index]+lc.localIndex(i).index();
		  }
	  }

      // global Indices from element, needs additional finite element lookup
	  void globalIndices (const Element& e,
						  std::vector<typename Traits::SizeType>& global) const
      {
        globalIndices(plfem->find(e),e,global);
      }

	  // update information, e.g. when grid has changed
	  void update ()
	  {
		std::cout << "GridFunctionSpace(general version):" << std::endl;

		// determine which geometry types are used
		// needs one traversal of the grid
		typedef std::set<Dune::GeometryType> GtUsedSetType;
		GtUsedSetType gtused;
		for (ElementIterator it = gv.template begin<0>();
			 it!=gv.template end<0>(); ++it)
		  {
			// check geometry type
			if ((plfem->find(*it)).type()!=it->type())
			  DUNE_THROW(Exception, "geometry type mismatch in GridFunctionSpace");

			// get local coefficients for this entity
			const typename Traits::LocalFiniteElementType::Traits::LocalCoefficientsType&
			  lc = (plfem->find(*it)).localCoefficients();

			// insert geometry type of all subentities into set
			for (int i=0; i<lc.size(); ++i)
			  {
				Dune::GeometryType gt=Dune::ReferenceElements<double,GV::Grid::dimension>
				  ::general(it->type()).type(lc.localIndex(i).subentity(),lc.localIndex(i).codim());
				gtused.insert(gt);
			  }
		  }

		// now we can allocate one number per entity that holds degrees of freedom
		typename Traits::SizeType nentities = 0;
		gtoffset.clear();
		const typename GV::IndexSet& is=gv.indexSet();
		for (typename GtUsedSetType::iterator i=gtused.begin(); i!=gtused.end(); ++i)
		  {
			gtoffset[*i] = nentities;
			std::cout << *i << ": " << is.size(*i) << " entries at "
					  << nentities << std::endl;
			nentities += is.size(*i);
		  }
		offset.resize(nentities);
		for (typename std::vector<typename Traits::SizeType>::iterator i=offset.begin(); i!=offset.end(); ++i)
		  *i = 0;
		std::cout << "allocated offset vector with size " << offset.size() << std::endl;

		// now compute the number of entries for each entity
		// requires second grid traversal
		nlocal = 0;
		for (ElementIterator it = gv.template begin<0>();
			 it!=gv.template end<0>(); ++it)
		  {
			// get local coefficients for this entity
			const typename Traits::LocalFiniteElementType::Traits::LocalCoefficientsType&
			  lc = (plfem->find(*it)).localCoefficients();

			// compute maximum number of degrees of freedom per element
			nlocal = std::max(nlocal,static_cast<typename Traits::SizeType>(lc.size()));

			// compute maximum size for each subentity
			for (int i=0; i<lc.size(); ++i)
			  {
				Dune::GeometryType gt=Dune::ReferenceElements<double,GV::Grid::dimension>
				  ::general(it->type()).type(lc.localIndex(i).subentity(),lc.localIndex(i).codim());
				unsigned int index = gtoffset[gt]
				  +eval_subindex<GV::Grid::dimension>(is,*it,lc.localIndex(i).subentity(),lc.localIndex(i).codim());
				offset[index] = std::max(offset[index],lc.localIndex(i).index()+1);
			  }
		  }

		// now count global number of dofs and compute offset
		nglobal = 0;
		for (typename std::vector<typename Traits::SizeType>::iterator i=offset.begin(); 
			 i!=offset.end(); ++i)
		  {
			typename Traits::SizeType size = *i;
			*i = nglobal;
			nglobal += size;
		  }
		std::cout << "total number of dofs is " << nglobal << std::endl;
	  }

	private:
	  const GV& gv;
	  CP<LFEM const> plfem;
	  typename Traits::SizeType nlocal;
	  typename Traits::SizeType nglobal;

	  std::map<Dune::GeometryType,typename Traits::SizeType> gtoffset; // offset in vector for given geometry type
	  std::vector<typename Traits::SizeType> offset; // offset into big vector for each entity;
	};



	// this class may be used to pass compile-time
	// parameters to the implementation 
	struct GridFunctionRestrictedMapper
	{
	  enum {dummy=1} ;
	};


	// specialization with restricted mapper
	// GV : Type implementing GridView
	// FEM  : Type implementing LocalFiniteElementMapInterface
	// B : Backend type
	template<typename GV, typename LFEM, typename B> 
	class GridFunctionSpace<GV,LFEM,B,GridFunctionRestrictedMapper> : 
	  public Countable, public LeafNode
	{
	public:
	  typedef GridFunctionSpaceTraits<GV,LFEM,B> Traits;
	  typedef typename GV::Traits::template Codim<0>::Entity Element;
	  typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;

	  // extract type of container storing Es
	  template<typename T>
	  struct VectorContainer
	  {
		//! \brief define Type as the Type of a container of E's
		typedef typename B::template VectorContainer<GridFunctionSpace,T> Type;	
	  };

	  // constructor
	  GridFunctionSpace (const GV& gridview, const LFEM& lfem) 
		: gv(gridview), plfem(&lfem)
	  {
		update();
	  }

	  // get grid view
	  const GV& gridview () const
	  {
		return gv;
	  }

	  // get finite element map, I think we dont need it
	  const LFEM& localFiniteElementMap () const
	  {
		return *plfem;
	  }

	  // get dimension of finite element space
	  typename Traits::SizeType globalSize () const
	  {
		return nglobal;
	  }

	  // get max dimension of shape function space
	  typename Traits::SizeType maxLocalSize () const
	  {
		return nlocal;
	  }

	  // map index [0,globalSize-1] to root index set
	  typename Traits::SizeType upMap (typename Traits::SizeType i) const
	  {
		return i;
	  }

	  // compute global indices for one element
	  void globalIndices (const typename Traits::LocalFiniteElementType& lfe,
                          const Element& e,
						  std::vector<typename Traits::SizeType>& global) const
	  {
		// get local coefficients for this entity
		const typename Traits::LocalFiniteElementType::Traits::LocalCoefficientsType&
		  lc = lfe.localCoefficients();
		global.resize(lc.size());

		for (unsigned int i=0; i<lc.size(); ++i)
		  {
			// get geometry type of subentity 
			Dune::GeometryType gt=Dune::ReferenceElements<double,GV::Grid::dimension>
			  ::general(lfe.type()).type(lc.localIndex(i).subentity(),lc.localIndex(i).codim());

			// evaluate consecutive index of subentity
			int index = eval_subindex<GV::Grid::dimension>(gv.indexSet(),e,
														   lc.localIndex(i).subentity(),
														   lc.localIndex(i).codim());
		
			// now compute 
			global[i] = offset.find(gt)->second+index*dofcountmap.find(gt)->second+lc.localIndex(i).index();
		  }
	  }

      // global Indices from element, needs additional finite element lookup
	  void globalIndices (const Element& e,
						  std::vector<typename Traits::SizeType>& global) const
      {
        globalIndices(plfem->find(e),e,global);
      }

	  // update information, e.g. when grid has changed
	  void update ()
	  {
		std::cout << "GridFunctionSpace(restricted version):" << std::endl;

		// clear counters
		dofcountmap.clear();
		nlocal = 0;

		// count number of dofs in each subentity
		for (ElementIterator it = gv.template begin<0>();
			 it!=gv.template end<0>(); ++it)
		  {
			// check geometry type
			if ((plfem->find(*it)).type()!=it->type())
			  DUNE_THROW(Exception, "geometry type mismatch in GridFunctionSpace");

			// get local coefficients for this entity
			const typename Traits::LocalFiniteElementType::Traits::LocalCoefficientsType&
			  lc = (plfem->find(*it)).localCoefficients();

			// compute maximum number of degrees of freedom per element
			nlocal = std::max(nlocal,static_cast<typename Traits::SizeType>(lc.size()));

			// store count for each subentity in a map
			typedef Dune::tuple<unsigned int, unsigned int> SubentityType;
			typedef std::map<SubentityType,unsigned int> CountMapType;
			CountMapType countmap;

			// assume that key within each subentity is unique
			for (unsigned int i=0; i<lc.size(); ++i)
			  {
				SubentityType subentity(lc.localIndex(i).subentity(),lc.localIndex(i).codim());
				if (countmap.find(subentity)==countmap.end())
				  countmap[subentity] = 1;
				else
				  (countmap[subentity])++;
			  }

			// traverse the map and compare #dofs per geometry type
			for (typename CountMapType::iterator i=countmap.begin(); i!=countmap.end(); ++i)
			  {
				Dune::GeometryType gt=Dune::ReferenceElements<double,GV::Grid::dimension>
				  ::general(it->type()).type(Dune::get<0>(i->first),Dune::get<1>(i->first));
				typename DofCountMapType::iterator j=dofcountmap.find(gt);
				if (j==dofcountmap.end())
				  {
					dofcountmap[gt] = i->second;
					continue;
				  }
				if (j->second != i->second)
				  {
					DUNE_THROW(Dune::NotImplemented, "non constant # dofs per geometry type, use general version instead");
				  }
			  }
		  }
		std::cout << "max local number of dofs = " << nlocal << std::endl;

		// print result
		for (typename DofCountMapType::iterator i=dofcountmap.begin(); i!=dofcountmap.end(); ++i)
		  {
			std::cout << i->first << " has " << i->second << " degrees of freedom" << std::endl;
		  }

		// compute offsets
		nglobal = 0;
		offset.clear();
		const typename GV::IndexSet& is=gv.indexSet();
		for (typename DofCountMapType::iterator i=dofcountmap.begin(); i!=dofcountmap.end(); ++i)
		  {
			offset[i->first] = nglobal;
			nglobal += is.size(i->first)*(i->second); 
			std::cout << i->first << " offset now " << nglobal << std::endl;
		  }
		std::cout << "total number of dofs = " << nglobal << std::endl;
	  }

	private:
	  const GV& gv;
	  CP<LFEM const> plfem;

	  typename Traits::SizeType nlocal;
	  typename Traits::SizeType nglobal;

	  typedef std::map<Dune::GeometryType,typename Traits::SizeType> DofCountMapType;
	  DofCountMapType dofcountmap; // number of degrees of freedom per geometry type
	  std::map<Dune::GeometryType,typename Traits::SizeType> offset; // offset in vector for given geometry type
	};


    //=======================================
    // local function space
    //=======================================

	template<typename T, bool isleaf, typename E, typename It, typename Int>
	struct LocalFunctionSpaceBaseVisitNodeMetaProgram;

	template<typename T, typename E, typename It, typename Int, int n, int i>
	struct LocalFunctionSpaceBaseVisitChildMetaProgram // visit child of inner node
	{
	  static void bind_localfunctionspace_to_element (T& t, const E& e, It begin, Int& offset)
	  {
        // vist children of node t in order
		typedef typename T::template Child<i>::Type C;
        LocalFunctionSpaceBaseVisitNodeMetaProgram<C,C::isLeaf,E,It,Int>::
          bind_localfunctionspace_to_element(t.template getChild<i>(),e,begin,offset);
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


    // implementation of local function space needs
    // to be specialized for interior and leaves
    template<typename GFS, bool isleaf> class LocalFunctionSpaceBase;

    // local function space description that can be bound to an element
    // depends on a grid function space
    template<typename GFS>
    class LocalFunctionSpace : public LocalFunctionSpaceBase<GFS,GFS::isLeaf>
    {
      typedef LocalFunctionSpaceBase<GFS,GFS::isLeaf> BaseT;

    public:
      typedef typename BaseT::Traits Traits;

      LocalFunctionSpace (const GFS& gfs) 
        : BaseT(gfs), global(gfs.maxLocalSize())
      {}

      //! \brief bind local function space to entity
      void bind (const typename Traits::Element& e)
      {
        typename BaseT::Traits::IndexContainer::size_type offset=0;

        // is implemented as a template metaprogram over BaseT
        LocalFunctionSpaceBaseVisitNodeMetaProgram<BaseT,BaseT::isLeaf,
          typename Traits::Element,
          typename BaseT::Traits::IndexContainer::iterator,
          typename BaseT::Traits::IndexContainer::size_type>::
          bind_localfunctionspace_to_element(*this,e,global.begin(),offset);
        global.resize(offset); // now the size is known
      }

    private:
      typename BaseT::Traits::IndexContainer global;
    };


	//! traits for multi component local function space
	template<typename GFS, bool isleaf>	
    struct LocalFunctionSpaceTraits
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

    // implementation for multi component local function space
    template<typename GFS, bool isleaf> 
    class LocalFunctionSpaceBase
    {
      template<typename T, bool b, typename E, typename It, typename Int> 
      friend struct LocalFunctionSpaceBaseVisitNodeMetaProgram;
      template<typename T, typename E, typename It, typename Int, int n, int i>
      friend struct LocalFunctionSpaceBaseVisitChildMetaProgram;

      typedef typename GFS::Traits::BackendType B;
 	  typedef typename GFS::Traits::GridViewType::Traits::template Codim<0>::Entity Element;

    public:
      typedef LocalFunctionSpaceTraits<GFS,false> Traits;

      //! \brief initialize with grid function space
      LocalFunctionSpaceBase (const GFS& gfs) : pgfs(&gfs)
      {
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
        for (int k=0; k<n; k++)
          std::cout << i[k] << " ";
        std::cout << ")" << std::endl;
      }

    private:
	  CP<GFS const> pgfs;
      typename Traits::IndexContainer::iterator i;
      typename Traits::IndexContainer::size_type n;
    };

	//! traits for single component local function space
	template<typename GFS>	
	struct LocalFunctionSpaceTraits<GFS,true> : public LocalFunctionSpaceTraits<GFS,false>
	{
	  //! \brief Type of local finite element
      typedef typename GFS::Traits::LocalFiniteElementType LocalFiniteElementType;
	};

    // implementation for leaf nodes
    template<typename GFS> 
    class LocalFunctionSpaceBase<GFS,true>
      : public LeafNode
    {
      template<typename T, bool b, typename E, typename It, typename Int> 
      friend struct LocalFunctionSpaceBaseVisitNodeMetaProgram;
      template<typename T, typename E, typename It, typename Int, int n, int i>
      friend struct LocalFunctionSpaceBaseVisitChildMetaProgram;

      typedef typename GFS::Traits::BackendType B;
 	  typedef typename GFS::Traits::GridViewType::Traits::template Codim<0>::Entity Element;

   public:
      typedef LocalFunctionSpaceTraits<GFS,true> Traits;

      //! \brief initialize with grid function space
      LocalFunctionSpaceBase (const GFS& gfs) : pgfs(&gfs)
      {
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
        for (int k=0; k<n; k++)
          std::cout << i[k] << " ";
        std::cout << ")" << std::endl;
      }

    private:
	  CP<GFS const> pgfs;
      typename Traits::IndexContainer::iterator i;
      typename Traits::IndexContainer::size_type n;
      const typename Traits::LocalFiniteElementType* plfem;
    };

  }
}

#endif
