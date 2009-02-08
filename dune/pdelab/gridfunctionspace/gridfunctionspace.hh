// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_HH

#include<vector>
#include<set>
#include<map>

#include<dune/common/exceptions.hh>
#include<dune/grid/common/referenceelements.hh>

#include"../common/countingptr.hh"
#include"../common/multitypetree.hh"
#include"../common/cpstoragepolicy.hh"

#include"localfunctionspace.hh"

namespace Dune {
  namespace PDELab {

    //! \addtogroup GridFunctionSpace
    //! \ingroup PDELab
    //! \{

    //=======================================
    // grid function space : single component case
    //=======================================

	//! collect types exported by a multi-component grid function space
	template<typename G, typename B>
	struct PowerCompositeGridFunctionSpaceTraits
	{
	  //! \brief the grid view where grid function is defined upon
	  typedef G GridViewType;

	  //! \brief vector backend
	  typedef B BackendType;

	  //! \brief short cut for size type exported by Backend
	  typedef typename B::size_type SizeType;
	};

	//! collect types exported by a leaf grid function space
	template<typename G, class L, typename B>
	struct GridFunctionSpaceTraits : public PowerCompositeGridFunctionSpaceTraits<G,B>
	{
	  //! \brief local finite element
	  typedef L LocalFiniteElementMapType;

	  //! \brief local finite element
	  typedef typename L::Traits::LocalFiniteElementType LocalFiniteElementType;

	  //! \brief type for constraints coefficients taken as range type from finite element
	  typedef typename L::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType CoefficientType;

	  //! \brief type to represent one row of constraints
	  typedef std::map<typename B::size_type,CoefficientType> RowType;

	  //! \brief type to store constraints
	  typedef std::map<typename B::size_type,RowType> TransformationType;
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

	  /** \brief get const_reference to container element
       *
       *  we can assume C to be std::vector<T>
       */
	  template<typename C>
	  static const typename C::value_type& const_access (const C& c, size_type i)
	  {
		return c.operator[](i);
	  }

	  /** \brief get non const_reference to container element 
       *
       *  note: this method does not depend on T!
       */
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
      typedef Dune::PDELab::LocalFunctionSpace<GridFunctionSpace> LocalFunctionSpace;

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

		for (int i=0; i<lc.size(); ++i)
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
			for (int i=0; i<lc.size(); ++i)
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
    // power grid function space
    //=======================================


	// this class may be used to pass compile-time
	// parameters to the implementation 
	struct PowerGridFunctionSpaceLexicographicMapper
	{
	  enum {dummy=0} ;
	};

	// this class may be used to pass compile-time
	// parameters to the implementation 
	struct PowerGridFunctionSpaceBlockwiseMapper
	{
	  enum {dummy=0} ;
	};

	template<typename T, int n, int i>
	struct PowerGridFunctionSpaceBaseVisitChildMetaProgram // visit child of inner node
	{
	};

	template<typename T, int n>
	struct PowerGridFunctionSpaceBaseVisitChildMetaProgram<T,n,n> // end of child recursion
	{
	  static void bind_localfunctionspace_to_element ()
	  {
	  }
	};

	template<typename T, int k, typename P>
	class PowerGridFunctionSpace;

    // product of identical grid function spaces
    // base class that holds implementation of the methods
    // this is the default version with lexicographic ordering
	template<typename T, int k, typename P>
	class PowerGridFunctionSpaceBase 
      : public PowerNode<T,k,CountingPointerStoragePolicy>,
        public Countable
	{
      friend class PowerGridFunctionSpace<T,k,P>;

	public:
      typedef PowerCompositeGridFunctionSpaceTraits<typename T::Traits::GridViewType, 
                                                    typename T::Traits::BackendType>
      Traits;

      // it is part of the trick to have a constructor without arguments
      // setting of the children is then done by the constructors
      // of the specialized derived classes
      PowerGridFunctionSpaceBase ()
      {
      }

	  // get grid view
	  const typename Traits::GridViewType& gridview () const
	  {
		return this->template getChild<0>().gridview();
	  }

	  // get dimension of finite element space
	  typename Traits::SizeType globalSize () const
	  {
        // this is bullshit all children may have different
        // size although they have the same type ...
		return offset[k];
	  }

	  // get max dimension of shape function space
	  typename Traits::SizeType maxLocalSize () const
	  {
        // this is bullshit !
		return maxlocalsize;
	  }

	  // map index [0,globalSize-1] to root index set
	  typename Traits::SizeType upMap (typename Traits::SizeType i) const
	  {
		return i;
	  }

	  // map index [0,globalSize-1] to root index set
      template<int i>
	  typename Traits::SizeType subMap (typename Traits::SizeType j) const
	  {
		return offset[i]+j;
	  }

      // recalculate sizes
      void update ()
      {
         for (int i=0; i<k; i++)
           (*this)[i].update();
         setup();
      }

    private:
      void setup ()
      {
        std::cout << "power grid function space(lexicographic version):" << std::endl;
        std::cout << "( ";
        offset[0] = 0;
        maxlocalsize = 0;
        for (int i=0; i<k; i++)
          {
            childSize[i] = (*this)[i].globalSize();
            std::cout << childSize[i] << " ";
            offset[i+1] = offset[i]+childSize[i];
            maxlocalsize += (*this)[i].maxLocalSize();
          }
        std::cout << ") total size = " << offset[k]
                  << " max local size = " << maxlocalsize 
                  << std::endl;
      }

      typename Traits::SizeType childSize[k];
      typename Traits::SizeType offset[k+1];
      typename Traits::SizeType maxlocalsize;
	};


    // product of identical grid function spaces
    // the specializations of this class just set the members
    // all the methods are generic in the implementation
	template<typename T, int k, typename P=PowerGridFunctionSpaceLexicographicMapper>
	class PowerGridFunctionSpace 
      : public PowerGridFunctionSpaceBase<T,k,P>
 	{
	public:
      typedef PowerCompositeGridFunctionSpaceTraits<typename T::Traits::GridViewType, 
                                                    typename T::Traits::BackendType>
      Traits;

      //! Construct a PowerGridFunction with k clones of the function t
	  PowerGridFunctionSpace (T& t)
      {
        for (int i=0; i<k; i++)
          (*this)[i] = t;
        PowerGridFunctionSpaceBase<T,k,P>::setup();
      }

	  PowerGridFunctionSpace (T** t)  
      {
         for (int i=0; i<k; i++)
           (*this)[i] = *(t[i]);
        PowerGridFunctionSpaceBase<T,k,P>::setup();
     }

	};

	template<typename T, typename P>
	class PowerGridFunctionSpace<T,2,P> 
      : public PowerGridFunctionSpaceBase<T,2,P>
	{
	public:
      typedef PowerCompositeGridFunctionSpaceTraits<typename T::Traits::GridViewType, 
                                                    typename T::Traits::BackendType>
      Traits;

      //! Construct a PowerGridFunction with k clones of the function t
	  PowerGridFunctionSpace (T& t)
      {
        this->template setChild<0>(t);
        this->template setChild<1>(t);
        PowerGridFunctionSpaceBase<T,2,P>::setup();
      }

      //! Construct a PowerGridFunction with k clones of the function t
	  PowerGridFunctionSpace (T& t0, T& t1)
      {
        this->template setChild<0>(t0);
        this->template setChild<1>(t1);
        PowerGridFunctionSpaceBase<T,2,P>::setup();
      }
	};

	template<typename T, typename P>
	class PowerGridFunctionSpace<T,3,P> 
      : public PowerGridFunctionSpaceBase<T,3,P>
	{
	public:
      typedef PowerCompositeGridFunctionSpaceTraits<typename T::Traits::GridViewType, 
                                                    typename T::Traits::BackendType>
      Traits;

      //! Construct a PowerGridFunction with k clones of the function t
	  PowerGridFunctionSpace (T& t)
      {
        this->template setChild<0>(t);
        this->template setChild<1>(t);
        this->template setChild<2>(t);
        PowerGridFunctionSpaceBase<T,3,P>::setup();
      }

      //! Construct a PowerGridFunction with k clones of the function t
	  PowerGridFunctionSpace (T& t0, T& t1, T& t2)
      {
        this->template setChild<0>(t0);
        this->template setChild<1>(t1);
        this->template setChild<2>(t2);
        PowerGridFunctionSpaceBase<T,3,P>::setup();
      }
	};


    //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif
