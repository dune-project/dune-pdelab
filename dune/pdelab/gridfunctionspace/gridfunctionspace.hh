// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_HH

#include <cstddef>
#include<map>
#include <ostream>
#include<set>
#include<vector>

#include<dune/common/exceptions.hh>
#include<dune/common/geometrytype.hh>
#include<dune/common/stdstreams.hh>
#include<dune/grid/common/genericreferenceelements.hh>

#include <dune/localfunctions/common/localkey.hh>

#include"../common/countingptr.hh"
#include"../common/multitypetree.hh"
#include"../common/cpstoragepolicy.hh"
#include"../common/geometrywrapper.hh"

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
	template<typename G, typename L, typename C, typename B>
	struct GridFunctionSpaceTraits : public PowerCompositeGridFunctionSpaceTraits<G,B>
	{
	  //! \brief local finite element map
	  typedef L LocalFiniteElementMapType;

	  //! \brief local finite element
	  typedef typename L::Traits::LocalFiniteElementType LocalFiniteElementType;

	  //! \brief type representing constraints
	  typedef C ConstraintsType;
	};


    //! a class holding transformation for constrained spaces
    template<typename S, typename T>
    class ConstraintsTransformation 
      : public std::map<S,std::map<S,T> >
    {
    public:
      //! export ElementType
      typedef T ElementType;
      //! export RowType
      typedef std::map<S,T> RowType;
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
        typedef StdVectorBackend Backend;

		VectorContainer (const T& t) : BaseT(t.globalSize()) {}
		VectorContainer (const T& t, const E& e) : BaseT(t.globalSize(),e) {}
		VectorContainer& operator= (const E& e) // set all elements to same value
		{
		  for (typename BaseT::size_type i=0; i<BaseT::size(); i++)
			(*this)[i] = e;
		  return *this;
		}

        template<typename X>
        void std_copy_to (std::vector<X>& x) const
        {
          size_type n = this->size();
          x.resize(n);
          for (size_t i=0; i<n; i++)
            x[i] = (*this)[i];
        }
        
        template<typename X>
        void std_copy_from (const std::vector<X>& x)
        {
          size_t n = this->size();
          x.resize(n);
          for (size_t i=0; i<n; i++)
            (*this)[i] = x[i];
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
              index = is.subIndex(e,i,C);
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

	/** \brief class used to pass compile-time parameters to the implementation
     *
     *  This is the dummy default class which does nothing
     */
	struct GridFunctionGeneralMapper
	{
	  enum {dummy=0} ;
	};


    // Empty constraints assembler class
    class NoConstraints
    {
    public:
      enum { doBoundary = false };
      enum { doProcessor = false }; // added ParallelStuff
      enum { doSkeleton = false };
      enum { doVolume = false }; // might be necessary for cell-centered in parallel

      // methods are here just to show interfaces; they are never called because doX are false above
      template<typename F, typename I, typename LFS, typename T>
      void boundary (const F& f, const IntersectionGeometry<I>& ig, const LFS& lfs, T& trafo) const
      {
      }

      template<typename I, typename LFS, typename T>
      void processor (const IntersectionGeometry<I>& ig, const LFS& lfs, T& trafo) const
      {
      }

      template<typename I, typename LFS, typename T>
      void skeleton (const IntersectionGeometry<I>& ig, const LFS& lfs, T& trafo) const
      {
      }

      template<typename E, typename LFS, typename T>
      void volume (const ElementGeometry<E>& eg, const LFS& lfs, T& trafo) const
      {
      }

    };


	/** \brief mapper for layouts with arbitrary number of entries per entity
     *
     *  \tparam GV   Type implementing GridView
     *  \tparam LFEM Type implementing LocalFiniteElementMapInterface
     *  \tparam CE   Type for constraints assembler
     *  \tparam B    Backend type
     *  \tparam P    Parameter type
     */
	template<typename GV, typename LFEM, typename CE=NoConstraints, 
             typename B=StdVectorBackend, typename P=GridFunctionGeneralMapper>
	class GridFunctionSpace : public Countable, public LeafNode
	{
	public:
      //! export Traits class
	  typedef GridFunctionSpaceTraits<GV,LFEM,CE,B> Traits;
	  typedef typename GV::Traits::template Codim<0>::Entity Element;
	  typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;

	  //! extract type of container storing Ts
	  template<typename T>
	  struct VectorContainer
	  {
		//! \brief define Type as the Type of a container of E's
		typedef typename B::template VectorContainer<GridFunctionSpace,T> Type;	
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

      //! define local function space parametrized by self 
      typedef Dune::PDELab::LocalFunctionSpace<GridFunctionSpace> LocalFunctionSpace;

	  //! constructor
	  GridFunctionSpace (const GV& gridview, const LFEM& lfem, const CE& ce_) 
		: defaultce(ce_), gv(gridview), plfem(&lfem), ce(ce_)
	  {
		update();
	  }

	  //! constructor
	  GridFunctionSpace (const GV& gridview, const LFEM& lfem) 
        : gv(gridview), plfem(&lfem), ce(defaultce)
	  {
		update();
	  }

	  //! get grid view
	  const GV& gridview () const
	  {
		return gv;
	  }

	  // get finite element map, I think we dont need it
	  const LFEM& localFiniteElementMap () const
	  {
		return *plfem;
	  }

	  //! get dimension of finite element space
	  typename Traits::SizeType globalSize () const
	  {
		return nglobal;
	  }

	  //! get max dimension of shape function space
	  typename Traits::SizeType maxLocalSize () const
	  {
		return nlocal;
	  }

	  //! map index [0,globalSize-1] to root index set
	  typename Traits::SizeType upMap (typename Traits::SizeType i) const
	  {
		return i;
	  }

      // return constraints engine
      const typename Traits::ConstraintsType& constraints () const
      {
        return ce;
      }

	  //! compute global indices for one element
	  void globalIndices (const typename Traits::LocalFiniteElementType& lfe, 
                          const Element& e, 
						  std::vector<typename Traits::SizeType>& global) const
	  {
		// get layout of entity
		const typename Traits::LocalFiniteElementType::Traits::LocalCoefficientsType&
		  lc = lfe.localCoefficients();
		global.resize(lc.size());

		for (std::size_t i=0; i<lc.size(); ++i)
		  {
			// get geometry type of subentity 
			Dune::GeometryType gt=Dune::GenericReferenceElements<double,GV::Grid::dimension>
			  ::general(lfe.type()).type(lc.localKey(i).subEntity(),lc.localKey(i).codim());

			// evaluate consecutive index of subentity
			int index = eval_subindex<GV::Grid::dimension>(gv.indexSet(),e,
														   lc.localKey(i).subEntity(),
														   lc.localKey(i).codim());
		
			// now compute 
			global[i] = offset[(gtoffset.find(gt)->second)+index]+lc.localKey(i).index();
		  }
	  }

      // global Indices from element, needs additional finite element lookup
	  void globalIndices (const Element& e,
						  std::vector<typename Traits::SizeType>& global) const
      {
        globalIndices(plfem->find(e),e,global);
      }

      //------------------------------
      // generic data handle interface
      //------------------------------

      //! returns true if data for this codim should be communicated
      bool dataHandleContains (int dim, int codim) const
      {
        return (codimUsed.find(codim)!=codimUsed.end());
      }
      
      //! returns true if size per entity of given dim and codim is a constant
      bool dataHandleFixedSize (int dim, int codim) const
      {
        return false;
      }
      
      /*! how many objects of type DataType have to be sent for a given entity
        
        Note: Only the sender side needs to know this size. 
      */
      template<class EntityType>
      size_t dataHandleSize (const EntityType& e) const
      {
        Dune::GeometryType gt=e.type();
        typename GV::IndexSet::IndexType index = gtoffset.find(gt)->second + gv.indexSet().index(e);
        return offset[index+1]-offset[index];
      }

      //! return vector of global indices associated with the given entity
      template<class EntityType>
      void dataHandleGlobalIndices (const EntityType& e, 
                                    std::vector<typename Traits::SizeType>& global) const
      {
        Dune::GeometryType gt=e.type();
        typename GV::IndexSet::IndexType index = gtoffset.find(gt)->second + gv.indexSet().index(e);
        unsigned int n = offset[index+1]-offset[index];
		global.resize(n);
        for (unsigned i=0; i<n; i++)
          global[i] = offset[index]+i;
      }

      //------------------------------


	  // update information, e.g. when grid has changed
	  void update ()
	  {
        Dune::dinfo << "GridFunctionSpace(general version):" << std::endl;

		// determine which geometry types are used
		// needs one traversal of the grid
		typedef std::set<Dune::GeometryType> GtUsedSetType;
		GtUsedSetType gtused;
        codimUsed.clear();
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
			for (std::size_t i=0; i<lc.size(); ++i)
			  {
				Dune::GeometryType gt=Dune::GenericReferenceElements<double,GV::Grid::dimension>
				  ::general(it->type()).type(lc.localKey(i).subEntity(),lc.localKey(i).codim());
				gtused.insert(gt);
                codimUsed.insert(GV::Grid::dimension-gt.dim());
			  }
		  }

		// now we can allocate one number per entity that holds degrees of freedom
		typename Traits::SizeType nentities = 0;
		gtoffset.clear();
		const typename GV::IndexSet& is=gv.indexSet();
		for (typename GtUsedSetType::iterator i=gtused.begin(); i!=gtused.end(); ++i)
		  {
			gtoffset[*i] = nentities;
            Dune::dinfo << *i << ": " << is.size(*i)
                        << " entries in offset vector at " << nentities
                        << std::endl;
			nentities += is.size(*i);
		  }
        nentities++; // add one additional dummy entry; this allows to compute size of last entity.
		offset.resize(nentities);
		for (typename std::vector<typename Traits::SizeType>::iterator i=offset.begin(); i!=offset.end(); ++i)
		  *i = 0;
        Dune::dinfo << "allocated offset vector with size " << offset.size()
                    << std::endl;

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
			for (std::size_t i=0; i<lc.size(); ++i)
			  {
				Dune::GeometryType gt=Dune::GenericReferenceElements<double,GV::Grid::dimension>
				  ::general(it->type()).type(lc.localKey(i).subEntity(),lc.localKey(i).codim());
				unsigned int index = gtoffset[gt]
				  +eval_subindex<GV::Grid::dimension>(is,*it,lc.localKey(i).subEntity(),lc.localKey(i).codim());
				offset[index] = std::max(offset[index], 
                                         typename Traits::SizeType(lc.localKey(i).index()+1));
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
        Dune::dinfo << "total number of dofs is " << nglobal << std::endl;
	  }

	private:
      CE defaultce;
	  const GV& gv;
	  CP<LFEM const> plfem;
	  typename Traits::SizeType nlocal;
	  typename Traits::SizeType nglobal;
      const CE& ce;

	  std::map<Dune::GeometryType,typename Traits::SizeType> gtoffset; // offset in vector for given geometry type
	  std::vector<typename Traits::SizeType> offset; // offset into big vector for each entity;
      std::set<unsigned int> codimUsed;
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
	template<typename GV, typename LFEM, typename CE, typename B> 
	class GridFunctionSpace<GV,LFEM,CE,B,GridFunctionRestrictedMapper> : 
	  public Countable, public LeafNode
	{
	public:
      //! export Traits class
	  typedef GridFunctionSpaceTraits<GV,LFEM,CE,B> Traits;
	  typedef typename GV::Traits::template Codim<0>::Entity Element;
	  typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;

      //! extract type of container storing Es
	  template<typename T>
	  struct VectorContainer
	  {
		//! \brief define Type as the Type of a container of E's
		typedef typename B::template VectorContainer<GridFunctionSpace,T> Type;	
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
      typedef Dune::PDELab::LocalFunctionSpace<GridFunctionSpace> LocalFunctionSpace;

	  // constructor
	  GridFunctionSpace (const GV& gridview, const LFEM& lfem, const CE& ce_) 
		: defaultce(ce_), gv(gridview), plfem(&lfem), ce(ce_)
	  {
		update();
	  }

	  // constructor
	  GridFunctionSpace (const GV& gridview, const LFEM& lfem) 
		: gv(gridview), plfem(&lfem), ce(defaultce)
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

      // return constraints engine
      const typename Traits::ConstraintsType& constraints () const
      {
        return ce;
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
			Dune::GeometryType gt=Dune::GenericReferenceElements<double,GV::Grid::dimension>
			  ::general(lfe.type()).type(lc.localKey(i).subEntity(),lc.localKey(i).codim());

			// evaluate consecutive index of subentity
			int index = eval_subindex<GV::Grid::dimension>(gv.indexSet(),e,
														   lc.localKey(i).subEntity(),
														   lc.localKey(i).codim());
		
			// now compute 
			global[i] = offset.find(gt)->second+index*dofcountmap.find(gt)->second+lc.localKey(i).index();
		  }
	  }

      // global Indices from element, needs additional finite element lookup
	  void globalIndices (const Element& e,
						  std::vector<typename Traits::SizeType>& global) const
      {
        globalIndices(plfem->find(e),e,global);
      }


      //------------------------------
      // generic data handle interface
      //------------------------------

      //! returns true if data for this codim should be communicated
      bool dataHandleContains (int dim, int codim) const
      {
        return (codimUsed.find(codim)!=codimUsed.end());
      }
      
      //! returns true if size per entity of given dim and codim is a constant
      bool dataHandleFixedSize (int dim, int codim) const
      {
        return true;
      }
      
      /*! how many objects of type DataType have to be sent for a given entity
        
        Note: Only the sender side needs to know this size. 
      */
      template<class EntityType>
      size_t dataHandleSize (const EntityType& e) const
      {
        Dune::GeometryType gt=e.type();
        return dofcountmap.find(gt)->second;
      }

      //! return vector of global indices associated with the given entity
      template<class EntityType>
      void dataHandleGlobalIndices (const EntityType& e, 
                                    std::vector<typename Traits::SizeType>& global) const
      {
        Dune::GeometryType gt=e.type();
        typename GV::IndexSet::IndexType index = gv.indexSet().index(e);
        unsigned int n = dofcountmap.find(gt)->second;
		global.resize(n);
        for (int i=0; i<n; i++) 
          global[i] = offset.find(gt)->second + index*n + i;
      }

      //------------------------------

	  // update information, e.g. when grid has changed
	  void update ()
	  {
        Dune::dinfo << "GridFunctionSpace(restricted version):" << std::endl;

		// clear counters
		dofcountmap.clear();
		nlocal = 0;
        codimUsed.clear();


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
				SubentityType subentity(lc.localKey(i).subEntity(),lc.localKey(i).codim());
				if (countmap.find(subentity)==countmap.end())
				  countmap[subentity] = 1;
				else
				  (countmap[subentity])++;
			  }

			// traverse the map and compare #dofs per geometry type
			for (typename CountMapType::iterator i=countmap.begin(); i!=countmap.end(); ++i)
			  {
				Dune::GeometryType gt=Dune::GenericReferenceElements<double,GV::Grid::dimension>
				  ::general(it->type()).type(Dune::get<0>(i->first),Dune::get<1>(i->first));
				typename DofCountMapType::iterator j=dofcountmap.find(gt);
				if (j==dofcountmap.end())
				  {
                    codimUsed.insert(GV::Grid::dimension-gt.dim());
					dofcountmap[gt] = i->second;
					continue;
				  }
				if (j->second != i->second)
				  {
					DUNE_THROW(Dune::NotImplemented, "non constant # dofs per geometry type, use general version instead");
				  }
			  }
		  }
        Dune::dinfo << "max local number of dofs = " << nlocal << std::endl;

		// print result
		for (typename DofCountMapType::iterator i=dofcountmap.begin(); i!=dofcountmap.end(); ++i)
          Dune::dinfo << i->first << " has " << i->second
                      << " degrees of freedom" << std::endl;

		// compute offsets
		nglobal = 0;
		offset.clear();
		const typename GV::IndexSet& is=gv.indexSet();
		for (typename DofCountMapType::iterator i=dofcountmap.begin(); i!=dofcountmap.end(); ++i)
		  {
			offset[i->first] = nglobal;
			nglobal += is.size(i->first)*(i->second); 
            Dune::dinfo << i->first << " offset now " << nglobal << std::endl;
		  }
        Dune::dinfo << "total number of dofs = " << nglobal << std::endl;
	  }

	private:
	  const GV& gv;
	  CP<LFEM const> plfem;

	  typename Traits::SizeType nlocal;
	  typename Traits::SizeType nglobal;
      CE defaultce;
      const CE& ce;

	  typedef std::map<Dune::GeometryType,typename Traits::SizeType> DofCountMapType;
	  DofCountMapType dofcountmap; // number of degrees of freedom per geometry type
	  std::map<Dune::GeometryType,typename Traits::SizeType> offset; // offset in vector for given geometry type
      std::set<unsigned int> codimUsed;
	};



	// Pass this class as last template argument to GridFunctionSpace
	// to select specialization for fixed number of degrees of freedom in intersections
    template<typename IIS>
	struct GridFunctionStaticSize
	{
	  enum {dummy=2} ;
      typedef IIS IntersectionIndexSet;
	};

    // dummy index set for intersection; used to have static size
    // grid function space without DOFs in intersections
    class DummyIntersectionIndexSet
    {
    public:
      typedef int IndexType;
 
      // number of intersections in index set 
      // (intersections do not have a geometry type)
      IndexType size () const
      {
        DUNE_THROW(Dune::Exception,"need IntersectionIndexSet for DOFs in intersections");
      }

      // number of intersections associated with given element
      template<typename Element>
      IndexType size (const Element& element) const
      {
        DUNE_THROW(Dune::Exception,"need IntersectionIndexSet for DOFs in intersections");
      }

      // get index assigned to intersection
      template<typename Intersection>
      IndexType index (const Intersection& intersection) const
      {
        DUNE_THROW(Dune::Exception,"need IntersectionIndexSet for DOFs in intersections");
      }

      // get index of i'th intersection of element 
      // (in order they are visited by intersection iterator)
      template<typename Element>
      IndexType subIndex (const Element& element, int i) const
      {
        DUNE_THROW(Dune::Exception,"need IntersectionIndexSet for DOFs in intersections");
      }
    };

    // define type that can be used for static sized GFS without DOFS in intersections
    typedef GridFunctionStaticSize<DummyIntersectionIndexSet> SimpleGridFunctionStaticSize;

	// specialization with restricted mapper
	// GV : Type implementing GridView
	// FEM  : Type implementing LocalFiniteElementMapInterface
	// B : Backend type
	template<typename GV, typename LFEM, typename CE, typename B, typename IIS> 
	class GridFunctionSpace<GV,LFEM,CE,B,GridFunctionStaticSize<IIS> > : 
	  public Countable, public LeafNode
	{
      typedef std::map<unsigned int,unsigned int> DofPerCodimMapType;
	public:
      //! export traits class
	  typedef GridFunctionSpaceTraits<GV,LFEM,CE,B> Traits;
	  typedef typename GV::Traits::template Codim<0>::Entity Element;
	  typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;

      //! extract type of container storing Es
	  template<typename T>
	  struct VectorContainer
	  {
		//! \brief define Type as the Type of a container of E's
		typedef typename B::template VectorContainer<GridFunctionSpace,T> Type;	
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
      typedef Dune::PDELab::LocalFunctionSpace<GridFunctionSpace> LocalFunctionSpace;

	  // constructors
	  GridFunctionSpace (const GV& gridview, const LFEM& lfem, const IIS& iis_, const CE& ce_) 
		: defaultce(ce_), gv(gridview), plfem(&lfem), iis(iis_), ce(ce_)
	  {
		update();
	  }

	  GridFunctionSpace (const GV& gridview, const LFEM& lfem, const IIS& iis_) 
		: gv(gridview), plfem(&lfem), iis(iis_), ce(defaultce)
	  {
		update();
	  }

	  GridFunctionSpace (const GV& gridview, const LFEM& lfem, const CE& ce_) 
		: defaultce(ce_), gv(gridview), plfem(&lfem), iis(dummyiis), ce(ce_)
	  {
		update();
	  }

	  GridFunctionSpace (const GV& gridview, const LFEM& lfem) 
		: gv(gridview), plfem(&lfem), iis(dummyiis), ce(defaultce)
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

      // return constraints engine
      const typename Traits::ConstraintsType& constraints () const
      {
        return ce;
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
            typename GV::IndexSet::IndexType index;
            unsigned int cd = lc.localKey(i).codim();
            unsigned int se = lc.localKey(i).subEntity();

			// evaluate consecutive index of subentity
            if (cd==Dune::intersectionCodim)
              index = iis.subIndex(e,se);
            else
              index = eval_subindex<GV::Grid::dimension>(gv.indexSet(),e,se,cd);

			// now compute 
			global[i] = offset.find(cd)->second + index * dofpercodim.find(cd)->second + lc.localKey(i).index();
		  }
	  }

      // global Indices from element, needs additional finite element lookup
	  void globalIndices (const Element& e,
						  std::vector<typename Traits::SizeType>& global) const
      {
        globalIndices(plfem->find(e),e,global);
      }


      //------------------------------
      // generic data handle interface
      //------------------------------

      //! returns true if data for this codim should be communicated
      bool dataHandleContains (int dim, int codim) const
      {
        return (dofpercodim.find(codim)!=dofpercodim.end());
      }
      
      //! returns true if size per entity of given dim and codim is a constant
      bool dataHandleFixedSize (int dim, int codim) const
      {
        return true;
      }
      
      /*! how many objects of type DataType have to be sent for a given entity
        
        Note: Only the sender side needs to know this size. 
      */
      template<class EntityType>
      size_t dataHandleSize (const EntityType& e) const
      {
        const int cd = EntityType::codimension;
        return dofpercodim.find(cd)->second;
      }

      //! return vector of global indices associated with the given entity
      template<class EntityType>
      void dataHandleGlobalIndices (const EntityType& e, 
                                    std::vector<typename Traits::SizeType>& global) const
      {
        const int cd = EntityType::codimension;
        typename GV::IndexSet::IndexType o = offset.find(cd)->second;
        typename GV::IndexSet::IndexType index = gv.indexSet().index(e);
        unsigned int n = dofpercodim.find(cd)->second;
		global.resize(n);
        for (int i=0; i<n; i++) 
          global[i] = o + index*n + i;
//         Dune::dinfo << "[" << gv.grid().comm().rank() << "]: "
//                     << " global indices "
//                     << " offset=" << o
//                     << " index=" << index
//                     << " n=" << n
//                     << std::endl;
      }

      //------------------------------

	  // update information, e.g. when grid has changed
	  void update ()
	  {
        Dune::dinfo << "GridFunctionSpace(static size version):" << std::endl;

        // analyse local coefficients of first element

        // check geometry type
        ElementIterator it = gv.template begin<0>();
        if ((plfem->find(*it)).type()!=it->type())
          DUNE_THROW(Exception, "geometry type mismatch in GridFunctionSpace");

        // get local coefficients for this entity
        const typename Traits::LocalFiniteElementType::Traits::LocalCoefficientsType&
          lc = (plfem->find(*it)).localCoefficients();

        // extract number of degrees of freedom per element
        nlocal = static_cast<typename Traits::SizeType>(lc.size());

        // count number of degrees of freedom per subentity (including intersections)
        typedef Dune::tuple<unsigned int, unsigned int> SubentityType;
        typedef std::map<SubentityType,unsigned int> CountMapType;
        CountMapType countmap;
        for (int i=0; i<lc.size(); ++i)
          {
            SubentityType subentity(lc.localKey(i).subEntity(),lc.localKey(i).codim());
            if (countmap.find(subentity)==countmap.end())
              countmap[subentity] = 1;
            else
              (countmap[subentity])++;
          }
       
        // compute number of degrees of freedom per codim
        dofpercodim.clear();
        for (typename CountMapType::iterator i=countmap.begin(); i!=countmap.end(); ++i)
          {
            unsigned int cd = Dune::get<1>(i->first);
            typename DofPerCodimMapType::iterator j=dofpercodim.find(cd);
            if (j==dofpercodim.end())
              {
                dofpercodim[cd] = i->second;
              }
            else if (j->second != i->second)
              {
                DUNE_THROW(Dune::NotImplemented, "non constant # dofs per codim in static size grid function space");
              }
          }
        for (typename DofPerCodimMapType::iterator j=dofpercodim.begin(); j!=dofpercodim.end(); ++j)
          Dune::dinfo << " " << j->second << " degrees of freedom in codim "
                      << j->first << std::endl;

        // now compute global size
        nglobal = 0;
        offset.clear();
        for (typename DofPerCodimMapType::iterator j=dofpercodim.begin(); j!=dofpercodim.end(); ++j)
          {
            typename Traits::SizeType n;
            offset[j->first] = nglobal;
            if (j->first==Dune::intersectionCodim)
              n = j->second*iis.size();
            else
              n = j->second*gv.size(j->first);
            Dune::dinfo << "codim=" << j->first << " offset=" << nglobal
                        << " size=" << n << std::endl;
            nglobal += n;
            if (j->first!=Dune::intersectionCodim)
              Dune::dwarn << "WARNING: cannot handle multiple geometry types "
                          << "in static size grid function space"
                          << std::endl;;
          }
	  }

	private:
      DummyIntersectionIndexSet dummyiis; // for version without intersection DOFs
	  const GV& gv;
	  CP<LFEM const> plfem;
      const IIS& iis;

	  typename Traits::SizeType nlocal;
	  typename Traits::SizeType nglobal;
      CE defaultce;
      const CE& ce;

      DofPerCodimMapType dofpercodim;
      std::map<unsigned int,typename Traits::SizeType> offset;
	};


    //=======================================
    // power grid function space
    //=======================================


	// this class may be used to pass compile-time
	// parameters to the implementation 
	struct GridFunctionSpaceLexicographicMapper
	{
	  enum {dummy=0} ;
	};

	// this class may be used to pass compile-time
	// parameters to the implementation 
	struct GridFunctionSpaceBlockwiseMapper
	{
	  enum {dummy=0} ;
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
      //! export traits class
      typedef PowerCompositeGridFunctionSpaceTraits<typename T::Traits::GridViewType, 
                                                    typename T::Traits::BackendType>
      Traits;

      //! extract type of container storing Es
	  template<typename E>
	  struct VectorContainer
	  {
		//! \brief define Type as the Type of a container of E's
		typedef typename Traits::BackendType::template VectorContainer<PowerGridFunctionSpaceBase,E> Type;	
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
      typedef Dune::PDELab::PowerLocalFunctionSpace<PowerGridFunctionSpaceBase> LocalFunctionSpace;


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

      //------------------------------
      // generic data handle interface
      //------------------------------

      //! returns true if data for this codim should be communicated
      bool dataHandleContains (int dim, int codim) const
      {
        for (int i=0; i<k; i++)
          if (this->getChild(i).dataHandleContains(dim,codim))
            return true;
        return false;
      }
      
      //! returns true if size per entity of given dim and codim is a constant
      bool dataHandleFixedSize (int dim, int codim) const
      {
        for (int i=0; i<k; i++)
          if (!this->getChild(i).dataHandleFixedSize(dim,codim))
            return false;
        return true;
      }
      
      /*! how many objects of type DataType have to be sent for a given entity
        
        Note: Only the sender side needs to know this size. 
      */
      template<class EntityType>
      size_t dataHandleSize (const EntityType& e) const
      {
        size_t n=0;
        for (int i=0; i<k; i++)
          n += this->getChild(i).dataHandleSize(e);
        return n;
      }

      //! return vector of global indices associated with the given entity
      template<class EntityType>
      void dataHandleGlobalIndices (const EntityType& e, 
                                    std::vector<typename Traits::SizeType>& global) const
      {
        size_t n=0;
        for (int i=0; i<k; i++)
          n += this->getChild(i).dataHandleSize(e);
        global.resize(n);
        n = 0;
        for (int i=0; i<k; i++)
          {
            this->getChild(i).dataHandleGlobalIndices(e,childglobal);
            for (size_t j=0; j<childglobal.size(); j++)
              global[n+j] = offset[i]+childglobal[j];
            n += childglobal.size();
          }          
      }

      //------------------------------

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
        Dune::dinfo << "power grid function space(lexicographic version):"
                    << std::endl;
        Dune::dinfo << "( ";
        offset[0] = 0;
        maxlocalsize = 0;
        for (int i=0; i<k; i++)
          {
            childSize[i] = this->getChild(i).globalSize();
            Dune::dinfo << childSize[i] << " ";
            offset[i+1] = offset[i]+childSize[i];
            maxlocalsize += this->getChild(i).maxLocalSize();
          }
        Dune::dinfo << ") total size = " << offset[k]
                    << " max local size = " << maxlocalsize
                    << std::endl;
        childglobal.resize(maxlocalsize);
      }

      typename Traits::SizeType childSize[k];
      typename Traits::SizeType offset[k+1];
      typename Traits::SizeType maxlocalsize;
      mutable std::vector<typename Traits::SizeType> childglobal;
	};


    // product of identical grid function spaces
    // base class that holds implementation of the methods
    // specialization for blockwise ordering
	template<typename T, int k>
	class PowerGridFunctionSpaceBase<T,k,GridFunctionSpaceBlockwiseMapper> 
      : public PowerNode<T,k,CountingPointerStoragePolicy>,
        public Countable
	{
      friend class PowerGridFunctionSpace<T,k,GridFunctionSpaceBlockwiseMapper>;

	public:
      //! export traits class
      typedef PowerCompositeGridFunctionSpaceTraits<typename T::Traits::GridViewType, 
                                                    typename T::Traits::BackendType>
      Traits;

      //! extract type of container storing Es
	  template<typename E>
	  struct VectorContainer
	  {
		//! \brief define Type as the Type of a container of E's
		typedef typename Traits::BackendType::template VectorContainer<PowerGridFunctionSpaceBase,E> Type;	
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
      typedef Dune::PDELab::PowerLocalFunctionSpace<PowerGridFunctionSpaceBase> LocalFunctionSpace;
      
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
		return j*k+i;
	  }

      //------------------------------
      // generic data handle interface
      //------------------------------

      //! returns true if data for this codim should be communicated
      bool dataHandleContains (int dim, int codim) const
      {
        for (int i=0; i<k; i++)
          if (this->getChild(i).dataHandleContains(dim,codim))
            return true;
        return false;
      }
      
      //! returns true if size per entity of given dim and codim is a constant
      bool dataHandleFixedSize (int dim, int codim) const
      {
        for (int i=0; i<k; i++)
          if (!this->getChild(i).dataHandleFixedSize(dim,codim))
            return false;
        return true;
      }
      
      /*! how many objects of type DataType have to be sent for a given entity
        
        Note: Only the sender side needs to know this size. 
      */
      template<class EntityType>
      size_t dataHandleSize (const EntityType& e) const
      {
        size_t n=0;
        for (int i=0; i<k; i++)
          n += this->getChild(i).dataHandleSize(e);
        return n;
      }

      //! return vector of global indices associated with the given entity
      template<class EntityType>
      void dataHandleGlobalIndices (const EntityType& e, 
                                    std::vector<typename Traits::SizeType>& global) const
      {
        size_t n=0;
        for (int i=0; i<k; i++)
          n += this->getChild(i).dataHandleSize(e);
        global.resize(n);
        n = 0;
        for (int i=0; i<k; i++)
          {
            this->getChild(i).dataHandleGlobalIndices(e,childglobal);
            for (size_t j=0; j<childglobal.size(); j++)
              global[n+j] = childglobal[j]*k+i;
            n += childglobal.size();
          }          
      }

      //------------------------------

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
        Dune::dinfo << "power grid function space(blockwise version):"
                    << std::endl;
        Dune::dinfo << "( ";
        offset[0] = 0;
        maxlocalsize = 0;
        for (int i=0; i<k; i++)
          {
            childSize[i] = this->getChild(i).globalSize();
            Dune::dinfo << childSize[i] << " ";
            offset[i+1] = offset[i]+childSize[i];
            maxlocalsize += this->getChild(i).maxLocalSize();
          }
        Dune::dinfo << ") total size = " << offset[k]
                    << " max local size = " << maxlocalsize
                    << std::endl;
        for (int i=1; i<k; i++)
          if (childSize[i]!=childSize[0])
            DUNE_THROW(Exception, "components must be of equal size");
        childglobal.resize(maxlocalsize);
      }

      typename Traits::SizeType childSize[k];
      typename Traits::SizeType offset[k+1];
      typename Traits::SizeType maxlocalsize;
      mutable std::vector<typename Traits::SizeType> childglobal;
	};


    // product of identical grid function spaces
    // the specializations of this class just set the members
    // all the methods are generic in the implementation
	template<typename T, int k, typename P=GridFunctionSpaceLexicographicMapper>
	class PowerGridFunctionSpace 
      : public PowerGridFunctionSpaceBase<T,k,P>
 	{
	public:
      //! export traits class
      typedef PowerCompositeGridFunctionSpaceTraits<typename T::Traits::GridViewType, 
                                                    typename T::Traits::BackendType>
      Traits;

      //! Construct a PowerGridFunction with k clones of the function t
	  PowerGridFunctionSpace (T& t)
      {
        for (int i=0; i<k; i++)
          setChild(i,t);
        PowerGridFunctionSpaceBase<T,k,P>::setup();
      }

	  PowerGridFunctionSpace (T** t)  
      {
         for (int i=0; i<k; i++)
           setChild(i,*(t[i]));
         PowerGridFunctionSpaceBase<T,k,P>::setup();
     }
	};

	template<typename T, typename P>
	class PowerGridFunctionSpace<T,2,P> 
      : public PowerGridFunctionSpaceBase<T,2,P>
	{
	public:
      //! export traits class
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
      //! export traits class
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

	template<typename T, typename P>
	class PowerGridFunctionSpace<T,4,P> 
      : public PowerGridFunctionSpaceBase<T,4,P>
	{
	public:
      //! export traits class
      typedef PowerCompositeGridFunctionSpaceTraits<typename T::Traits::GridViewType, 
                                                    typename T::Traits::BackendType>
      Traits;

      //! Construct a PowerGridFunction with k clones of the function t
	  PowerGridFunctionSpace (T& t)
      {
        this->template setChild<0>(t);
        this->template setChild<1>(t);
        this->template setChild<2>(t);
        this->template setChild<3>(t);
        PowerGridFunctionSpaceBase<T,4,P>::setup();
      }

      //! Construct a PowerGridFunction with k clones of the function t
	  PowerGridFunctionSpace (T& t0, T& t1, T& t2, T& t3)
      {
        this->template setChild<0>(t0);
        this->template setChild<1>(t1);
        this->template setChild<2>(t2);
        this->template setChild<3>(t3);
        PowerGridFunctionSpaceBase<T,4,P>::setup();
      }
	};

	template<typename T, typename P>
	class PowerGridFunctionSpace<T,5,P> 
      : public PowerGridFunctionSpaceBase<T,5,P>
	{
	public:
      //! export traits class
      typedef PowerCompositeGridFunctionSpaceTraits<typename T::Traits::GridViewType, 
                                                    typename T::Traits::BackendType>
      Traits;

      //! Construct a PowerGridFunction with k clones of the function t
	  PowerGridFunctionSpace (T& t)
      {
        this->template setChild<0>(t);
        this->template setChild<1>(t);
        this->template setChild<2>(t);
        this->template setChild<3>(t);
        this->template setChild<4>(t);
        PowerGridFunctionSpaceBase<T,5,P>::setup();
      }

      //! Construct a PowerGridFunction with k clones of the function t
	  PowerGridFunctionSpace (T& t0, T& t1, T& t2, T& t3, T& t4)
      {
        this->template setChild<0>(t0);
        this->template setChild<1>(t1);
        this->template setChild<2>(t2);
        this->template setChild<3>(t3);
        this->template setChild<4>(t4);
        PowerGridFunctionSpaceBase<T,5,P>::setup();
      }
	};

	template<typename T, typename P>
	class PowerGridFunctionSpace<T,6,P> 
      : public PowerGridFunctionSpaceBase<T,6,P>
	{
	public:
      //! export traits class
      typedef PowerCompositeGridFunctionSpaceTraits<typename T::Traits::GridViewType, 
                                                    typename T::Traits::BackendType>
      Traits;

      //! Construct a PowerGridFunction with k clones of the function t
	  PowerGridFunctionSpace (T& t)
      {
        this->template setChild<0>(t);
        this->template setChild<1>(t);
        this->template setChild<2>(t);
        this->template setChild<3>(t);
        this->template setChild<4>(t);
        this->template setChild<5>(t);
        PowerGridFunctionSpaceBase<T,6,P>::setup();
      }

      //! Construct a PowerGridFunction with k clones of the function t
	  PowerGridFunctionSpace (T& t0, T& t1, T& t2, T& t3, T& t4, T& t5)
      {
        this->template setChild<0>(t0);
        this->template setChild<1>(t1);
        this->template setChild<2>(t2);
        this->template setChild<3>(t3);
        this->template setChild<4>(t4);
        this->template setChild<5>(t5);
        PowerGridFunctionSpaceBase<T,6,P>::setup();
      }
	};

	template<typename T, typename P>
	class PowerGridFunctionSpace<T,7,P> 
      : public PowerGridFunctionSpaceBase<T,7,P>
	{
	public:
      //! export traits class
      typedef PowerCompositeGridFunctionSpaceTraits<typename T::Traits::GridViewType, 
                                                    typename T::Traits::BackendType>
      Traits;

      //! Construct a PowerGridFunction with k clones of the function t
	  PowerGridFunctionSpace (T& t)
      {
        this->template setChild<0>(t);
        this->template setChild<1>(t);
        this->template setChild<2>(t);
        this->template setChild<3>(t);
        this->template setChild<4>(t);
        this->template setChild<5>(t);
        this->template setChild<6>(t);
        PowerGridFunctionSpaceBase<T,7,P>::setup();
      }

      //! Construct a PowerGridFunction with k clones of the function t
	  PowerGridFunctionSpace (T& t0, T& t1, T& t2, T& t3, T& t4, T& t5, T& t6)
      {
        this->template setChild<0>(t0);
        this->template setChild<1>(t1);
        this->template setChild<2>(t2);
        this->template setChild<3>(t3);
        this->template setChild<4>(t4);
        this->template setChild<5>(t5);
        this->template setChild<6>(t6);
        PowerGridFunctionSpaceBase<T,7,P>::setup();
      }
	};

	template<typename T, typename P>
	class PowerGridFunctionSpace<T,8,P> 
      : public PowerGridFunctionSpaceBase<T,8,P>
	{
	public:
      //! export traits class
      typedef PowerCompositeGridFunctionSpaceTraits<typename T::Traits::GridViewType, 
                                                    typename T::Traits::BackendType>
      Traits;

      //! Construct a PowerGridFunction with k clones of the function t
	  PowerGridFunctionSpace (T& t)
      {
        this->template setChild<0>(t);
        this->template setChild<1>(t);
        this->template setChild<2>(t);
        this->template setChild<3>(t);
        this->template setChild<4>(t);
        this->template setChild<5>(t);
        this->template setChild<6>(t);
        this->template setChild<7>(t);
        PowerGridFunctionSpaceBase<T,8,P>::setup();
      }

      //! Construct a PowerGridFunction with k clones of the function t
	  PowerGridFunctionSpace (T& t0, T& t1, T& t2, T& t3, T& t4, T& t5, T& t6, T& t7)
      {
        this->template setChild<0>(t0);
        this->template setChild<1>(t1);
        this->template setChild<2>(t2);
        this->template setChild<3>(t3);
        this->template setChild<4>(t4);
        this->template setChild<5>(t5);
        this->template setChild<6>(t6);
        this->template setChild<7>(t7);
        PowerGridFunctionSpaceBase<T,8,P>::setup();
      }
	};

	template<typename T, typename P>
	class PowerGridFunctionSpace<T,9,P> 
      : public PowerGridFunctionSpaceBase<T,9,P>
	{
	public:
      //! export traits class
      typedef PowerCompositeGridFunctionSpaceTraits<typename T::Traits::GridViewType, 
                                                    typename T::Traits::BackendType>
      Traits;

      //! Construct a PowerGridFunction with k clones of the function t
	  PowerGridFunctionSpace (T& t)
      {
        this->template setChild<0>(t);
        this->template setChild<1>(t);
        this->template setChild<2>(t);
        this->template setChild<3>(t);
        this->template setChild<4>(t);
        this->template setChild<5>(t);
        this->template setChild<6>(t);
        this->template setChild<7>(t);
        this->template setChild<8>(t);
        PowerGridFunctionSpaceBase<T,9,P>::setup();
      }

      //! Construct a PowerGridFunction with k clones of the function t
	  PowerGridFunctionSpace (T& t0, T& t1, T& t2, T& t3, T& t4, T& t5, T& t6, T& t7, T& t8)
      {
        this->template setChild<0>(t0);
        this->template setChild<1>(t1);
        this->template setChild<2>(t2);
        this->template setChild<3>(t3);
        this->template setChild<4>(t4);
        this->template setChild<5>(t5);
        this->template setChild<6>(t6);
        this->template setChild<7>(t7);
        this->template setChild<8>(t8);
        PowerGridFunctionSpaceBase<T,9,P>::setup();
      }
	};

	template<typename T, typename P>
	class PowerGridFunctionSpace<T,10,P> 
      : public PowerGridFunctionSpaceBase<T,10,P>
	{
	public:
      //! export traits class
      typedef PowerCompositeGridFunctionSpaceTraits<typename T::Traits::GridViewType, 
                                                    typename T::Traits::BackendType>
      Traits;

      //! Construct a PowerGridFunction with k clones of the function t
	  PowerGridFunctionSpace (T& t)
      {
        this->template setChild<0>(t);
        this->template setChild<1>(t);
        this->template setChild<2>(t);
        this->template setChild<3>(t);
        this->template setChild<4>(t);
        this->template setChild<5>(t);
        this->template setChild<6>(t);
        this->template setChild<7>(t);
        this->template setChild<8>(t);
        this->template setChild<9>(t);
        PowerGridFunctionSpaceBase<T,10,P>::setup();
      }

      //! Construct a PowerGridFunction with k clones of the function t
	  PowerGridFunctionSpace (T& t0, T& t1, T& t2, T& t3, T& t4, T& t5, T& t6, T& t7, T& t8, T& t9)
      {
        this->template setChild<0>(t0);
        this->template setChild<1>(t1);
        this->template setChild<2>(t2);
        this->template setChild<3>(t3);
        this->template setChild<4>(t4);
        this->template setChild<5>(t5);
        this->template setChild<6>(t6);
        this->template setChild<7>(t7);
        this->template setChild<8>(t8);
        this->template setChild<9>(t9);
        PowerGridFunctionSpaceBase<T,10,P>::setup();
      }
	};


    //=======================================
    // composite grid function space
    //=======================================


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

	template<typename P, typename T0, typename T1, typename T2, typename T3,
			 typename T4, typename T5, typename T6, typename T7, typename T8>
	class CompositeGridFunctionSpace;


    // tupel of grid function spaces
    // base class that holds implementation of the methods
    // this is the default version with lexicographic ordering
    // P is the ordering parameter
    // Ti are all grid function spaces
	template<typename P, typename T0, typename T1, typename T2=EmptyChild, typename T3=EmptyChild,
			 typename T4=EmptyChild, typename T5=EmptyChild, typename T6=EmptyChild,
			 typename T7=EmptyChild, typename T8=EmptyChild>
	class CompositeGridFunctionSpaceBase
	  : public CompositeNode<CountingPointerStoragePolicy,T0,T1,T2,T3,T4,T5,T6,T7,T8>,
		public Countable
	{
      friend class CompositeGridFunctionSpace<P,T0,T1,T2,T3,T4,T5,T6,T7,T8>; // for setup

      typedef CompositeNode<CountingPointerStoragePolicy,T0,T1,T2,T3,T4,T5,T6,T7,T8> BaseT;

	public:
      //! export traits class
      typedef PowerCompositeGridFunctionSpaceTraits<typename T0::Traits::GridViewType, 
                                                    typename T0::Traits::BackendType>
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
      typedef Dune::PDELab::CompositeLocalFunctionSpace<CompositeGridFunctionSpaceBase> LocalFunctionSpace;


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

	  // get dimension of finite element space
	  typename Traits::SizeType globalSize () const
	  {
        // this is bullshit all children may have different
        // size although they have the same type ...
		return offset[BaseT::CHILDREN];
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

    private:
      void setup ()
      {
        Dune::dinfo << "composite grid function space(lexicographic version):"
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
			 typename T4, typename T5, typename T6,
			 typename T7, typename T8>
	class CompositeGridFunctionSpaceBase<GridFunctionSpaceBlockwiseMapper,T0,T1,T2,T3,T4,T5,T6,T7,T8>
	  : public CompositeNode<CountingPointerStoragePolicy,T0,T1,T2,T3,T4,T5,T6,T7,T8>,
		public Countable
	{
      friend class CompositeGridFunctionSpace<GridFunctionSpaceBlockwiseMapper,
                                              T0,T1,T2,T3,T4,T5,T6,T7,T8>; // for setup

      typedef CompositeNode<CountingPointerStoragePolicy,T0,T1,T2,T3,T4,T5,T6,T7,T8> BaseT;

	public:
      //! export traits class
      typedef PowerCompositeGridFunctionSpaceTraits<typename T0::Traits::GridViewType, 
                                                    typename T0::Traits::BackendType>
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
      typedef Dune::PDELab::CompositeLocalFunctionSpace<CompositeGridFunctionSpaceBase> LocalFunctionSpace;

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

	  // get dimension of finite element space
	  typename Traits::SizeType globalSize () const
	  {
        // this is bullshit all children may have different
        // size although they have the same type ...
		return offset[BaseT::CHILDREN];
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
		return j*BaseT::CHILDREN+i;
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

    private:
      void setup ()
      {
        Dune::dinfo << "composite grid function space(blockwise version):"
                    << std::endl;

        CompositeGridFunctionSpaceBaseVisitChildMetaProgram<CompositeGridFunctionSpaceBase,BaseT::CHILDREN,0>::
          setup(*this,childGlobalSize,childLocalSize);

        for (int i=1; i<BaseT::CHILDREN; i++)
          if (childGlobalSize[i]!=childGlobalSize[0])
            DUNE_THROW(Exception, "components must be of equal size");

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

      typename Traits::SizeType childGlobalSize[BaseT::CHILDREN];
      typename Traits::SizeType childLocalSize[BaseT::CHILDREN];
      typename Traits::SizeType offset[BaseT::CHILDREN+1];
      typename Traits::SizeType maxlocalsize;
      mutable std::vector<typename Traits::SizeType> childglobal;
	};




    // composite grid function space primary template
	template<typename P, typename T0, typename T1, typename T2=EmptyChild, typename T3=EmptyChild,
			 typename T4=EmptyChild, typename T5=EmptyChild, typename T6=EmptyChild,
			 typename T7=EmptyChild, typename T8=EmptyChild>
	class CompositeGridFunctionSpace : 
      public CompositeGridFunctionSpaceBase<P,T0,T1,T2,T3,T4,T5,T6,T7,T8>
    {
      typedef CompositeGridFunctionSpaceBase<P,T0,T1,T2,T3,T4,T5,T6,T7,T8> BaseT;

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


    //=======================================
    // Subspace construction
    //=======================================

    template<typename GFS, int k, typename CGFS> // primary template, only specializations are used !
    class GridFunctionSubSpaceBase
    {
    };


    // CGFS is a composite
	template<typename GFS, int k, typename P, typename T0, typename T1, typename T2, typename T3,
			 typename T4, typename T5, typename T6, typename T7, typename T8>
    class GridFunctionSubSpaceBase<GFS,k, CompositeGridFunctionSpace<P,T0,T1,T2,T3,T4,T5,T6,T7,T8> >
      : public Countable // behave like child k of gfs which is a composite grid function space
    {
      typedef CompositeGridFunctionSpace<P,T0,T1,T2,T3,T4,T5,T6,T7,T8> CGFS;

    public:
      //! export traits class
	  typedef typename CGFS::Traits Traits;

      GridFunctionSubSpaceBase (const GFS& gfs)
        : pgfs(&gfs), pcgfs(&gfs.template getChild<k>())
      {
      }

	  enum { isLeaf = CGFS::isLeaf };
	  enum { isPower = CGFS::isPower /**< */ };
	  enum { isComposite = CGFS::isComposite /**< */ };
	  enum { CHILDREN = CGFS::CHILDREN };

	  template<int i>
	  struct Child
	  {
		typedef typename CGFS::template Child<i>::Type Type;
	  };

	  template<int i>
	  const typename CGFS::template Child<i>::Type& getChild () const
	  {
		return pcgfs->template getChild<i>();
	  }

	  // extract type of container storing Es
	  template<typename E>
	  struct VectorContainer
	  {
		//! \brief define Type as the Type of a container of E's
		typedef typename Traits::BackendType::template VectorContainer<GridFunctionSubSpaceBase,E> Type;	
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
      typedef CompositeLocalFunctionSpace<GridFunctionSubSpaceBase> LocalFunctionSpace;

	  // get grid view
	  const typename Traits::GridViewType& gridview () const
	  {
		return pgfs->gridview();
	  }

	  // get dimension of finite element space
	  typename Traits::SizeType globalSize () const
	  {
        return pgfs->globalSize();
	  }

	  // get max dimension of shape function space
	  typename Traits::SizeType maxLocalSize () const
	  {
		return pcgfs->maxLocalSize();
	  }

	  // map from gfs.child<k> to gfs
	  typename Traits::SizeType upMap (typename Traits::SizeType i) const
	  {
		return pgfs->upMap(pgfs->template subMap<k>(i));
	  }

	  // map from gfs.child<k>.child<i> to gfs.child<k>
      template<int i>
	  typename Traits::SizeType subMap (typename Traits::SizeType j) const
	  {
		return pcgfs->template subMap<i>(j);
	  }

    private:
      CP<GFS const> pgfs;
      CP<CGFS const> pcgfs;
    };


    // CGFS is a power
	template<typename GFS, int k, typename T, int l, typename P>
    class GridFunctionSubSpaceBase<GFS,k, PowerGridFunctionSpace<T,l,P> >
      : public Countable // behave like child k of gfs which is a composite grid function space
    {
      typedef PowerGridFunctionSpace<T,l,P> CGFS;

    public:
      //! export traits class
	  typedef typename CGFS::Traits Traits;

      GridFunctionSubSpaceBase (const GFS& gfs)
        : pgfs(&gfs), pcgfs(&gfs.template getChild<k>())
      {
      }

	  enum { isLeaf = CGFS::isLeaf };
	  enum { isPower = CGFS::isPower /**< */ };
	  enum { isComposite = CGFS::isComposite /**< */ };
	  enum { CHILDREN = CGFS::CHILDREN };

	  template<int i>
	  struct Child
	  {
		typedef typename CGFS::template Child<i>::Type Type;
	  };

	  template<int i>
	  const typename CGFS::template Child<i>::Type& getChild () const
	  {
		return pcgfs->template getChild<i>();
	  }

	  const T& getChild (int i) const
	  {
		return pcgfs->getChild(i);
	  }

	  // extract type of container storing Es
	  template<typename E>
	  struct VectorContainer
	  {
		//! \brief define Type as the Type of a container of E's
		typedef typename Traits::BackendType::template VectorContainer<GridFunctionSubSpaceBase,E> Type;	
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
      typedef PowerLocalFunctionSpace<GridFunctionSubSpaceBase> LocalFunctionSpace;

	  // get grid view
	  const typename Traits::GridViewType& gridview () const
	  {
		return pgfs->gridview();
	  }

	  // get dimension of finite element space
	  typename Traits::SizeType globalSize () const
	  {
        return pgfs->globalSize();
	  }

	  // get max dimension of shape function space
	  typename Traits::SizeType maxLocalSize () const
	  {
		return pcgfs->maxLocalSize();
	  }

	  // map from gfs.child<k> to gfs
	  typename Traits::SizeType upMap (typename Traits::SizeType i) const
	  {
		return pgfs->upMap(pgfs->template subMap<k>(i));
	  }

	  // map from gfs.child<k>.child<i> to gfs.child<k>
      template<int i>
	  typename Traits::SizeType subMap (typename Traits::SizeType j) const
	  {
		return pcgfs->template subMap<i>(j);
	  }

    private:
      CP<GFS const> pgfs;
      CP<CGFS const> pcgfs;
    };


    // CGFS is a leaf
	template<typename GFS, int k, typename GV, typename LFEM, typename CE, typename B, typename P>
    class GridFunctionSubSpaceBase<GFS,k, GridFunctionSpace<GV,LFEM,CE,B,P> >
      : public Countable // behave like child k of GFS which is a grid function space
    {
      typedef GridFunctionSpace<GV,LFEM,CE,B,P> CGFS;

    public:
      //! export traits class
	  typedef typename CGFS::Traits Traits;
	  typedef typename GV::Traits::template Codim<0>::Entity Element;

      GridFunctionSubSpaceBase (const GFS& gfs)
        : pgfs(&gfs), pcgfs(&gfs.template getChild<k>())
      {
      }

	  enum { isLeaf = CGFS::isLeaf };
	  enum { isPower = CGFS::isPower /**< */ };
	  enum { isComposite = CGFS::isComposite /**< */ };
	  enum { CHILDREN = CGFS::CHILDREN };

      //! extract type of container storing Es
	  template<typename E>
	  struct VectorContainer
	  {
		//! \brief define Type as the Type of a container of E's
		typedef typename Traits::BackendType::template VectorContainer<GridFunctionSubSpaceBase,E> Type;	
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
      typedef Dune::PDELab::LocalFunctionSpace<GridFunctionSubSpaceBase> LocalFunctionSpace;

	  // get grid view
	  const typename Traits::GridViewType& gridview () const
	  {
		return pcgfs->gridview();
	  }

	  // get finite element map
	  const LFEM& localFiniteElementMap () const
	  {
		return pcgfs->localFiniteElementMap();
	  }

	  // get dimension of finite element space
	  typename Traits::SizeType globalSize () const
	  {
        return pgfs->globalSize();
	  }

	  // get max dimension of shape function space
	  typename Traits::SizeType maxLocalSize () const
	  {
		return pcgfs->maxLocalSize();
	  }

	  // map from gfs.child<k> to gfs
	  typename Traits::SizeType upMap (typename Traits::SizeType i) const
	  {
		return pgfs->upMap(pgfs->template subMap<k>(i));
	  }

	  // compute global indices for one element
	  void globalIndices (const typename Traits::LocalFiniteElementType& lfe, 
                          const Element& e, 
						  std::vector<typename Traits::SizeType>& global) const
	  {
        pcgfs->globalIndices(lfe,e,global);
	  }

      // global Indices from element, needs additional finite element lookup
	  void globalIndices (const Element& e,
						  std::vector<typename Traits::SizeType>& global) const
      {
        pcgfs->globalIndices(e,global);
      }

    private:
      CP<GFS const> pgfs;
      CP<CGFS const> pcgfs;
    };

    // ensure that GFS is not a leaf
    template<typename GFS, int k, int isleaf>
    class GridFunctionSubSpaceIntermediateBase 
      : public GridFunctionSubSpaceBase<GFS,k,typename GFS::template Child<k>::Type>
    {
      typedef GridFunctionSubSpaceBase<GFS,k,typename GFS::template Child<k>::Type> BaseT;
    public:
      GridFunctionSubSpaceIntermediateBase (const GFS& gfs) : BaseT(gfs)
      {
      }
    };


    // compilation error if subspace from a leaf is taken
    template<typename GFS, int k>
    class GridFunctionSubSpaceIntermediateBase<GFS,k,true> : public Countable
    {
    public:
      GridFunctionSubSpaceIntermediateBase (const GFS& gfs)
      {
 		dune_static_assert((static_cast<int>(GFS::isLeaf)==0),"subspace cannot be taken from a leaf");
      }
    };




    template<typename GFS, int k>
    class GridFunctionSubSpace : public GridFunctionSubSpaceIntermediateBase<GFS,k,GFS::isLeaf>
    {
    public:
      GridFunctionSubSpace (const GFS& gfs) 
        : GridFunctionSubSpaceIntermediateBase<GFS,k,GFS::isLeaf>(gfs)
      {
        Dune::dinfo << "grid function subspace:" << std::endl;
        Dune::dinfo << "root space size = " << gfs.globalSize()
                    << " max local size = " << this->maxLocalSize()
                    << std::endl;
      }
    };

    //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif
