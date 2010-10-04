// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_HH

#include <cstddef>
#include <map>
#include <ostream>
#include <set>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/geometrytype.hh>
#include <dune/common/static_assert.hh>
#include <dune/common/stdstreams.hh>
#include <dune/common/typetraits.hh>

#include <dune/grid/common/genericreferenceelements.hh>

#include <dune/localfunctions/common/localkey.hh>

#include "../common/countingptr.hh"
#include "../common/multitypetree.hh"
#include "../common/cpstoragepolicy.hh"
#include "../common/geometrywrapper.hh"
#include <dune/pdelab/finiteelement/traits.hh>

#include"localfunctionspace.hh"

namespace Dune {
  namespace PDELab {

    //! \addtogroup GridFunctionSpace grid function space
    //! \ingroup PDELab
    //! \{

    //=======================================
    // grid function space : single component case
    //=======================================

	template<typename G, typename B, typename M, int k>
	struct PowerCompositeGridFunctionSpaceTraits
	{
      enum{ 
        //! \brief True if this grid function space is composed of others.
        isComposite = 1,
        //! \brief number of child spaces
        noChilds = k
      };
      
	  //! \brief the grid view where grid function is defined upon
	  typedef G GridViewType;

	  //! \brief vector backend
	  typedef B BackendType;

      //! \brief mapper
      typedef M MapperType;
      
	  //! \brief short cut for size type exported by Backend
	  typedef typename B::size_type SizeType;
	};

    //! collect types exported by a leaf grid function space
    /**
     * This is based on a global FiniteElementMap
     */
    template<typename G, typename L, typename C, typename B, class = void>
    struct GridFunctionSpaceTraits
    {
      //! True if this grid function space is composed of others.
      static const bool isComposite = false;

      //! the grid view where grid function is defined upon
      typedef G GridViewType;

      //! vector backend
      typedef B BackendType;

      //! short cut for size type exported by Backend
      typedef typename B::size_type SizeType;

      //! local finite element map
      typedef L LocalFiniteElementMapType DUNE_DEPRECATED;

      //! finite element map
      typedef L FiniteElementMapType;

      //! local finite element
      typedef typename L::Traits::FiniteElementType LocalFiniteElementType
        DUNE_DEPRECATED;

      //! finite element
      typedef typename L::Traits::FiniteElementType FiniteElementType;

      //! type representing constraints
      typedef C ConstraintsType;
    };

	//! \brief collect types exported by a leaf grid function space
    /**
     * This is based on LocalFiniteElementMap
     */
	template<typename G, typename L, typename C, typename B>
    struct GridFunctionSpaceTraits<G, L, C, B,
             typename enable_if<AlwaysTrue<typename
                   L::Traits::LocalFiniteElementType>::value>::type>
	{
      enum{ 
        //! \brief True if this grid function space is composed of others.
        isComposite = 0
      };

	  //! \brief the grid view where grid function is defined upon
	  typedef G GridViewType;

	  //! \brief vector backend
	  typedef B BackendType;

	  //! \brief short cut for size type exported by Backend
	  typedef typename B::size_type SizeType;

	  //! \brief local finite element map
	  typedef L LocalFiniteElementMapType DUNE_DEPRECATED;

      //! \brief finite element map
      typedef L FiniteElementMapType;

	  //! \brief local finite element
      typedef typename L::Traits::LocalFiniteElementType LocalFiniteElementType
        DUNE_DEPRECATED;

      //! \brief finite element
      typedef typename L::Traits::LocalFiniteElementType FiniteElementType;

	  //! \brief type representing constraints
	  typedef C ConstraintsType;
	};

    //! \brief a class holding transformation for constrained spaces
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

	//! \brief Simple Backend for std::vector
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
	  static const typename C::value_type& access (const C& c, size_type i)
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

	/** \brief Tag indicating an arbitrary number of unkowns per entity.
     *
     * class used to pass compile-time parameter to the GridFunctionSpace.
     */
	struct GridFunctionGeneralMapper {};


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


	/** \brief A grid function space.
     *
     *  \tparam GV   Type implementing GridView
     *  \tparam FEM  Type implementing FiniteElementMapInterface
     *  \tparam CE   Type for constraints assembler
     *  \tparam B    Backend type
     *  \tparam P    Parameter type. Possible types are 
     * \link GridFunctionGeneralMapper \endlink (arbitrary number of unknowns per
     * entity) or \link GridFunctionRestrictedMapper \endlink (fixed number of unknowns per
     * entity) or \link GridFunctionStaticSize \endlink (number of unknowns per
     * entity, known at compile-time)
     */
    template<typename GV, typename FEM, typename CE=NoConstraints,
             typename B=StdVectorBackend, typename P=GridFunctionGeneralMapper>
	class GridFunctionSpace : public Countable, public LeafNode
	{
	public:
      //! export Traits class
	  typedef GridFunctionSpaceTraits<GV,FEM,CE,B> Traits;
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
      typedef Dune::PDELab::LeafLocalFunctionSpaceNode<GridFunctionSpace> LocalFunctionSpace;

	  //! constructor
	  GridFunctionSpace (const GV& gridview, const FEM& fem, const CE& ce_) 
		: defaultce(ce_), gv(gridview), pfem(&fem), ce(ce_)
	  {
		update();
	  }

	  //! constructor
	  GridFunctionSpace (const GV& gridview, const FEM& fem) 
        : gv(gridview), pfem(&fem), ce(defaultce)
	  {
		update();
	  }

	  //! get grid view
	  const GV& gridview () const
	  {
		return gv;
	  }

      // get finite element map, I think we dont need it
      const FEM& finiteElementMap () const
      {
        return *pfem;
      }

	  // get finite element map, I think we dont need it
      const FEM& localFiniteElementMap () const DUNE_DEPRECATED
	  {
        return *pfem;
	  }

      //! get dimension of root finite element space
	  typename Traits::SizeType globalSize () const
	  {
		return nglobal;
	  }

      //! get dimension of this finite element space
      typename Traits::SizeType size () const
      {
        return nglobal;
      }

	  //! get max dimension of shape function space
      //! \todo What are the exact semantics of maxLocalSize?
	  typename Traits::SizeType maxLocalSize () const
	  {
		return nlocal;
	  }

      //! map index from our index set [0,size()-1] to root index set
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
      void globalIndices (const typename Traits::FiniteElementType& fe,
                          const Element& e, 
						  std::vector<typename Traits::SizeType>& global) const
	  {
        typedef FiniteElementTraits<typename Traits::FiniteElementType>
          FETraits;
		// get layout of entity
        const typename FETraits::Coefficients &coeffs =
          FETraits::coefficients(fe);
        global.resize(coeffs.size());

        for (std::size_t i=0; i<coeffs.size(); ++i)
		  {
			// get geometry type of subentity 
			Dune::GeometryType gt=Dune::GenericReferenceElements<double,GV::Grid::dimension>
              ::general(fe.type()).type(coeffs.localKey(i).subEntity(),
                                        coeffs.localKey(i).codim());

			// evaluate consecutive index of subentity
            int index = gv.indexSet().subIndex(e,
                                               coeffs.localKey(i).subEntity(),
                                               coeffs.localKey(i).codim());
		
			// now compute 
            global[i] = offset[(gtoffset.find(gt)->second)+index]+
              coeffs.localKey(i).index();
		  }
	  }

      // global Indices from element, needs additional finite element lookup
	  void globalIndices (const Element& e,
						  std::vector<typename Traits::SizeType>& global) const
      {
        globalIndices(pfem->find(e),e,global);
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
            const typename Traits::FiniteElementType &fe = pfem->find(*it);
			// check geometry type
            if (fe.type()!=it->type())
			  DUNE_THROW(Exception, "geometry type mismatch in GridFunctionSpace");

			// get local coefficients for this entity
            typedef FiniteElementTraits<typename Traits::FiniteElementType>
              FETraits;
            const typename FETraits::Coefficients& coeffs =
              FETraits::coefficients(fe);

			// insert geometry type of all subentities into set
            for (std::size_t i=0; i<coeffs.size(); ++i)
			  {
				Dune::GeometryType gt=Dune::GenericReferenceElements<double,GV::Grid::dimension>
                  ::general(it->type()).type(coeffs.localKey(i).subEntity(),
                                             coeffs.localKey(i).codim());
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
            const typename Traits::FiniteElementType &fe = pfem->find(*it);

			// get local coefficients for this entity
            typedef FiniteElementTraits<typename Traits::FiniteElementType>
              FETraits;
            const typename FETraits::Coefficients& coeffs =
              FETraits::coefficients(fe);

			// compute maximum number of degrees of freedom per element
            nlocal = std::max(nlocal, static_cast<typename Traits::SizeType>
                                        (coeffs.size()));

			// compute maximum size for each subentity
            for (std::size_t i=0; i<coeffs.size(); ++i)
			  {
				Dune::GeometryType gt=Dune::GenericReferenceElements<double,GV::Grid::dimension>
                  ::general(it->type()).type(coeffs.localKey(i).subEntity(),
                                             coeffs.localKey(i).codim());
                unsigned int index = gtoffset[gt] +
                  is.subIndex(*it, coeffs.localKey(i).subEntity(),
                              coeffs.localKey(i).codim());
				offset[index] = std::max(offset[index], 
                                         typename Traits::SizeType
                                            (coeffs.localKey(i).index()+1));
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
	  CountingPointer<FEM const> pfem;
	  typename Traits::SizeType nlocal;
	  typename Traits::SizeType nglobal;
      const CE& ce;

	  std::map<Dune::GeometryType,typename Traits::SizeType> gtoffset; // offset in vector for given geometry type
	  std::vector<typename Traits::SizeType> offset; // offset into big vector for each entity;
      std::set<unsigned int> codimUsed;
	};

	/** \brief Tag indicating a fixed number of unkowns per entity (known at compile time).
     *
     * class used to pass compile-time parameter to the GridFunctionSpace.
     *
     */
	struct GridFunctionRestrictedMapper {};

    //! \}

	// specialization with restricted mapper
	// GV : Type implementing GridView
    // FEM  : Type implementing FiniteElementMapInterface
	// B : Backend type
    template<typename GV, typename FEM, typename CE, typename B>
    class GridFunctionSpace<GV,FEM,CE,B,GridFunctionRestrictedMapper> :
	  public Countable, public LeafNode
	{
	public:
      //! export Traits class
	  typedef GridFunctionSpaceTraits<GV,FEM,CE,B> Traits;
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
      typedef Dune::PDELab::LeafLocalFunctionSpaceNode<GridFunctionSpace> LocalFunctionSpace;

	  // constructor
      GridFunctionSpace (const GV& gridview, const FEM& fem, const CE& ce_)
        : gv(gridview), pfem(&fem), defaultce(ce_), ce(ce_)
	  {
		update();
	  }

	  // constructor
      GridFunctionSpace (const GV& gridview, const FEM& fem)
        : gv(gridview), pfem(&fem), ce(defaultce)
	  {
		update();
	  }

	  // get grid view
	  const GV& gridview () const
	  {
		return gv;
	  }

      // get finite element map, I think we dont need it
      const FEM& finiteElementMap () const
      {
        return *pfem;
      }

	  // get finite element map, I think we dont need it
      const FEM& localFiniteElementMap () const DUNE_DEPRECATED
	  {
        return *pfem;
	  }

      //! get dimension of root finite element space
	  typename Traits::SizeType globalSize () const
	  {
		return nglobal;
	  }

      //! get dimension of this finite element space
      typename Traits::SizeType size () const
      {
        return nglobal;
      }

	  // get max dimension of shape function space
	  typename Traits::SizeType maxLocalSize () const
	  {
		return nlocal;
	  }

      // map index [0,size()-1] to root index set
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
      void globalIndices (const typename Traits::FiniteElementType& fe,
                          const Element& e,
						  std::vector<typename Traits::SizeType>& global) const
	  {
		// get local coefficients for this entity
        typedef FiniteElementTraits<typename Traits::FiniteElementType>
          FETraits;
        const typename FETraits::Coefficients &coeffs =
          FETraits::coefficients(fe);
        global.resize(coeffs.size());

        for (unsigned int i=0; i<coeffs.size(); ++i)
		  {
			// get geometry type of subentity 
			Dune::GeometryType gt=Dune::GenericReferenceElements<double,GV::Grid::dimension>
              ::general(fe.type()).type(coeffs.localKey(i).subEntity(),
                                        coeffs.localKey(i).codim());

			// evaluate consecutive index of subentity
            int index = gv.indexSet().subIndex(e,
                                               coeffs.localKey(i).subEntity(),
                                               coeffs.localKey(i).codim());
		
			// now compute 
            global[i] =
              offset.find(gt)->second+index*dofcountmap.find(gt)->second
              + coeffs.localKey(i).index();
		  }
	  }

      // global Indices from element, needs additional finite element lookup
	  void globalIndices (const Element& e,
						  std::vector<typename Traits::SizeType>& global) const
      {
        globalIndices(pfem->find(e),e,global);
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
        for(unsigned i=0; i<n; i++)
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
            const typename Traits::FiniteElementType &fe = pfem->find(*it);
			// check geometry type
            if (fe.type()!=it->type())
			  DUNE_THROW(Exception, "geometry type mismatch in GridFunctionSpace");

			// get local coefficients for this entity
            typedef FiniteElementTraits<typename Traits::FiniteElementType>
              FETraits;
            const typename FETraits::Coefficients &coeffs =
              FETraits::coefficients(fe);

			// compute maximum number of degrees of freedom per element
            nlocal = std::max(nlocal, static_cast<typename Traits::SizeType>
                                       (coeffs.size()));

			// store count for each subentity in a map
			typedef Dune::tuple<unsigned int, unsigned int> SubentityType;
			typedef std::map<SubentityType,unsigned int> CountMapType;
			CountMapType countmap;

			// assume that key within each subentity is unique
            for (unsigned i=0; i<coeffs.size(); ++i)
			  {
                SubentityType subentity(coeffs.localKey(i).subEntity(),
                                        coeffs.localKey(i).codim());
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
      CountingPointer<FEM const> pfem;

	  typename Traits::SizeType nlocal;
	  typename Traits::SizeType nglobal;
      CE defaultce;
      const CE& ce;

	  typedef std::map<Dune::GeometryType,typename Traits::SizeType> DofCountMapType;
	  DofCountMapType dofcountmap; // number of degrees of freedom per geometry type
	  std::map<Dune::GeometryType,typename Traits::SizeType> offset; // offset in vector for given geometry type
      std::set<unsigned int> codimUsed;
	};

    //! \addtogroup GridFunctionSpace
    //! \{

    /** \brief Tag indicating a fixed number of unkowns per entity (known at compile time).
     *
     * class used to pass compile-time parameter to the GridFunctionSpace.
     *
     * \tparam IIS type of the index set for the intersections.
     * Use \link DummyIntersectionIndexSet \endlink if no unknowns are associated with intersections.
     */
    template<typename IIS>
	struct GridFunctionStaticSize
	{
      typedef IIS IntersectionIndexSet;
	};

    //!! \brief dummy index set for intersection for grid function spaces with static size
    //!  but without DOFs in intersections
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

    //! \brief type that can be used for static sized GFS without DOFS in intersections
    typedef GridFunctionStaticSize<DummyIntersectionIndexSet> SimpleGridFunctionStaticSize;
    
    //! \}

	// specialization with restricted mapper
	// GV : Type implementing GridView
    // FEM  : Type implementing FiniteElementMapInterface
	// B : Backend type
    template<typename GV, typename FEM, typename CE, typename B, typename IIS>
    class GridFunctionSpace<GV,FEM,CE,B,GridFunctionStaticSize<IIS> > :
	  public Countable, public LeafNode
	{
      typedef std::map<unsigned int,unsigned int> DofPerCodimMapType;
	public:
      //! export traits class
      typedef GridFunctionSpaceTraits<GV,FEM,CE,B> Traits;
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
      typedef Dune::PDELab::LeafLocalFunctionSpaceNode<GridFunctionSpace> LocalFunctionSpace;

	  // constructors
      GridFunctionSpace (const GV& gridview, const FEM& fem, const IIS& iis_,
                         const CE& ce_)
        : gv(gridview), pfem(&fem), iis(iis_), defaultce(ce_), ce(ce_)
	  {
		update();
	  }

      GridFunctionSpace (const GV& gridview, const FEM& fem, const IIS& iis_)
        : gv(gridview), pfem(&fem), iis(iis_), ce(defaultce)
	  {
		update();
	  }

      GridFunctionSpace (const GV& gridview, const FEM& fem, const CE& ce_)
        : gv(gridview), pfem(&fem), iis(dummyiis), defaultce(ce_), ce(ce_)
	  {
		update();
	  }

      GridFunctionSpace (const GV& gridview, const FEM& fem)
        : gv(gridview), pfem(&fem), iis(dummyiis), ce(defaultce)
	  {
		update();
	  }

	  // get grid view
	  const GV& gridview () const
	  {
		return gv;
	  }

      // get finite element map, I think we dont need it
      const FEM& finiteElementMap() const
      {
        return *pfem;
      }

	  // get finite element map, I think we dont need it
      const FEM& localFiniteElementMap() const DUNE_DEPRECATED
	  {
        return *pfem;
	  }

      //! get dimension of root finite element space
	  typename Traits::SizeType globalSize () const
	  {
		return nglobal;
	  }

      //! get dimension of this finite element space
      typename Traits::SizeType size () const
      {
        return nglobal;
      }

	  // get max dimension of shape function space
	  typename Traits::SizeType maxLocalSize () const
	  {
		return nlocal;
	  }

      //! map index from our index set [0,size()-1] to root index set
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
      void globalIndices (const typename Traits::FiniteElementType& fe,
                          const Element& e,
						  std::vector<typename Traits::SizeType>& global) const
	  {
        typedef FiniteElementTraits<typename Traits::FiniteElementType>
          FETraits;
		// get local coefficients for this entity
        const typename FETraits::Coefficients &coeffs =
          FETraits::coefficients(fe);
        global.resize(coeffs.size());

        for (unsigned int i=0; i<coeffs.size(); ++i)
		  {
            typename GV::IndexSet::IndexType index;
            unsigned int cd = coeffs.localKey(i).codim();
            unsigned int se = coeffs.localKey(i).subEntity();

			// evaluate consecutive index of subentity
            if (cd==Dune::LocalKey::intersectionCodim)
              index = iis.subIndex(e,se);
            else
              index = gv.indexSet().subIndex(e,se,cd);

			// now compute 
            global[i] =
              offset.find(cd)->second + index * dofpercodim.find(cd)->second
              + coeffs.localKey(i).index();
		  }
	  }

      // global Indices from element, needs additional finite element lookup
	  void globalIndices (const Element& e,
						  std::vector<typename Traits::SizeType>& global) const
      {
        globalIndices(pfem->find(e),e,global);
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
        for (unsigned int i=0; i<n; i++) 
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
        typename Traits::FiniteElementType &fe = pfem->find(*it);
        if (fe.type()!=it->type())
          DUNE_THROW(Exception, "geometry type mismatch in GridFunctionSpace");

        // get local coefficients for this entity
        typedef FiniteElementTraits<typename Traits::FiniteElementType>
          FETraits;
        const typename FETraits::Coefficients &coeffs =
          FETraits::coefficients(fe);

        // extract number of degrees of freedom per element
        nlocal = static_cast<typename Traits::SizeType>(coeffs.size());

        // count number of degrees of freedom per subentity (including intersections)
        typedef Dune::tuple<unsigned int, unsigned int> SubentityType;
        typedef std::map<SubentityType,unsigned int> CountMapType;
        CountMapType countmap;
        for (unsigned int i=0; i<coeffs.size(); ++i)
          {
            SubentityType subentity(coeffs.localKey(i).subEntity(),
                                    coeffs.localKey(i).codim());
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
            if (j->first==Dune::LocalKey::intersectionCodim)
              n = j->second*iis.size();
            else
              n = j->second*gv.size(j->first);
            Dune::dinfo << "codim=" << j->first << " offset=" << nglobal
                        << " size=" << n << std::endl;
            nglobal += n;
//             if (j->first!=Dune::LocalKey::intersectionCodim)
//               Dune::dwarn << "WARNING: cannot handle multiple geometry types "
//                           << "in static size grid function space"
//                           << std::endl;;
          }
	  }

	private:
      DummyIntersectionIndexSet dummyiis; // for version without intersection DOFs
	  const GV& gv;
      CountingPointer<FEM const> pfem;
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

    //! \addtogroup GridFunctionSpace
    //! \{

    //! \brief Indicates lexicographics ordering of the unknowns for composite
    //! grid function spaces.
    //!
    //! this class may be used to pass compile-time
    //! parameters to the implementation of 
    //! \link PowerGridFunctionSpace PowerGridFunctionSpace \endlink or
    //! \link CompositeGridFunctionSpace CompositeGridFunctionSpace \endlink
    struct GridFunctionSpaceLexicographicMapper {};

    //! \brief Indicates using block-wise ordering of the unknowns for composite
    //! grid function spaces.
    //!
    //! The exact blocking structure can be passed as template parameters
    //!
    //! this class may be used to pass compile-time
    //! parameters to the implementation of 
    //! \link PowerGridFunctionSpace PowerGridFunctionSpace \endlink or
    //! \link CompositeGridFunctionSpace CompositeGridFunctionSpace \endlink
    template<int s0 = 1, int s1 = 1, int s2 = 1, int s3 = 1, int s4 = 1, int s5 = 1, int s6 = 1, int s7 = 1, int s8 = 1, int s9 = 1>
    struct GridFunctionSpaceComponentBlockwiseMapper
    {
      static const int size[];
      static const int offset[];
    };
    template<int s0, int s1, int s2, int s3, int s4, int s5, int s6, int s7, int s8, int s9>
    const int GridFunctionSpaceComponentBlockwiseMapper<s0,s1,s2,s3,s4,s5,s6,s7,s8,s9>::
      size[] = { s0, s1, s2, s3, s4, s5, s6, s7, s8, s9 };
    template<int s0, int s1, int s2, int s3, int s4, int s5, int s6, int s7, int s8, int s9>
    const int GridFunctionSpaceComponentBlockwiseMapper<s0,s1,s2,s3,s4,s5,s6,s7,s8,s9>::
      offset[] = { 0, s0, s0+s1, s0+s1+s2, s0+s1+s2+s3, s0+s1+s2+s3+s4,
                   s0+s1+s2+s3+s4+s5, s0+s1+s2+s3+s4+s5+s6, s0+s1+s2+s3+s4+s5+s6+s7,
                   s0+s1+s2+s3+s4+s5+s6+s7+s8, s0+s1+s2+s3+s4+s5+s6+s7+s8+s9 };
    
    //! \brief Indicates using block-wise ordering of the unknowns for composite
    //! grid function spaces.
    //!
    //! this class may be used to pass compile-time
    //! parameters to the implementation of 
    //! \link PowerGridFunctionSpace PowerGridFunctionSpace \endlink or
    //! \link CompositeGridFunctionSpace CompositeGridFunctionSpace \endlink
    struct GridFunctionSpaceBlockwiseMapper : GridFunctionSpaceComponentBlockwiseMapper<> {};

    //! \}

	template<typename T, int k, typename P>
	class PowerGridFunctionSpace;

    //! product of identical grid function spaces
    //! base class that holds implementation of the methods
    //! this is the default version with lexicographic ordering
    //!
    //! \tparam T PLEASE DOCUMENT
    //! \tparam k PLEASE DOCUMENT
    //! \tparam P PLEASE DOCUMENT
	template<typename T, int k, typename P>
	class PowerGridFunctionSpaceBase 
      : public PowerNode<T,k,CountingPointerStoragePolicy>,
        public Countable
	{
	public:
      //! export traits class
      typedef PowerCompositeGridFunctionSpaceTraits<typename T::Traits::GridViewType, 
                                                    typename T::Traits::BackendType,
                                                    P, k>
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
      typedef Dune::PDELab::PowerLocalFunctionSpaceNode<PowerGridFunctionSpaceBase> LocalFunctionSpace;


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

      //! get dimension of root finite element space
	  typename Traits::SizeType globalSize () const
	  {
        // this is bullshit all children may have different
        // size although they have the same type ...
        // [Jö] well, it does happen for the elements I use for the Yee FDTD
        //      scheme, at least
		return offset[k];
	  }

      //! get dimension of this finite element space
      typename Traits::SizeType size () const
      {
        return offset[k];
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

    protected:
      void setup ()
      {
        Dune::dinfo << "PowerGridFunctionSpace(lexicographic version):"
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

    private:
      typename Traits::SizeType childSize[k];
      typename Traits::SizeType offset[k+1];
      typename Traits::SizeType maxlocalsize;
      mutable std::vector<typename Traits::SizeType> childglobal;
	};


    // product of identical grid function spaces
    // base class that holds implementation of the methods
    // specialization for blockwise ordering
	template<typename T, int k, int s>
	class PowerGridFunctionSpaceBase<T,k,GridFunctionSpaceComponentBlockwiseMapper<s> >
      : public PowerNode<T,k,CountingPointerStoragePolicy>,
        public Countable
	{
	public:
      //! export traits class
      typedef PowerCompositeGridFunctionSpaceTraits<typename T::Traits::GridViewType, 
                                                    typename T::Traits::BackendType,
                                                    GridFunctionSpaceComponentBlockwiseMapper<s>, k>
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
      typedef Dune::PDELab::PowerLocalFunctionSpaceNode<PowerGridFunctionSpaceBase> LocalFunctionSpace;
      
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

      //! get dimension of root finite element space
	  typename Traits::SizeType globalSize () const
	  {
        // this is bullshit all children may have different
        // size although they have the same type ...
        // [Jö] well, it does happen for the elements I use for the Yee FDTD
        //      scheme, at least
		return offset[k];
	  }

      //! get dimension of this finite element space
      typename Traits::SizeType size () const
      {
        return offset[k];
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
        return (j%s)+(j/s)*k*s+i*s;
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

    protected:
      void setup ()
      {
        Dune::dinfo << "PowerGridFunctionSpace(blockwise version):"
                    << std::endl;
        Dune::dinfo << "( ";
        offset[0] = 0;
        maxlocalsize = 0;
        for (int i=0; i<k; i++)
        {
          childSize[i] = this->getChild(i).globalSize();
          offset[i+1] = offset[i]+childSize[i];
          Dune::dinfo << childSize[i] << "[" << offset[i] << "] ";
          maxlocalsize += this->getChild(i).maxLocalSize();
        }
        Dune::dinfo << ") total size = " << offset[k]
                    << " max local size = " << maxlocalsize
                    << std::endl;
        /* check the local block size */
        if (this->getChild(0).maxLocalSize()%s != 0)
          DUNE_THROW(Exception,
            "number of DOFs (" << this->getChild(0).maxLocalSize() << ") per component "
            "must be a multiple of the BlockSize (" << s << ")");
        for (int i=1; i<k; i++)
          if (childSize[i]!=childSize[0])
            DUNE_THROW(Exception, "components must be of equal size");
        childglobal.resize(maxlocalsize);
      }

    private:
      typename Traits::SizeType childSize[k];
      typename Traits::SizeType offset[k+1];
      typename Traits::SizeType maxlocalsize;
      mutable std::vector<typename Traits::SizeType> childglobal;
	};

    template<typename T, int k>
    class PowerGridFunctionSpaceBase<T,k,GridFunctionSpaceBlockwiseMapper >
      : public PowerGridFunctionSpaceBase<T,k,GridFunctionSpaceComponentBlockwiseMapper<1> >
    {
    protected:
      using PowerGridFunctionSpaceBase<T,k,GridFunctionSpaceComponentBlockwiseMapper<1> >::setup;
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
	template<typename T, int k, typename P=GridFunctionSpaceLexicographicMapper>
	class PowerGridFunctionSpace 
      : public PowerGridFunctionSpaceBase<T,k,P>
 	{
	public:
      //! export traits class
      typedef PowerCompositeGridFunctionSpaceTraits<typename T::Traits::GridViewType, 
                                                    typename T::Traits::BackendType,
                                                    P, k>
      Traits;

      //! \brief Construct a PowerGridFunction with k clones of the function t
      //! \param t grid function space to clone.
	  PowerGridFunctionSpace (T& t)
      {
        for (int i=0; i<k; i++)
          setChild(i,t);
        PowerGridFunctionSpaceBase<T,k,P>::setup();
      }

      //! \brief Construct a PowerGridFunction with k clones of the function t
      //! \param t grid function space to clone.
	  PowerGridFunctionSpace (T** t)  
      {
         for (int i=0; i<k; i++)
           setChild(i,*(t[i]));
         PowerGridFunctionSpaceBase<T,k,P>::setup();
     }
	};

    //! \}

#ifndef DOXYGEN
	template<typename T, typename P>
	class PowerGridFunctionSpace<T,2,P> 
      : public PowerGridFunctionSpaceBase<T,2,P>
	{
	public:
      //! export traits class
      typedef PowerCompositeGridFunctionSpaceTraits<typename T::Traits::GridViewType, 
                                                    typename T::Traits::BackendType,
                                                    P, 2>
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
                                                    typename T::Traits::BackendType,
                                                    P, 3>
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
                                                    typename T::Traits::BackendType,
                                                    P, 4>
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
                                                    typename T::Traits::BackendType,
                                                    P, 5>
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
                                                    typename T::Traits::BackendType,
                                                    P, 6>
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
                                                    typename T::Traits::BackendType,
                                                    P, 7>
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
                                                    typename T::Traits::BackendType,
                                                    P, 8>
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
                                                    typename T::Traits::BackendType,
                                                    P, 9>
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
                                                    typename T::Traits::BackendType,
                                                    P, 10>
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
#endif // DOXYGEN

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
    

	template<typename P, typename T0, typename T1, typename T2, typename T3,
			 typename T4, typename T5, typename T6, typename T7, typename T8>
	class CompositeGridFunctionSpace;


    // \brief base classe for tuples of grid function spaces
    //
    // base class that holds implementation of the methods
    // this is the default version with lexicographic ordering
    // \tparam P is the ordering parameter. Use e.g. 
    // \link GridFunctionSpaceLexicographicMapper GridFunctionSpaceLexicographicMapper \endlink
    // or \link  GridFunctionSpaceComponentBlockwiseMapper  GridFunctionSpaceComponentBlockwiseMapper \endlink
    // or \link  GridFunctionSpaceBlockwiseMapper  GridFunctionSpaceBlockwiseMapper \endlink
    // \tparam Ti are all grid function spaces
	template<typename P, typename T0, typename T1, typename T2=EmptyChild, typename T3=EmptyChild,
			 typename T4=EmptyChild, typename T5=EmptyChild, typename T6=EmptyChild,
			 typename T7=EmptyChild, typename T8=EmptyChild>
	class CompositeGridFunctionSpaceBase
	  : public CompositeNode<CountingPointerStoragePolicy,T0,T1,T2,T3,T4,T5,T6,T7,T8>,
		public Countable
	{
      typedef CompositeNode<CountingPointerStoragePolicy,T0,T1,T2,T3,T4,T5,T6,T7,T8> BaseT;
	public:
      //! export traits class
      typedef PowerCompositeGridFunctionSpaceTraits<typename T0::Traits::GridViewType, 
                                                    typename T0::Traits::BackendType,
                                                    P,
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
    
    //! \addtogroup GridFunctionSpace
    //! \{

    //! \brief grid function space composed of other grid function spaces
    //!
    //! Composes a tuple of arbitray grid function spaces into a grid function space.
    //! The ordering of the resulting unknowns can be done lexicographically or block-wise.
    //! \tparam P is the ordering parameter. Use e.g. 
    /// \link GridFunctionSpaceLexicographicMapper GridFunctionSpaceLexicographicMapper \endlink
    /// or \link  GridFunctionSpaceBlockwiseMapper  GridFunctionSpaceBlockwiseMapper \endlink
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
      typedef Dune::PDELab::CompositeLocalFunctionSpaceNode<GridFunctionSubSpaceBase> LocalFunctionSpace;

	  // get grid view
	  const typename Traits::GridViewType& gridview () const
	  {
		return pgfs->gridview();
	  }

      //! get dimension of root finite element space
	  typename Traits::SizeType globalSize () const
	  {
        return pgfs->globalSize();
	  }

      //! get dimension of this finite element space
      typename Traits::SizeType size () const
      {
        return pcgfs->globalSize();
      }

	  // get max dimension of shape function space
	  typename Traits::SizeType maxLocalSize () const
	  {
		return pcgfs->maxLocalSize();
	  }

      //! map index from our index set [0,size()-1] to root index set
	  typename Traits::SizeType upMap (typename Traits::SizeType i) const
	  {
		return pgfs->upMap(pgfs->template subMap<k>(i));
	  }

      //! map index from child i's index set into our index set
      template<int i>
	  typename Traits::SizeType subMap (typename Traits::SizeType j) const
	  {
		return pcgfs->template subMap<i>(j);
	  }

    private:
      CountingPointer<GFS const> pgfs;
      CountingPointer<CGFS const> pcgfs;
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
      typedef PowerLocalFunctionSpaceNode<GridFunctionSubSpaceBase> LocalFunctionSpace;

	  // get grid view
	  const typename Traits::GridViewType& gridview () const
	  {
		return pgfs->gridview();
	  }

      //! get dimension of root finite element space
	  typename Traits::SizeType globalSize () const
	  {
        return pgfs->globalSize();
	  }

      //! get dimension of this finite element space
      typename Traits::SizeType size () const
      {
        return pcgfs->globalSize();
      }

	  // get max dimension of shape function space
	  typename Traits::SizeType maxLocalSize () const
	  {
		return pcgfs->maxLocalSize();
	  }

      //! map index from our index set [0,size()-1] to root index set
	  typename Traits::SizeType upMap (typename Traits::SizeType i) const
	  {
		return pgfs->upMap(pgfs->template subMap<k>(i));
	  }

      //! map index from child i's index set into our index set
      template<int i>
	  typename Traits::SizeType subMap (typename Traits::SizeType j) const
	  {
		return pcgfs->template subMap<i>(j);
	  }

    private:
      CountingPointer<GFS const> pgfs;
      CountingPointer<CGFS const> pcgfs;
    };


    // CGFS is a leaf
    template<typename GFS, int k, typename GV, typename FEM, typename CE,
             typename B, typename P>
    class GridFunctionSubSpaceBase<GFS,k, GridFunctionSpace<GV,FEM,CE,B,P> >
      : public Countable // behave like child k of GFS which is a grid function space
    {
      typedef GridFunctionSpace<GV,FEM,CE,B,P> CGFS;

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
      typedef Dune::PDELab::LeafLocalFunctionSpaceNode<GridFunctionSubSpaceBase> LocalFunctionSpace;

	  // get grid view
	  const typename Traits::GridViewType& gridview () const
	  {
		return pcgfs->gridview();
	  }

      // get finite element map
      const FEM& finiteElementMap () const
      {
        return pcgfs->finiteElementMap();
      }

	  // get finite element map
      const FEM& localFiniteElementMap () const DUNE_DEPRECATED
	  {
        return pcgfs->finiteElementMap();
	  }

      //! get dimension of root finite element space
	  typename Traits::SizeType globalSize () const
	  {
        return pgfs->globalSize();
	  }

      //! get dimension of this finite element space
      typename Traits::SizeType size () const
      {
        return pcgfs->globalSize();
      }

	  // get max dimension of shape function space
	  typename Traits::SizeType maxLocalSize () const
	  {
		return pcgfs->maxLocalSize();
	  }

      //! map from our index set [0..size()-1] into root index set
	  typename Traits::SizeType upMap (typename Traits::SizeType i) const
	  {
		return pgfs->upMap(pgfs->template subMap<k>(i));
	  }

	  // compute global indices for one element
      void globalIndices (const typename Traits::FiniteElementType& fe,
                          const Element& e, 
						  std::vector<typename Traits::SizeType>& global) const
	  {
        pcgfs->globalIndices(fe,e,global);
	  }

      // global Indices from element, needs additional finite element lookup
	  void globalIndices (const Element& e,
						  std::vector<typename Traits::SizeType>& global) const
      {
        pcgfs->globalIndices(e,global);
      }

    private:
      CountingPointer<GFS const> pgfs;
      CountingPointer<CGFS const> pcgfs;
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
        Dune::dinfo << "GridFunctionSubSpace:" << std::endl;
        Dune::dinfo << "root space size = " << gfs.globalSize()
                    << " max local size = " << this->maxLocalSize()
                    << std::endl;
      }
    };

  } // namespace PDELab
} // namespace Dune

#endif
