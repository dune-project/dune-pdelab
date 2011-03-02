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

#include <dune/localfunctions/common/interfaceswitch.hh>
#include <dune/localfunctions/common/localkey.hh>

#include "../common/typetree.hh"
#include "../common/geometrywrapper.hh"

#include "../backend/backendselector.hh"

#include"localfunctionspace.hh"
#include"gridfunctionspaceutilities.hh"
#include"powergridfunctionspace.hh"
#include"compositegridfunctionspace.hh"

namespace Dune {
  namespace PDELab {

    //! \addtogroup GridFunctionSpace grid function space
    //! \ingroup PDELab
    //! \{

    //=======================================
    // grid function space : single component case
    //=======================================

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

      //! finite element map
      typedef L FiniteElementMapType;

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

      //! \brief finite element map
      typedef L FiniteElementMapType;

      //! \brief finite element
      typedef typename L::Traits::LocalFiniteElementType FiniteElementType;

	  //! \brief type representing constraints
	  typedef C ConstraintsType;
	};

    class StdVectorBackend;

    //! container construction
    template<typename T, typename E>
    class StdVectorContainer
    {
    public:
      typedef std::vector<E> ContainerType;
      typedef typename ContainerType::iterator iterator;
      typedef typename ContainerType::const_iterator const_iterator;
      typedef E ElementType;
      typedef StdVectorBackend Backend;
      typedef typename std::vector<E>::size_type size_type;
      
      StdVectorContainer (const T& t) : container(t.globalSize()) {}
      StdVectorContainer (const T& t, const E& e) : container(t.globalSize(),e) {}
      StdVectorContainer& operator= (const E& e) // set all elements to same value
      {
        for (typename ContainerType::size_type i=0; i<container.size(); i++)
          container[i] = e;
        return *this;
      }

      ContainerType& base ()
      {
        return container;
      }
      
      const ContainerType& base () const
      {
        return container;
      }

      
      iterator begin()
      {
        return container.begin();
      }
      

      const_iterator begin() const
      {
        return container.begin();
      }

      iterator end()
      {
        return container.end();
      }
      

      const_iterator end() const
      {
        return container.end();
      }


      size_t flatsize() const
      {
        return container.size();
      }

      E& operator[](size_type i)
      {
        return container[i];
      }
      
      const E& operator[](size_type i) const
      {
        return container[i];
      }

      template<typename X>
      void std_copy_to (std::vector<X>& x) const
      {
        typename std::vector<X>::size_type n = flatsize();
        x.resize(n);
        for (size_t i=0; i<n; i++)
          x[i] = container[i];
      }

      template<typename X>
      void std_copy_from (const std::vector<X>& x)
      {
        typename std::vector<X>::size_t n = flatsize();
        x.resize(n);
        for (size_t i=0; i<n; i++)
          container[i] = x[i];
      }

    private:
      ContainerType container;
    };


	//! \brief Simple Backend for std::vector
	class StdVectorBackend
	{      
	public:

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
	  template<typename C, typename E>
	  static const typename StdVectorContainer<C,E>::ContainerType::value_type& 
      access (const StdVectorContainer<C,E>& c, size_type i)
	  {
		return c[i];
	  }
      
	  /** \brief get non const_reference to container element
       *
       *  note: this method does not depend on T!
       */
	  template<typename C, typename E>
	  static typename StdVectorContainer<C,E>::ContainerType::value_type& 
      access (StdVectorContainer<C,E>& c, size_type i)
	  {
		return c[i];
	  }
	};

    template<typename T, typename E>
    struct BackendVectorSelectorHelper<StdVectorBackend,T,E>
    {
      typedef StdVectorContainer<T,E> Type;
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
	class GridFunctionSpace : public TypeTree::LeafNode
	{
	public:
      //! export Traits class
	  typedef GridFunctionSpaceTraits<GV,FEM,CE,B> Traits;
	  typedef typename GV::Traits::template Codim<0>::Entity Element;
	  typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;

 	  //! extract type for storing constraints
	  template<typename E>
	  struct ConstraintsContainer
	  {
		//! \brief define Type as the Type of a container of E's
		typedef ConstraintsTransformation<typename Traits::SizeType,E> Type;
      private:
        ConstraintsContainer () {}
	  };

      typedef LeafGridFunctionSpaceTag ImplementationTag;

	  //! constructor
	  GridFunctionSpace (const GV& gridview, const FEM& fem, const CE& ce_)
		: defaultce(ce_), gv(gridview), pfem(stackobject_to_shared_ptr(fem)), ce(ce_)
	  {
		update();
	  }

	  //! constructor
	  GridFunctionSpace (const GV& gridview, const FEM& fem)
        : gv(gridview), pfem(stackobject_to_shared_ptr(fem)), ce(defaultce)
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
      template<typename StorageIterator>
      void globalIndices (const typename Traits::FiniteElementType& fe,
                          const Element& e, StorageIterator it, StorageIterator endit) const
      {
        typedef FiniteElementInterfaceSwitch<
          typename Traits::FiniteElementType
          > FESwitch;
		// get layout of entity
        const typename FESwitch::Coefficients &coeffs =
          FESwitch::coefficients(fe);

        for (std::size_t i=0; i<std::size_t(coeffs.size()); ++i, ++it)
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
            (*it) = offset[(gtoffset.find(gt)->second)+index]+
              coeffs.localKey(i).index();

            // make sure we don't write past the end of the iterator range
            assert(it != endit);
		  }
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

        typename GTOffsetMap::const_iterator git = gtoffset.find(gt);
        if (git == gtoffset.end())
          return 0;

        typename GV::IndexSet::IndexType index = git->second + gv.indexSet().index(e);
        return offset[index+1]-offset[index];
      }

      //! return vector of global indices associated with the given entity
      template<class EntityType>
      std::size_t dataHandleGlobalIndices (const EntityType& e,
                                           std::vector<typename Traits::SizeType>& global) const
      {
        return dataHandleGlobalIndices(e,global,0,true);
      }

#ifndef DOXYGEN

      template<class EntityType>
      std::size_t dataHandleGlobalIndices (const EntityType& e,
                                           std::vector<typename Traits::SizeType>& global,
                                           std::size_t pos,
                                           bool resize) const
      {
        Dune::GeometryType gt=e.type();

        typename GTOffsetMap::const_iterator git = gtoffset.find(gt);
        if (git == gtoffset.end())
          return 0;

        typename GV::IndexSet::IndexType index = git->second + gv.indexSet().index(e);
        unsigned int n = offset[index+1]-offset[index];
        if (resize)
          global.resize(n+pos);
        for (unsigned i=0; i<n; i++)
          global[pos+i] = offset[index]+i;
        return n;
      }

#endif // DOXYGEN

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
            typedef FiniteElementInterfaceSwitch<
              typename Traits::FiniteElementType
              > FESwitch;
            const typename FESwitch::Coefficients& coeffs =
              FESwitch::coefficients(fe);

			// insert geometry type of all subentities into set
            for (std::size_t i=0; i<std::size_t(coeffs.size()); ++i)
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
            typedef FiniteElementInterfaceSwitch<
              typename Traits::FiniteElementType
              > FESwitch;
            const typename FESwitch::Coefficients& coeffs =
              FESwitch::coefficients(fe);

			// compute maximum number of degrees of freedom per element
            nlocal = std::max(nlocal, static_cast<typename Traits::SizeType>
                                        (coeffs.size()));

			// compute maximum size for each subentity
            for (std::size_t i=0; i<std::size_t(coeffs.size()); ++i)
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
	  shared_ptr<FEM const> pfem;
	  typename Traits::SizeType nlocal;
	  typename Traits::SizeType nglobal;
      const CE& ce;

      typedef std::map<Dune::GeometryType,typename Traits::SizeType> GTOffsetMap;
	  GTOffsetMap gtoffset; // offset in vector for given geometry type
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
	  public TypeTree::LeafNode
	{
	public:
      //! export Traits class
	  typedef GridFunctionSpaceTraits<GV,FEM,CE,B> Traits;
	  typedef typename GV::Traits::template Codim<0>::Entity Element;
	  typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;

 	  //! extract type for storing constraints
	  template<typename E>
	  struct ConstraintsContainer
	  {
		//! \brief define Type as the Type of a container of E's
		typedef ConstraintsTransformation<typename Traits::SizeType,E> Type;
      private:
        ConstraintsContainer () {}
	  };

      typedef LeafGridFunctionSpaceTag ImplementationTag;

	  // constructor
      GridFunctionSpace (const GV& gridview, const FEM& fem, const CE& ce_)
        : gv(gridview), pfem(stackobject_to_shared_ptr(fem)), defaultce(ce_), ce(ce_)
	  {
		update();
	  }

	  // constructor
      GridFunctionSpace (const GV& gridview, const FEM& fem)
        : gv(gridview), pfem(stackobject_to_shared_ptr(fem)), ce(defaultce)
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
      template<typename StorageIterator>
      void globalIndices (const typename Traits::FiniteElementType& fe,
                          const Element& e, StorageIterator it, StorageIterator endit) const
	  {
		// get local coefficients for this entity
        typedef FiniteElementInterfaceSwitch<
          typename Traits::FiniteElementType
          > FESwitch;
        const typename FESwitch::Coefficients &coeffs =
          FESwitch::coefficients(fe);

        for (unsigned int i=0; i<coeffs.size(); ++i, ++it)
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
            (*it) = offset.find(gt)->second+index*dofcountmap.find(gt)->second
              + coeffs.localKey(i).index();

            // make sure we don't write past the end of the iterator range
            assert(it != endit);
		  }
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
        typename DofCountMapType::const_iterator git = dofcountmap.find(gt);
        return git != dofcountmap.end() ? dofcountmap.find(gt)->second : 0;
      }

      //! return vector of global indices associated with the given entity
      template<class EntityType>
      std::size_t dataHandleGlobalIndices (const EntityType& e,
                                           std::vector<typename Traits::SizeType>& global) const
      {
        return dataHandleGlobalIndices(e,global,0,true);
      }

#ifndef DOXYGEN

      template<class EntityType>
      std::size_t dataHandleGlobalIndices (const EntityType& e,
                                           std::vector<typename Traits::SizeType>& global,
                                           std::size_t pos,
                                           bool resize) const
      {
        Dune::GeometryType gt=e.type();
        typename GV::IndexSet::IndexType index = gv.indexSet().index(e);
        typename DofCountMapType::const_iterator git = dofcountmap.find(gt);
        if (git == dofcountmap.end())
          return 0;
        unsigned int n = git->second;
        if (resize)
          global.resize(n);
        for(unsigned i=0; i<n; i++)
          global[pos+i] = offset.find(gt)->second + index*n + i;
        return n;
      }

#endif // DOXYGEN


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
            typedef FiniteElementInterfaceSwitch<
              typename Traits::FiniteElementType
              > FESwitch;
            const typename FESwitch::Coefficients &coeffs =
              FESwitch::coefficients(fe);

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
      shared_ptr<FEM const> pfem;

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
	  public TypeTree::LeafNode
	{
      typedef std::map<unsigned int,unsigned int> DofPerCodimMapType;
	public:
      //! export traits class
      typedef GridFunctionSpaceTraits<GV,FEM,CE,B> Traits;
	  typedef typename GV::Traits::template Codim<0>::Entity Element;
	  typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;

 	  //! extract type for storing constraints
	  template<typename E>
	  struct ConstraintsContainer
	  {
		//! \brief define Type as the Type of a container of E's
		typedef ConstraintsTransformation<typename Traits::SizeType,E> Type;
      private:
        ConstraintsContainer () {}
	  };

      typedef LeafGridFunctionSpaceTag ImplementationTag;

	  // constructors
      GridFunctionSpace (const GV& gridview, const FEM& fem, const IIS& iis_,
                         const CE& ce_)
        : gv(gridview), pfem(stackobject_to_shared_ptr(fem)), iis(iis_), defaultce(ce_), ce(ce_)
	  {
		update();
	  }

      GridFunctionSpace (const GV& gridview, const FEM& fem, const IIS& iis_)
        : gv(gridview), pfem(stackobject_to_shared_ptr(fem)), iis(iis_), ce(defaultce)
	  {
		update();
	  }

      GridFunctionSpace (const GV& gridview, const FEM& fem, const CE& ce_)
        : gv(gridview), pfem(stackobject_to_shared_ptr(fem)), iis(dummyiis), defaultce(ce_), ce(ce_)
	  {
		update();
	  }

      GridFunctionSpace (const GV& gridview, const FEM& fem)
        : gv(gridview), pfem(stackobject_to_shared_ptr(fem)), iis(dummyiis), ce(defaultce)
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
      template<typename StorageIterator>
      void globalIndices (const typename Traits::FiniteElementType& fe,
                          const Element& e, StorageIterator it, StorageIterator endit) const
	  {
        typedef FiniteElementInterfaceSwitch<
          typename Traits::FiniteElementType
          > FESwitch;
		// get local coefficients for this entity
        const typename FESwitch::Coefficients &coeffs =
          FESwitch::coefficients(fe);

        for (unsigned int i=0; i<coeffs.size(); ++i, ++it)
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
            (*it) = offset.find(cd)->second + index * dofpercodim.find(cd)->second
              + coeffs.localKey(i).index();

            // make sure we don't write past the end of the iterator range
            assert(it != endit);
		  }
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
        typename DofPerCodimMapType::const_iterator git = dofpercodim.find(cd);
        return git != dofpercodim.end() ? dofpercodim.find(cd)->second : 0;
      }

      //! return vector of global indices associated with the given entity
      template<class EntityType>
      std::size_t dataHandleGlobalIndices (const EntityType& e,
                                           std::vector<typename Traits::SizeType>& global) const
      {
        return dataHandleGlobalIndices(e,global,0,true);
      }

#ifndef DOXYGEN

      template<class EntityType>
      std::size_t dataHandleGlobalIndices (const EntityType& e,
                                           std::vector<typename Traits::SizeType>& global,
                                           std::size_t pos,
                                           bool resize) const
      {
        const int cd = EntityType::codimension;
        typename DofPerCodimMapType::const_iterator git = dofpercodim.find(cd);
        if (git == dofpercodim.end())
          return 0;
        typename GV::IndexSet::IndexType o = offset.find(cd)->second;
        typename GV::IndexSet::IndexType index = gv.indexSet().index(e);
        unsigned int n = git->second;
        if (resize)
          global.resize(n);
        for (unsigned int i=0; i<n; i++)
          global[pos+i] = o + index*n + i;
//         Dune::dinfo << "[" << gv.grid().comm().rank() << "]: "
//                     << " global indices "
//                     << " offset=" << o
//                     << " index=" << index
//                     << " n=" << n
//                     << std::endl;
        return n;
      }

#endif // DOXYGEN

      //------------------------------

	  // update information, e.g. when grid has changed
	  void update ()
	  {
        Dune::dinfo << "GridFunctionSpace(static size version):" << std::endl;

        // analyse local coefficients of first element

        // check geometry type
        ElementIterator it = gv.template begin<0>();
        const typename Traits::FiniteElementType &fe = pfem->find(*it);
        if (fe.type()!=it->type())
          DUNE_THROW(Exception, "geometry type mismatch in GridFunctionSpace");

        // get local coefficients for this entity
        typedef FiniteElementInterfaceSwitch<
          typename Traits::FiniteElementType
          > FESwitch;
        const typename FESwitch::Coefficients &coeffs =
          FESwitch::coefficients(fe);

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
      shared_ptr<FEM const> pfem;
      const IIS& iis;

	  typename Traits::SizeType nlocal;
	  typename Traits::SizeType nglobal;
      CE defaultce;
      const CE& ce;

      DofPerCodimMapType dofpercodim;
      std::map<unsigned int,typename Traits::SizeType> offset;
	};


    //=======================================
    // Subspace construction
    //=======================================

    template<typename GFS, std::size_t, typename Tag> // primary template, only specializations are used !
    class GridFunctionSubSpaceBase
    {
    };

    template<typename GFS>
    class CompositeGridFunctionSubSpaceNode;


    template<typename Mapper, DUNE_TYPETREE_COMPOSITENODE_TEMPLATE_CHILDREN_FOR_SPECIALIZATION>
    class CompositeGridFunctionSubSpaceNode<CompositeGridFunctionSpace<Mapper,DUNE_TYPETREE_COMPOSITENODE_CHILDTYPES> >
      : public DUNE_TYPETREE_COMPOSITENODE_BASETYPE
    {

      typedef DUNE_TYPETREE_COMPOSITENODE_BASETYPE NodeType;

    public:

      CompositeGridFunctionSubSpaceNode(const typename NodeType::NodeStorage& nodeStorage)
        : NodeType(nodeStorage)
      {}

    };


    // CGFS is a composite
	template<typename GFS, std::size_t k>
    class GridFunctionSubSpaceBase<GFS,k,CompositeGridFunctionSpaceTag>
      : public CompositeGridFunctionSubSpaceNode<typename GFS::template Child<k>::Type>
    {
      typedef typename GFS::template Child<k>::Type CGFS;

    public:
      //! export traits class
	  typedef typename CGFS::Traits Traits;
      typedef CompositeGridFunctionSpaceTag ImplementationTag;

      GridFunctionSubSpaceBase (const GFS& gfs)
        : CompositeGridFunctionSubSpaceNode<CGFS>(gfs.template child<k>().nodeStorage())
        , pgfs(stackobject_to_shared_ptr(gfs))
        , pcgfs(gfs.template childStorage<k>())
      {
      }

  	  //! extract type for storing constraints
	  template<typename E>
	  struct ConstraintsContainer
	  {
		//! \brief define Type as the Type of a container of E's
		typedef ConstraintsTransformation<typename Traits::SizeType,E> Type;
      private:
        ConstraintsContainer () {}
	  };

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
		return pgfs->upMap(pgfs->subMap(k,i));
	  }

      //! map index from child i's index set into our index set
      template<int i>
	  typename Traits::SizeType subMap (typename Traits::SizeType j) const
	  {
		return this->subMap(i,j);
	  }

      typename Traits::SizeType subMap (typename Traits::SizeType i, typename Traits::SizeType j) const
      {
		return pcgfs->subMap(i,j);
      }

    private:
      shared_ptr<GFS const> pgfs;
      shared_ptr<CGFS const> pcgfs;
    };


    // CGFS is a power
	template<typename GFS, std::size_t k>
    class GridFunctionSubSpaceBase<GFS,k,PowerGridFunctionSpaceTag>
      : public TypeTree::PowerNode<typename GFS::template Child<k>::Type::ChildType,GFS::template Child<k>::Type::CHILDREN>
    {
      typedef typename GFS::template Child<k>::Type CGFS;
      typedef TypeTree::PowerNode<typename GFS::template Child<k>::Type::ChildType,GFS::template Child<k>::Type::CHILDREN> NodeType;

    public:
      //! export traits class
	  typedef typename CGFS::Traits Traits;
      typedef PowerGridFunctionSpaceTag ImplementationTag;

      GridFunctionSubSpaceBase (const GFS& gfs)
        : NodeType(gfs.template child<k>().nodeStorage())
        , pgfs(stackobject_to_shared_ptr(gfs))
        , pcgfs(gfs.template childStorage<k>())
      {
      }

 	  //! extract type for storing constraints
	  template<typename E>
	  struct ConstraintsContainer
	  {
		//! \brief define Type as the Type of a container of E's
		typedef ConstraintsTransformation<typename Traits::SizeType,E> Type;
      private:
        ConstraintsContainer () {}
	  };

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
		return pgfs->upMap(pgfs->subMap(k,i));
	  }

      //! map index from child i's index set into our index set
      template<int i>
	  typename Traits::SizeType subMap (typename Traits::SizeType j) const
	  {
		return this->subMap(i,j);
	  }

	  typename Traits::SizeType subMap (typename Traits::SizeType i, typename Traits::SizeType j) const
	  {
		return pcgfs->subMap(i,j);
	  }

    private:
      shared_ptr<GFS const> pgfs;
      shared_ptr<CGFS const> pcgfs;
    };


    // CGFS is a leaf
    template<typename GFS, std::size_t k>
    class GridFunctionSubSpaceBase<GFS,k,LeafGridFunctionSpaceTag>
      : public TypeTree::LeafNode
    {
      typedef typename GFS::template Child<k>::Type CGFS;

    public:
      //! export traits class
	  typedef typename CGFS::Traits Traits;
      typedef typename Traits::GridViewType GV;
      typedef typename Traits::FiniteElementMapType FEM;
	  typedef typename GV::Traits::template Codim<0>::Entity Element;

      typedef LeafGridFunctionSpaceTag ImplementationTag;

      GridFunctionSubSpaceBase (const GFS& gfs)
        : pgfs(stackobject_to_shared_ptr(gfs))
        , pcgfs(gfs.template childStorage<k>())
      {
      }

  	  //! extract type for storing constraints
	  template<typename E>
	  struct ConstraintsContainer
	  {
		//! \brief define Type as the Type of a container of E's
		typedef ConstraintsTransformation<typename Traits::SizeType,E> Type;
      private:
        ConstraintsContainer () {}
	  };

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
		return pgfs->upMap(pgfs->subMap(k,i));
	  }

	  // compute global indices for one element
      template<typename StorageIterator>
      void globalIndices (const typename Traits::FiniteElementType& fe,
                          const Element& e, StorageIterator it, StorageIterator endit) const
	  {
        pcgfs->globalIndices(fe,e,it,endit);
	  }

    private:
      shared_ptr<GFS const> pgfs;
      shared_ptr<CGFS const> pcgfs;
    };

    template<typename GFS, std::size_t k>
    class GridFunctionSubSpace : public GridFunctionSubSpaceBase<GFS,k,typename GFS::template Child<k>::Type::ImplementationTag>
    {

      typedef GridFunctionSubSpaceBase<GFS,k,typename GFS::template Child<k>::Type::ImplementationTag> BaseT;

    public:

      typedef typename Dune::PDELab::TypeTree::TransformTree<GridFunctionSubSpace,gfs_to_lfs>::Type LocalFunctionSpace;

      GridFunctionSubSpace (const GFS& gfs)
        : BaseT(gfs)
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
