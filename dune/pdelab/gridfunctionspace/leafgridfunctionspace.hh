#ifndef DUNE_PDELAB_LEAFGRIDFUNCTIONSPACE_HH
#define DUNE_PDELAB_LEAFGRIDFUNCTIONSPACE_HH

#include<vector>
#include<set>
#include<map>

#include<dune/grid/common/referenceelements.hh>
#include<dune/grid/common/exception.hh>

#include"../common/countingptr.hh"
#include"../common/multitypetree.hh"
#include"../common/cpstoragepolicy.hh"

namespace Dune {
  namespace PDELab {

	// collect types exported by a leaf grid function space
	template<typename G, class L, typename B>
	struct LeafGridFunctionSpaceTraits
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


	// Simple Backend for std::vector
	class StdVectorBackend
	{
	public:
	  // container construction
	  template<typename T, typename E>
	  class VectorContainer : public std::vector<E>
	  {
		typedef std::vector<E> BaseT;
	  public:
		typedef E ElementType;

		VectorContainer (const T& t) : BaseT(t.globaldimension()) {}
		VectorContainer (const T& t, const E& e) : BaseT(t.globaldimension(),e) {}
		VectorContainer& operator= (const E& e)
		{
		  for (typename BaseT::size_type i=0; i<BaseT::size(); i++)
			(*this)[i] = e;
		  return *this;
		}
	  };

	  // extract type of container element 
	  template<class C>
	  struct Value
	  {
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
	struct LeafGridFunctionGeneralMapper
	{
	  enum {dummy=0} ;
	};


	// mapper for layouts with arbitrary number of entries per entity
	// GV : Type implementing GridView
	// FEM  : Type implementing LocalFiniteElementMapInterface
	// B : Backend type
	// P : Parameter type
	template<typename GV, typename LFEM, typename B=StdVectorBackend, 
			 typename P=LeafGridFunctionGeneralMapper>
	class LeafGridFunctionSpace : public Countable, public LeafNode
	{
	public:
	  typedef LeafGridFunctionSpaceTraits<GV,LFEM,B> Traits;
	  typedef typename GV::Traits::template Codim<0>::Entity Element;
	  typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;

	  // constructor
	  LeafGridFunctionSpace (const GV& gridview, const LFEM& lfem) 
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
// 	  const LFEM& localFiniteElementMap () const
// 	  {
// 		return *plfem;
// 	  }

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
	  void globalIndices (const Element& e, 
						  std::vector<typename Traits::SizeType>& global) const
	  {
		// get layout of entity
// 		LocalLayout layout;
// 		llm.find(entity,layout);
// 		global.resize(layout.size());

// 		for (unsigned int i=0; i<layout.size(); ++i)
// 		  {
// 			// get geometry type of subentity 
// 			Dune::GeometryType gt=Dune::ReferenceElements<double,GV::Grid::dimension>
// 			  ::general(entity.type()).type(layout[i].subentity(),layout[i].codim());

// 			// evaluate consecutive index of subentity
// 			int index = eval_subindex<GV::Grid::dimension>(gv.indexSet(),
// 														   entity,layout[i].subentity(),layout[i].codim());
		
// 			// now compute 
// 			global[i] = offset[ (gtoffset.find(gt)->second)+index ] + layout[i].index();
//		  }
	  }

	private:
	  // update information needed to do mapping
	  void update ()
	  {
		std::cout << "LeafGridFunctionSpace(general version):" << std::endl;

		// determine which geometry types are used
		// needs one traversal of the grid
		typedef std::set<Dune::GeometryType> GtUsedSetType;
		GtUsedSetType gtused;
		for (ElementIterator it = gv.template begin<0>();
			 it!=gv.template end<0>(); ++it)
		  {
			// check geometry type
			if ((plfem->find(*it)).type!=it->type())
			  throw CountableException(counter);

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
				std::cout << "offset["<<index<<"] = " << offset[index] << std::endl;
			  }
		  }

		// now count global number of dofs and compute offset
		nglobal = 0;
		for (typename std::vector<typename Traits::SizeType>::iterator i=offset.begin(); 
			 i!=offset.end(); ++i)
		  {
			std::cout << "nglobal_in=" << nglobal;
			typename Traits::SizeType size = *i;
			std::cout << " size=" << size;
			*i = nglobal;
			nglobal += size;
			std::cout << " nglobal_out=" << nglobal << std::endl;
		  }
		std::cout << "total number of dofs is " << nglobal << std::endl;
	  }

	  const GV& gv;
	  CP<LFEM const> plfem;
	  typename Traits::SizeType nlocal;
	  typename Traits::SizeType nglobal;

	  std::map<Dune::GeometryType,typename Traits::SizeType> gtoffset; // offset in vector for given geometry type
	  std::vector<typename Traits::SizeType> offset; // offset into big vector for each entity;
	};



	
  }
}

#endif
