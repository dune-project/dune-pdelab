#ifndef DUNE_PDELAB_FINITELEMENTMAP_HH
#define DUNE_PDELAB_FINITELEMENTMAP_HH

namespace Dune {
  namespace PDELab {

	// collect types exported by a finite element map
	template<class T>
	struct LocalFiniteElementMapTraits
	{
	  //! \brief Type of finite element from local functions
	  typedef T LocalFiniteElementType;
	};

	// interface for a finite element map
	template<class T, class Imp>
	class LocalFiniteElementMapInterface
	{
	public:
	  //! \brief Export traits
	  typedef T Traits;  

	  /** \brief Return local basis for the given entity.

		  The return value is a reference to Traits::LocalBasisType. If 
		  there is a different local basis for two elements then this
		  type must be polymorphic.
	  */ 
	  template<class EntityType>
	  const typename Traits::LocalFiniteElementType& find (const EntityType& e) const
	  {
		return asImp().find(e);
	  }

	private:
	  Imp& asImp () {return static_cast<Imp &> (*this);}
	  const Imp& asImp () const {return static_cast<const Imp &>(*this);}
	};

	// simple implementation where all entities have the same finite element
	template<class Imp>
	class SimpleLocalFiniteElementMap : 
	  public LocalFiniteElementMapInterface<LocalFiniteElementMapTraits<Imp>, 
											SimpleLocalFiniteElementMap<Imp> >
	{
	public:
	  //! \brief export type of the signature
	  typedef LocalFiniteElementMapTraits<Imp> Traits;  

	  //! \brief Use when Imp has a standard constructor
	  SimpleLocalFiniteElementMap ()
	  {}

	  //! \brief Constructor where an instance of Imp can be provided
	  SimpleLocalFiniteElementMap (const Imp& imp_) : imp(imp_)
	  {}

	  //! \brief get local basis functions for entity
	  template<class EntityType>
	  const typename Traits::LocalFiniteElementType& find (const EntityType& e) const
	  {
		return imp;
	  }

	private:
	  Imp imp; // create once
	};

  }
}

#endif
