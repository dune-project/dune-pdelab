// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_FINITELEMENTMAP_HH
#define DUNE_PDELAB_FINITELEMENTMAP_HH

#include <dune/grid/common/genericreferenceelements.hh>

#include"../common/countingptr.hh"

namespace Dune {
  namespace PDELab {

    //! \addtogroup FiniteElementMap
    //! \ingroup PDELab
    //! \{

	//! collect types exported by a finite element map
	template<class T>
	struct LocalFiniteElementMapTraits
	{
	  //! \brief Type of finite element from local functions
	  typedef T LocalFiniteElementType;
	};

	//! interface for a finite element map
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

	//! simple implementation where all entities have the same finite element
	template<class Imp>
	class SimpleLocalFiniteElementMap : 
	  public LocalFiniteElementMapInterface<LocalFiniteElementMapTraits<Imp>, 
											SimpleLocalFiniteElementMap<Imp> >,
	  public Countable
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

    /** \brief implementation for finite elements requiring oriented edges
     *
     * This is for edge elements. It works for one type of Geometry only, and
     * the requirements fot the local finite element are:
     *
     *  - The default orientation of the shape functions on the edges must be
     *    from the vertex with the lower index within the reference element to
     *    the vertex with the higher index.
     *  - The local finite element must allow assignment.
     *  - The local finite element must have a default constructor.
     *  - the local finite element hust have a constructor of the form
     *    FE(unsigned int orientation), where orientation is a bitfield.  If
     *    bit i in orientation is set it means that the orientation of the
     *    shape function corresponding to the edge with id i in the reference
     *    element is inverted.
     *
     * \tparam GV  Type of gridview to work with
     * \tparam FE  Type of local finite element
     * \tparam Imp Type of the final LocalFiniteElementMap implementation
     */
    template<typename GV, typename FE, typename Imp>
    class EdgeS0LocalFiniteElementMap
      : public LocalFiniteElementMapInterface<
          LocalFiniteElementMapTraits<FE>,
          Imp
        >
      , public Countable
    {
      typedef typename GV::IndexSet IndexSet;
      static const int dim = GV::dimension;
     
    public:
	  //! \brief export type of the signature
	  typedef LocalFiniteElementMapTraits<FE> Traits;  

	  //! \brief construct EdgeSLocalFiniteElementMap
	  EdgeS0LocalFiniteElementMap (const GV& gv_) 
        : gv(gv_), is(gv_.indexSet()), orient(gv_.size(0))
	  {
        typedef typename GV::Grid::ctype ct;
        const GenericReferenceElement<ct, dim> &refElem =
          GenericReferenceElements<ct, dim>::general(FE().type());

        // create all variants
        variant.resize(1 << refElem.size(dim-1));
        for (unsigned int i=0; i<variant.size(); i++)
          variant[i] = FE(i);

        // compute orientation for all elements
        typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
        typedef typename GV::IntersectionIterator IntersectionIterator;

        // loop once over the grid
        for (ElementIterator it = gv.template begin<0>(); it!=gv.template end<0>(); ++it)
          {
            unsigned int elemid = is.template index<0>(*it);
            orient[elemid] = 0;

            std::vector<unsigned int> vid(refElem.size(dim));
            for(unsigned int i = 0; i < vid.size(); ++i)
              vid[i] = is.subIndex(*it, i, dim);

            // loop over all edges of the element
            for(int i = 0; i < refElem.size(dim-1); ++i) {
              int v0 = refElem.subEntity(i, dim-1, 0, dim);
              int v1 = refElem.subEntity(i, dim-1, 1, dim);
              // if (edge orientation in refelement) != (edge orientation in indexset)
              if((v0 > v1) != (vid[v0] > vid[v1]))
                orient[elemid] |= 1 << i;
            }
          }
      }

	  //! \brief get local basis functions for entity
	  template<class EntityType>
	  const typename Traits::LocalFiniteElementType& find (const EntityType& e) const
	  {
        return variant[orient[is.template index<0>(e)]];
	  }

	private:
      const GV& gv;
      std::vector<FE> variant;
      const IndexSet& is;
      std::vector<unsigned char> orient;
    };

    //! \} group FiniteElementMap

  } // namespace PDELab
} // namespace Dune

#endif
