// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_Pk2DFEM_HH
#define DUNE_PDELAB_Pk2DFEM_HH

#include<dune/common/exceptions.hh>
#include<dune/grid/common/referenceelements.hh>
#include<dune/common/geometrytype.hh>
#include<dune/grid/common/quadraturerules.hh>

#include<dune/finiteelements/pk2d.hh>
#include"finiteelementmap.hh"

namespace Dune {
  namespace PDELab {

	//! wrap up element from local functions
    //! \ingroup FiniteElementMap

    template<typename GV, typename D, typename R, unsigned int k>
    class Pk2DLocalFiniteElementMap : 
	  public LocalFiniteElementMapInterface<LocalFiniteElementMapTraits< Dune::Pk2DLocalFiniteElement<D,R,k> >, 
											Pk2DLocalFiniteElementMap<GV,D,R,k> >,
      public Countable
    {
      typedef Dune::Pk2DLocalFiniteElement<D,R,k> FE;
      typedef typename GV::IndexSet IndexSet;
     
    public:
	  //! \brief export type of the signature
	  typedef LocalFiniteElementMapTraits<FE> Traits;  

	  //! \brief Use when Imp has a standard constructor
	  Pk2DLocalFiniteElementMap (const GV& gv_) : is(gv_.indexSet())
	  {
        // create all variants 
        for (int i=0; i<8; i++)
          variant[i] = FE(i);
      }

	  //! \brief get local basis functions for entity
	  template<class EntityType>
	  const typename Traits::LocalFiniteElementType& find (const EntityType& e) const
	  {
        unsigned int n0,n1,n2;
        n0 = is.template subIndex<2>(e,0);
        n1 = is.template subIndex<2>(e,1);
        n2 = is.template subIndex<2>(e,2);
        unsigned int j=0;
        if (n1>n2) j += 1;
        if (n0>n2) j += 2;
        if (n0>n1) j += 4;
		return variant[j];
	  }

	private:
      FE variant[8];
      const IndexSet& is;
    };


    //! Constraints construction for P1 elements on triangles
    class Pk2DConstraints
    {
    public:
      enum { doBoundary = true };
      enum { doSkeleton = false };
      enum { doVolume = false };

      // boundary constraints
      // F : grid function returning boundary condition type
      // IG : intersection geometry
      // LFS : local function space
      // T : TransformationType
      template<typename F, typename I, typename LFS, typename T>
      void boundary (const F& f, const IntersectionGeometry<I>& ig, 
                     const LFS& lfs, T& trafo) const
      {
        // 2D here, get midpoint of edge
        typename F::Traits::DomainType ip(0.5);

        // determine type of boundary condition FOR WHOLE FACE IN THE MIDPOINT
        typename F::Traits::RangeType bctype;
        f.evaluate(ig,ip,bctype);

        // if dirichlet boundary constrain all dofs on that face
        if (bctype>0)
          {
            // empty map indicates Dirichlet constraint
            typename T::RowType empty;

            // find all local indices of this face
            Dune::GeometryType gt = ig.inside()->type();
            typedef typename IntersectionGeometry<I>::ctype DT;
            const int dim = IntersectionGeometry<I>::Entity::Geometry::dimension;
            const int dimw = IntersectionGeometry<I>::Entity::Geometry::dimensionworld;
            const Dune::ReferenceElement<DT,dim>& refelem = Dune::ReferenceElements<DT,dim>::general(gt);
            int face = ig.numberInSelf();
            for (size_t i=0; i<lfs.localFiniteElement().localCoefficients().size(); i++)
              {
                unsigned int codim = lfs.localFiniteElement().localCoefficients().localIndex(i).codim();
                if (codim==0) continue;
                for (int j=0; j<refelem.size(face,1,codim); j++)
                  if (lfs.localFiniteElement().localCoefficients().localIndex(i).subentity()==refelem.subEntity(face,1,j,codim))
                    trafo[i] = empty;
              }
          }
      }
    };


  }
}

#endif
