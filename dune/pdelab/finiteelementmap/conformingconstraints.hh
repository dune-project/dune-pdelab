// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_CONFORMINGCONSTRAINTS_HH
#define DUNE_PDELAB_CONFORMINGCONSTRAINTS_HH

#include<dune/common/exceptions.hh>
#include<dune/grid/common/genericreferenceelements.hh>
#include<dune/common/geometrytype.hh>
#include<dune/pdelab/common/geometrywrapper.hh>

namespace Dune {
  namespace PDELab {

    //! Constraints construction
    // works in any dimension and on all element types
    class ConformingDirichletConstraints
    {
    public:
      enum { doBoundary = true };
      enum { doProcessor = false };
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
        typename F::Traits::RangeType bctype;

        const int face = ig.indexInInside();

        // find all local indices of this face
        Dune::GeometryType gt = ig.inside()->type();
        typedef typename IntersectionGeometry<I>::ctype DT;
        const int dim = IntersectionGeometry<I>::Entity::Geometry::dimension;
        const Dune::GenericReferenceElement<DT,dim>& refelem = Dune::GenericReferenceElements<DT,dim>::general(gt);

	    const Dune::GenericReferenceElement<DT,dim-1> & 
	      face_refelem = Dune::GenericReferenceElements<DT,dim-1>::general(ig.geometry().type()); 

        // empty map means Dirichlet constraint
        typename T::RowType empty;

        for (int i=0; i<lfs.localFiniteElement().localCoefficients().size(); i++)
          {
            // The codim to which this dof is attached to
            unsigned int codim = lfs.localFiniteElement().localCoefficients().localKey(i).codim();

            if (codim==0) continue;

            for (int j=0; j<refelem.size(face,1,codim); j++){
              
              // test point to check whether we have dirichlet or neumann
              const typename F::Traits::DomainType testpoint 
                = face_refelem.position(j,codim-1);
              f.evaluate(ig,testpoint,bctype);

              if (bctype > 0 && (int) lfs.localFiniteElement().localCoefficients().localKey(i).subEntity()
                  ==
                  refelem.subEntity(face,1,j,codim))
                trafo[i] = empty;
            }
          }

      }
    };


  }
}

#endif
