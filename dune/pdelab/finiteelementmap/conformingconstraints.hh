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
        typename F::Traits::DomainType ip(0.5); // OK, its 2D here

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
            const Dune::GenericReferenceElement<DT,dim>& refelem = Dune::GenericReferenceElements<DT,dim>::general(gt);
            int face = ig.indexInInside();
            for (int i=0; i<lfs.localFiniteElement().localCoefficients().size(); i++)
              {
                unsigned int codim = lfs.localFiniteElement().localCoefficients().localKey(i).codim();
                if (codim==0) continue;
                for (int j=0; j<refelem.size(face,1,codim); j++)
                    if ((int)lfs.localFiniteElement().localCoefficients().localKey(i).subEntity()==refelem.subEntity(face,1,j,codim))
                    trafo[i] = empty;
              }
          }
      }
    };


  }
}

#endif
