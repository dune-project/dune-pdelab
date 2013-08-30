// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_MFDCONSTRAINTS_HH
#define DUNE_PDELAB_MFDCONSTRAINTS_HH

#include<dune/geometry/referenceelements.hh>
#include<dune/geometry/type.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup Constraints
    //! \ingroup FiniteElementMap
    //! \{

    //! Mimetic Constraints
    class MimeticConstraints
    {
    public:
      enum{doVolume=false};
      enum{doSkeleton=false};
      enum{doBoundary=true};
      enum{doProcessor=false};

      //! boundary constraints
      /**
       * \tparam F   grid function returning boundary condition type
       * \tparam IG  intersection geometry
       * \tparam LFS local function space
       * \tparam T   TransformationType
       */
      template<typename B, typename IG, typename LFS, typename T>
      void boundary(const B& b, const IG& ig, const LFS& lfs, T& trafo) const
      {
        static const unsigned int dimIntersection = IG::dimension - 1;
        typedef typename IG::ctype ctype;

        GeometryType gt = ig.intersection().type();

        Dune::FieldVector<ctype,dimIntersection> center
          = ReferenceElements<ctype,dimIntersection>::general(gt).position(0,0);
        if(b.isDirichlet(ig, center) )
          {
            typename T::RowType empty;
            trafo[ig.intersectionIndex()] = empty;
          }
      }
    };
    //! \}

  }
}

#endif
