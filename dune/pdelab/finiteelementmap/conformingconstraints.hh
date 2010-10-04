// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_CONFORMINGCONSTRAINTS_HH
#define DUNE_PDELAB_CONFORMINGCONSTRAINTS_HH

#include <cstddef>

#include<dune/common/exceptions.hh>
#include<dune/grid/common/genericreferenceelements.hh>
#include<dune/grid/common/grid.hh>
#include<dune/common/geometrytype.hh>

#include<dune/pdelab/common/geometrywrapper.hh>
#include <dune/pdelab/finiteelement/traits.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/genericdatahandle.hh>

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

      //! boundary constraints
      /**
       * \tparam F   Boundary grid function returning boundary condition type.
       *             If the value of the boundary grid function evaluated at
       *             the local face center is >0, a dirichlet boundary is
       *             assumed.
       * \tparam IG  intersection geometry
       * \tparam LFS local function space
       * \tparam T   TransformationType
       */
      template<typename F, typename I, typename LFS, typename T>
      void boundary (const F& f, const IntersectionGeometry<I>& ig, 
                     const LFS& lfs, T& trafo) const
      {
        typedef FiniteElementTraits<typename LFS::Traits::FiniteElementType>
          FETraits;

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

        const typename F::Traits::DomainType testpoint = face_refelem.position(0,0);
        f.evaluate(ig,testpoint,bctype);

        for (std::size_t i=0;
             i<FETraits::coefficients(lfs.finiteElement()).size(); i++)
          {
            // The codim to which this dof is attached to
            unsigned int codim =
              FETraits::coefficients(lfs.finiteElement()).localKey(i).codim();

            if (codim==0) continue;

            for (int j=0; j<refelem.size(face,1,codim); j++){
              
              // test point to check whether we have dirichlet or neumann
              const typename F::Traits::DomainType testpoint 
                = face_refelem.position(j,codim-1);
              //      f.evaluate(ig,testpoint,bctype);

              if (bctype > 0 &&
                  static_cast<int>(FETraits::coefficients(lfs.finiteElement()).
                                   localKey(i).subEntity())
                  == refelem.subEntity(face,1,j,codim))
                trafo[i] = empty;
            }
          }
      }
    };

    // extend constraints class by processor boundary 
    class OverlappingConformingDirichletConstraints : public ConformingDirichletConstraints
    {
    public:
      enum { doProcessor = true };
      
      // boundary constraints
      // IG : intersection geometry
      // LFS : local function space
      // T : TransformationType
      template<typename I, typename LFS, typename T>
      void processor (const Dune::PDELab::IntersectionGeometry<I>& ig, 
                      const LFS& lfs, T& trafo) const
      {
        typedef FiniteElementTraits<typename LFS::Traits::FiniteElementType>
          FETraits;

        // determine face
        const int face = ig.indexInInside();

        // find all local indices of this face
        Dune::GeometryType gt = ig.inside()->type();
        typedef typename Dune::PDELab::IntersectionGeometry<I>::ctype DT;
        const int dim = Dune::PDELab::IntersectionGeometry<I>::Entity::Geometry::dimension;


        const Dune::GenericReferenceElement<DT,dim>& refelem = Dune::GenericReferenceElements<DT,dim>::general(gt);

        // empty map means Dirichlet constraint
        typename T::RowType empty;

        // loop over all degrees of freedom and check if it is on given face
        for (size_t i=0; i<FETraits::coefficients(lfs.finiteElement()).size();
             i++)
          {
            // The codim to which this dof is attached to
            unsigned int codim =
              FETraits::coefficients(lfs.finiteElement()).localKey(i).codim();

            if (codim==0) continue;

            for (int j=0; j<refelem.size(face,1,codim); j++)
              if (FETraits::coefficients(lfs.finiteElement()).localKey(i).
                  subEntity() == refelem.subEntity(face,1,j,codim))
                trafo[i] = empty;
          }
      }
    };

    // extend constraints class by processor boundary
    class NonoverlappingConformingDirichletConstraints : public ConformingDirichletConstraints
    {
    public:
      enum { doVolume = true };

      template<typename E, typename LFS, typename T>
      void volume (const Dune::PDELab::ElementGeometry<E>& eg, const LFS& lfs, T& trafo) const
      {
        typedef FiniteElementTraits<typename LFS::Traits::FiniteElementType>
          FETraits;

        // nothing to do for interior entities
        if (eg.entity().partitionType()==Dune::InteriorEntity)
          return;

        // empty map means Dirichlet constraint
        typename T::RowType empty;

		typedef typename LFS::Traits::GridFunctionSpaceType::Traits::BackendType B;

        // loop over all degrees of freedom and check if it is not owned by this processor
        for (size_t i=0; i<FETraits::coefficients(lfs.finiteElement()).size();
             i++)
          {
            if (gh[lfs.globalIndex(i)]!=0)
              {
                trafo[i] = empty;
              }
          }
      }

      template<class GFS>
      void compute_ghosts (const GFS& gfs)
      {
        typedef typename GFS::template VectorContainer<int>::Type V;
        V ighost(gfs);
        Dune::PDELab::GhostDataHandle<GFS,V> gdh(gfs,ighost);
        if (gfs.gridview().comm().size()>1)
          gfs.gridview().communicate(gdh,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);
        ighost.std_copy_to(gh);
        rank = gfs.gridview().comm().rank();
      }

      void print ()
      {
        std::cout << "/" << rank << "/ " << "ghost size=" 
                  << gh.size() << std::endl; 
        for (std::size_t i=0; i<gh.size(); i++)
          std::cout << "/" << rank << "/ " << "ghost[" << i << "]=" 
                    << gh[i] << std::endl;
      } 
      
    private:
      int rank;
      std::vector<int> gh;
    };

  }
}

#endif
