// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_CONFORMINGCONSTRAINTS_HH
#define DUNE_PDELAB_CONFORMINGCONSTRAINTS_HH

#include <cstddef>

#include <dune/common/exceptions.hh>
#include <dune/common/geometrytype.hh>

#include <dune/grid/common/genericreferenceelements.hh>
#include <dune/grid/common/grid.hh>

#include <dune/localfunctions/common/interfaceswitch.hh>

#include <dune/pdelab/common/geometrywrapper.hh>
#include <dune/pdelab/common/typetree.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/localfunctionspacetags.hh>
#include <dune/pdelab/gridfunctionspace/localvector.hh>

namespace Dune {
  namespace PDELab {

    //! Interface for the constraints parameters describing dirichlet constraints
    struct DirichletConstraintsParameters :
      public TypeTree::LeafNode
    {
      template<typename I>
      bool isDirichlet(const I & intersection, const FieldVector<typename I::ctype, I::dimension-1> & coord)
      {
        return true;
      }
    };
    
    //! Dirichlet Constraints construction
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
       * \tparam P   Parameter class, wich fulfills the DirichletConstraintsParameters interface
       * \tparam IG  intersection geometry
       * \tparam LFS local function space
       * \tparam T   TransformationType
       */
      template<typename P, typename IG, typename LFS, typename T>
      void boundary (const P& param, const IG& ig, const LFS& lfs, T& trafo) const
      {
        typedef FiniteElementInterfaceSwitch<
          typename LFS::Traits::FiniteElementType
          > FESwitch;
        typedef FieldVector<typename IG::ctype, IG::dimension-1> FaceCoord;
   
        const int face = ig.indexInInside();

        // find all local indices of this face
        Dune::GeometryType gt = ig.inside()->type();
        typedef typename IG::ctype DT;
        const int dim = IG::Entity::Geometry::dimension;
        const Dune::GenericReferenceElement<DT,dim>& refelem = Dune::GenericReferenceElements<DT,dim>::general(gt);

	    const Dune::GenericReferenceElement<DT,dim-1> & 
	      face_refelem = Dune::GenericReferenceElements<DT,dim-1>::general(ig.geometry().type()); 

        // empty map means Dirichlet constraint
        typename T::RowType empty;

        const FaceCoord testpoint = face_refelem.position(0,0);
        bool isDirichlet = param.isDirichlet(ig,testpoint);

        for (std::size_t i=0;
             i<std::size_t(FESwitch::coefficients(lfs.finiteElement()).size());
             i++)
          {
            // The codim to which this dof is attached to
            unsigned int codim =
              FESwitch::coefficients(lfs.finiteElement()).localKey(i).codim();

            if (codim==0) continue;

            for (int j=0; j<refelem.size(face,1,codim); j++){
              
              // test point to check whether we have dirichlet or neumann
              const FaceCoord testpoint 
                = face_refelem.position(j,codim-1);

              if (isDirichlet &&
                  static_cast<int>(FESwitch::coefficients(lfs.finiteElement()).
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
      template<typename IG, typename LFS, typename T>
      void processor (const IG& ig, const LFS& lfs, T& trafo) const
      {
        typedef FiniteElementInterfaceSwitch<
          typename LFS::Traits::FiniteElementType
          > FESwitch;

        // determine face
        const int face = ig.indexInInside();

        // find all local indices of this face
        Dune::GeometryType gt = ig.inside()->type();
        typedef typename IG::ctype DT;
        const int dim = IG::Entity::Geometry::dimension;

        const Dune::GenericReferenceElement<DT,dim>& refelem = Dune::GenericReferenceElements<DT,dim>::general(gt);

        // empty map means Dirichlet constraint
        typename T::RowType empty;

        // loop over all degrees of freedom and check if it is on given face
        for (size_t i=0; i<FESwitch::coefficients(lfs.finiteElement()).size();
             i++)
          {
            // The codim to which this dof is attached to
            unsigned int codim =
              FESwitch::coefficients(lfs.finiteElement()).localKey(i).codim();

            if (codim==0) continue;

            for (int j=0; j<refelem.size(face,1,codim); j++)
              if (FESwitch::coefficients(lfs.finiteElement()).localKey(i).
                  subEntity() == std::size_t(refelem.subEntity(face,1,j,codim)))
                trafo[i] = empty;
          }
      }
    };

    // extend constraints class by processor boundary
    class NonoverlappingConformingDirichletConstraints : public ConformingDirichletConstraints
    {
    public:
      enum { doVolume = true };

      template<typename EG, typename LFS, typename T>
      void volume (const EG& eg, const LFS& lfs, T& trafo) const
      {
        typedef FiniteElementInterfaceSwitch<
          typename LFS::Traits::FiniteElementType
          > FESwitch;

        // nothing to do for interior entities
        if (eg.entity().partitionType()==Dune::InteriorEntity)
          return;

        // empty map means Dirichlet constraint
        typename T::RowType empty;

		typedef typename LFS::Traits::GridFunctionSpaceType::Traits::BackendType B;

        // loop over all degrees of freedom and check if it is not owned by this processor
        for (size_t i=0; i<FESwitch::coefficients(lfs.finiteElement()).size();
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
        typedef typename GFS::Traits::GridViewType GV;
        typedef typename GV::template Codim<0>::
          template Partition<Interior_Partition>::Iterator Iterator;
        typedef typename GFS::template VectorContainer<int>::Type V;

        gh.assign(gfs.globalSize(), 1);
        V ighost(gfs, 1);
        typename GFS::LocalFunctionSpace lfs(gfs);
        LocalVector<int,AnySpaceTag> lv(gfs.maxLocalSize(), 0);

        const GV &gv = gfs.gridview();
        const Iterator &end = gv.template end<0, Interior_Partition>();
        for(Iterator it = gv.template begin<0, Interior_Partition>();
            it != end; ++it)
        {
          lfs.bind(*it);
          lfs.vwrite(lv, ighost);
        }
        ighost.std_copy_to(gh);

        rank = gv.comm().rank();
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
