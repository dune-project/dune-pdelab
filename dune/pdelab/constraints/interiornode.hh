// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_INTERIORNODECONSTRAINTS_HH
#define DUNE_PDELAB_INTERIORNODECONSTRAINTS_HH

#include <dune/common/array.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/localfunctions/common/interfaceswitch.hh>
#include <dune/localfunctions/common/localkey.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup Constraints
    //! \ingroup FiniteElementMap
    //! \{

    //! \brief constraints all DOFs associated with interior vertices
    //! This allows to implement surface FEM using standard first order FEM
    class InteriorNodeConstraints
    {
      std::vector<bool> interior;
    public:
      enum{doBoundary=false};
      enum{doProcessor=false};
      enum{doSkeleton=false};
      enum{doVolume=true};

      //! volume constraints
      /**
       * \tparam EG  element geometry
       * \tparam LFS local function space
       * \tparam T   TransformationType
       */

      template<typename EG, typename LFS, typename T>
      void volume (const EG& eg, const LFS& lfs, T& trafo) const
      {
        typedef typename EG::Entity Entity;
        typedef typename EG::Geometry Geometry;
        typedef typename Geometry::ctype ctype;
        enum { dim = EG::Entity::dimension, dimw = EG::Entity::dimensionworld };

        // update component
        typename T::RowType empty;
        typedef typename LFS::Traits::SizeType size_type;
        typedef FiniteElementInterfaceSwitch<
          typename LFS::Traits::FiniteElementType
          > FESwitch;
        for (size_type i=0; i<lfs.size(); i++){
          const LocalKey& key = FESwitch::coefficients(lfs.finiteElement()).localKey(i);
          assert(key.codim() == dim && "InteriorNodeConstraints only work for vertex DOFs");
          assert(key.index() == 0   && "InteriorNodeConstraints only work for P1 shape functions");
          // subentity index
          unsigned int local_idx = key.subEntity();
          // global idx

          unsigned int idx = lfs.gridFunctionSpace().gridView().indexSet().subIndex(eg.entity(), local_idx, dim);

          // update constraints
          if (interior[idx])
            trafo[i] = empty;
        }
      }

      const std::vector<bool> & interiorNodes() const
      {
        return interior;
      }

      template<typename GV>
      void updateInteriorNodes(const GV & gv)
      {
        // update vector size
        const int dim = GV::dimension;
        typedef typename GV::Grid::ctype ctype;


        interior.resize(gv.indexSet().size(dim));
        for(int i=0; i< interior.size(); i++)
          interior[i] = true;

        typedef typename GV::template Codim<0>::Iterator ElementIterator;
        typedef typename ElementIterator::Entity Entity;

        // loop over all cells
        for(ElementIterator it = gv.template begin<0>(); it != gv.template end<0>(); ++it)
        {
          const Entity &entity = *it;

          // find boundary faces & associated vertices
          typedef typename GV::IntersectionIterator IntersectionIterator;
          IntersectionIterator iit = gv.ibegin(entity);
          const IntersectionIterator iend = gv.iend(entity);
          for (; iit != iend; ++iit)
          {
            if (iit->boundary())
            {
              // boundary face
              unsigned int f = iit->indexInInside();
              // remember associated vertices
              const ReferenceElement<ctype,dim> & refelem =
                ReferenceElements<ctype,dim>::simplex();
              assert(entity.geometry().type().isSimplex() && "InteriorNodeConstraints only work for simplicial meshes");
              unsigned int sz = refelem.size(f,1, dim);
              assert(sz == dim);
              for (unsigned int v = 0; v < sz; ++v)
              {
                unsigned int local_idx = refelem.subEntity (f,1, v,dim);
                unsigned int idx = gv.indexSet().subIndex(entity, local_idx, dim);
                interior[idx] = false;
              }
            }
          }
        }
      }
    };
    //! \}

  } // end namespace PDELab
} // end namespace Dune

#endif // DUNE_PDELAB_INTERIORNODECONSTRAINTS_HH
