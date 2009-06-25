// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_HANGINGNODECONSTRAINTS_HH
#define DUNE_PDELAB_HANGINGNODECONSTRAINTS_HH

#include<dune/common/exceptions.hh>
#include<dune/grid/common/genericreferenceelements.hh>
#include<dune/common/geometrytype.hh>
#include<dune/pdelab/common/geometrywrapper.hh>
#include<dune/pdelab/finiteelementmap/conformingconstraints.hh>
#include<dune/pdelab/finiteelementmap/hangingnodemanager.hh>

namespace Dune {
  namespace PDELab {

    class HangingNodesConstraintsAssemblers
    {
    public:
      class CubeGridQ1Assembler
      {
      public:
        template<typename I, typename LFS, typename T, typename FlagVector>
        static void assembleConstraints(const FlagVector & nodeState_e, const FlagVector & nodeState_f,
                                        const bool & e_has_hangingnodes, const bool & f_has_hangingnodes,
                                        const LFS & lfs_e, const LFS & lfs_f,
                                        T& trafo_e, T& trafo_f,
                                        const IntersectionGeometry<I>& ig) 
        {
          typedef IntersectionGeometry<I> Intersection;
          typedef typename Intersection::EntityPointer CellEntityPointer;
          typedef typename Intersection::Geometry FaceGeometry;
          typedef typename FaceGeometry::ctype DT;
          typedef typename LFS::Traits::LocalFiniteElementType LocalFiniteElementType; 
          typedef typename LocalFiniteElementType::Traits::LocalCoefficientsType LocalCoefficientType;
          typedef typename LocalFiniteElementType::Traits::LocalBasisType::Traits::DomainFieldType DFT;
          typedef typename LocalFiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType RFT;
          typedef typename LFS::Traits::SizeType SizeType;

          const CellEntityPointer e = ig.inside();
          const CellEntityPointer f = ! ig.boundary() ? ig.outside() : ig.inside();

          const int dimension = Intersection::dimension;

          typedef Dune::GenericReferenceElement<DT,dimension> GRE;
          const GRE& refelement_e = Dune::GenericReferenceElements<DT,dimension>::general(e->type());
          const GRE& refelement_f = Dune::GenericReferenceElements<DT,dimension>::general(f->type());

          // If both entities have hangingnodes, then the face is
          // conforming and no constraints have to be applied.
          if(e_has_hangingnodes && f_has_hangingnodes)
            return;

          // Choose local function space etc for element with hanging nodes
          const LFS & lfs = e_has_hangingnodes ? lfs_e : lfs_f;
          const LocalCoefficientType & localCoefficients = lfs.localFiniteElement().localCoefficients();
          assert(localCoefficients.size()==8);
          const int faceindex = e_has_hangingnodes ? ig.indexInInside() : ig.indexInOutside();
          const GRE & refelement = e_has_hangingnodes ? refelement_e : refelement_f;
          const FlagVector & nodeState = e_has_hangingnodes ? nodeState_e : nodeState_f;
          T & trafo = e_has_hangingnodes ? trafo_e : trafo_f;

          // A map mapping the local entity index to the local coefficient index
          std::vector<SizeType> mapEntityCoeff(localCoefficients.size());
          for (SizeType i=0; i<localCoefficients.size(); i++){
            mapEntityCoeff[localCoefficients.localKey(i).subEntity()] = localCoefficients.localKey(i).index() + i;
            if( localCoefficients.localKey(i).codim() != dimension){
              DUNE_THROW(Dune::InvalidStateException,
                         "Local coefficients are expected to be bound to vertex entities");
            }
          }//i

          // A map mapping the local indices from the face to local
          // indices of the cell
          std::vector<int> m(refelement.size(faceindex,1,dimension));
          for (int j=0; j<refelement.size(faceindex,1,dimension); j++)
            m[j] = refelement.subEntity(faceindex,1,j,dimension);
            
          // Find the corresponding entity in the reference element
          for (int j=0; j<refelement.size(faceindex,1,dimension); j++){

            // The contribution factors of the base function bound to entity j
            typename T::RowType contribution;

            if(dimension == 3){
              assert(nodeState.size() == 8);
              
              const SizeType i = 4*j;

              // Neigbor relations in local indices on a quadrilateral face: 
              // {Node, Direct Neighbor, Direct Neighbor, Diagonal Neighbor, Node, ...}
              const unsigned int fi[16] = {0,1,2,3, 1,0,3,2, 2,0,3,1, 3,1,2,0};

              // Only hanging nodes have contribution to other nodes
              if(nodeState[m[j]].isHanging() && !nodeState[m[j]].isBoundary()){

                const SizeType node_coeff_index = mapEntityCoeff[m[j]];
                
                // If both neighbors are hanging nodes, then this node
                // is diagonal to the target of the contribution
                if(nodeState[m[fi[i+1]]].isHanging() && nodeState[m[fi[i+2]]].isHanging()){
                  //if(!nodeState[m[fi[i+3]]].isBoundary())
                    contribution[mapEntityCoeff[m[fi[i+3]]]] = 0.25;
                    trafo[node_coeff_index] = contribution;
                }
                // Direct neigbor
                else if(!nodeState[m[fi[i+1]]].isHanging()){
                  //if(!nodeState[m[fi[i+1]]].isBoundary())
                    contribution[mapEntityCoeff[m[fi[i+1]]]] = 0.5;
                    trafo[node_coeff_index] = contribution;
                }
                // Direct neigbor
                else if(!nodeState[m[fi[i+2]]].isHanging()){
                  //if(!nodeState[m[fi[i+2]]].isBoundary())
                    contribution[mapEntityCoeff[m[fi[i+2]]]] = 0.5;
                    trafo[node_coeff_index] = contribution;
                }
                // Write into local constraints container
              }

            } // dimension == 3
            else if(dimension == 2){
              
            }

          }//j          

        }// assembleConstraints
      };
    };

    //! Constraints construction
    // works in any dimension and on all element types
    template <class Grid, class HangingNodesConstraintsAssemblerType, class BoundaryFunction>
    class HangingNodesDirichletConstraints : public ConformingDirichletConstraints
    {
    private:
      typedef Dune::PDELab::HangingNodeManager<Grid,BoundaryFunction> HangingNodeManager;
      HangingNodeManager manager;

    public:
      enum { doBoundary = true };
      enum { doSkeleton = true };
      enum { doVolume = false };
      enum { dimension = Grid::dimension };

      HangingNodesDirichletConstraints(Grid & grid, bool adaptToIsolatedHangingNotes, const BoundaryFunction & _boundaryFunction)
        : manager(grid, _boundaryFunction)
      {
        if(adaptToIsolatedHangingNotes)
          manager.adaptToIsolatedHangingNodes();
      }


      // boundary constraints
      // F : grid function returning boundary condition type
      // IG : intersection geometry
      // LFS : local function space
      // T : TransformationType
      template<typename F, typename I, typename LFS, typename T>
      void boundary (const F& f, const IntersectionGeometry<I>& ig, 
                     const LFS& lfs, T& trafo) const
      {
        ConformingDirichletConstraints::boundary(f,ig,lfs,trafo);
        return ;

        typedef IntersectionGeometry<I> Intersection;
        typedef typename Intersection::EntityPointer CellEntityPointer;
        typedef typename Intersection::Geometry FaceGeometry;
        typedef typename FaceGeometry::ctype DT;
        
        const CellEntityPointer e = ig.inside();

        const Dune::GenericReferenceElement<DT,dimension>& refelem_e
          = Dune::GenericReferenceElements<DT,dimension>::general(e->type());

        // the return values of the hanging node manager
        typedef typename std::vector<typename HangingNodeManager::NodeState> FlagVector;
        const FlagVector isHangingNode_e(manager.hangingNodes(e));

        // just to make sure that the hanging node manager is doing
        // what is expected of him
        assert(refelem_e.size(dimension) == isHangingNode_e.size());

        // the LOCAL indices of the faces in the reference element
        const int faceindex_e = ig.indexInInside();
        
        typedef typename LFS::Traits::LocalFiniteElementType LocalFiniteElementType; 
        typedef typename LocalFiniteElementType::Traits::LocalCoefficientsType LocalCoefficientType;
        typedef typename LocalFiniteElementType::Traits::LocalBasisType::Traits::DomainFieldType DFT;
        typedef typename LFS::Traits::SizeType SizeType;

        bool e_has_hangingnodes = false;
        {
          for (int j=0; j<refelem_e.size(faceindex_e,1,dimension); j++){
            const int index = refelem_e.subEntity(faceindex_e,1,j,dimension);
            if(isHangingNode_e[index].isHanging())
              e_has_hangingnodes = true;
          }
        }
        bool f_has_hangingnodes = false;

        if(! e_has_hangingnodes)
          return;

        HangingNodesConstraintsAssemblerType::
          assembleConstraints(isHangingNode_e, isHangingNode_e,
                              e_has_hangingnodes, f_has_hangingnodes,
                              lfs,lfs,
                              trafo, trafo,
                              ig);
        
      }

      // boundary constraints
      // I : intersection 
      // LFS : local function space
      // T : TransformationType
      template<typename I, typename LFS, typename T>
      void skeleton (const IntersectionGeometry<I>& ig, 
                     const LFS& lfs_e, const LFS& lfs_f, 
                     T& trafo_e, T& trafo_f) const
      {
        typedef IntersectionGeometry<I> Intersection;
        typedef typename Intersection::EntityPointer CellEntityPointer;
        typedef typename Intersection::Geometry FaceGeometry;
        typedef typename FaceGeometry::ctype DT;
        
        const CellEntityPointer e = ig.inside();
        const CellEntityPointer f = ig.outside();

        const Dune::GenericReferenceElement<DT,dimension>& refelem_e
          = Dune::GenericReferenceElements<DT,dimension>::general(e->type());
        const Dune::GenericReferenceElement<DT,dimension>& refelem_f
          = Dune::GenericReferenceElements<DT,dimension>::general(f->type());

        // the return values of the hanging node manager
        typedef typename std::vector<typename HangingNodeManager::NodeState> FlagVector;
        const FlagVector isHangingNode_e(manager.hangingNodes(e));
        const FlagVector isHangingNode_f(manager.hangingNodes(f));

        // just to make sure that the hanging node manager is doing
        // what is expected of him
        assert(refelem_e.size(dimension) == isHangingNode_e.size());
        assert(refelem_f.size(dimension) == isHangingNode_f.size());

        // the LOCAL indices of the faces in the reference element
        const int faceindex_e = ig.indexInInside();
        const int faceindex_f = ig.indexInOutside();
        
        typedef typename LFS::Traits::LocalFiniteElementType LocalFiniteElementType; 
        typedef typename LocalFiniteElementType::Traits::LocalCoefficientsType LocalCoefficientType;
        typedef typename LocalFiniteElementType::Traits::LocalBasisType::Traits::DomainFieldType DFT;
        typedef typename LFS::Traits::SizeType SizeType;

        bool e_has_hangingnodes = false;
        {
          for (int j=0; j<refelem_e.size(faceindex_e,1,dimension); j++){
            const int index = refelem_e.subEntity(faceindex_e,1,j,dimension);
            if(isHangingNode_e[index].isHanging())
              e_has_hangingnodes = true;
          }
        }
        bool f_has_hangingnodes = false;
        {
          for (int j=0; j<refelem_f.size(faceindex_f,1,dimension); j++){
            const int index = refelem_f.subEntity(faceindex_f,1,j,dimension);
            if(isHangingNode_f[index].isHanging())
              f_has_hangingnodes = true;
          }
        }

        if(! e_has_hangingnodes && ! f_has_hangingnodes)
          return;

        HangingNodesConstraintsAssemblerType::
          assembleConstraints(isHangingNode_e, isHangingNode_f,
                              e_has_hangingnodes, f_has_hangingnodes,
                              lfs_e,lfs_f,
                              trafo_e, trafo_f,
                              ig);
      } // skeleton

    };


  }
}

#endif
