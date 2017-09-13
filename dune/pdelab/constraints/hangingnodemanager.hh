// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_CONSTRAINTS_HANGINGNODEMANAGER_HH
#define DUNE_PDELAB_CONSTRAINTS_HANGINGNODEMANAGER_HH

#include<dune/grid/common/grid.hh>
#include<dune/grid/common/mcmgmapper.hh>
#include<dune/common/float_cmp.hh>

#include"../common/geometrywrapper.hh"

namespace Dune {
  namespace PDELab {

    template<class Grid, class BoundaryFunction>
    class HangingNodeManager
    {
    private:
#ifdef DEBUG
      enum{ verbosity = 2 };
#else
      enum{ verbosity = 0 };
#endif
      typedef typename Grid::LeafIndexSet::IndexType IndexType;

    private:
      class NodeInfo
      {
      public:
        // Minimum level of elements containing this node
        unsigned short minimum_level;

        // Maximum level of elements containing this node
        unsigned short maximum_level;

        // Minimum level of elements touching this node
        unsigned short minimum_touching_level;

        bool is_boundary;

        void addLevel(unsigned short level)
        {
          minimum_level = std::min(minimum_level,level);
          maximum_level = std::max(maximum_level,level);
        }

        void addTouchingLevel(unsigned short level)
        {
          minimum_touching_level = std::min(minimum_touching_level,level);
        }

        NodeInfo() : minimum_level(std::numeric_limits<unsigned short>::max()),
                     maximum_level(0),
                     minimum_touching_level(std::numeric_limits<unsigned short>::max()),
                     is_boundary(false)
        {}
      };

      std::vector<NodeInfo> node_info;

    public:

      class NodeState
      {
        friend class HangingNodeManager;
      private:
        bool is_boundary;
        bool is_hanging;

      public:

        inline bool isBoundary() const { return is_boundary; }
        inline bool isHanging() const { return is_hanging; }

        NodeState() :is_boundary(false),
                     is_hanging(false)
        {}
      };


      typedef typename Grid::LeafGridView GridView;
      enum {dim = GridView::dimension};
      typedef typename GridView::template Codim<0>::Entity Cell;

      typedef typename GridView::template Codim<0>::Iterator Iterator;
      typedef typename GridView::IntersectionIterator IntersectionIterator;
      typedef typename GridView::Grid::ctype ctype;
      typedef typename Dune::FieldVector<ctype,dim> Point;
      typedef typename Dune::FieldVector<ctype,dim-1> FacePoint;

      typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView,
                                                        MCMGElementLayout> CellMapper;

      Grid & grid;
      const BoundaryFunction & boundaryFunction;
      CellMapper cell_mapper;

    public:

      void analyzeView()
      {
        cell_mapper.update();
        const typename GridView::IndexSet& indexSet = grid.leafGridView().indexSet();

        node_info = std::vector<NodeInfo>(indexSet.size(dim));

        GridView gv = grid.leafGridView();

        // loop over all codim<0> leaf elements of the partially refined grid
        for(const auto& cell : elements(gv)) {

          auto reference_element = referenceElement(cell.geometry());

          // level of this element
          const unsigned short level = cell.level();

          // number of vertices in this element
          const IndexType v_size = reference_element.size(dim);

          // update minimum_level and maximum_level for vertices in this
          // cell
          // loop over all vertices of the element
          for(IndexType i=0; i<v_size; ++i){
            const IndexType v_globalindex = indexSet.subIndex(cell,i,dim);
            NodeInfo& ni = node_info[v_globalindex];
            ni.addLevel(level);

            if(verbosity>10){
              // This will produce a lot of output on the screen!
              std::cout << "   cell-id=" << cell_mapper.index(cell);
              std::cout << "   level=" << level;
              std::cout << "   v_size=" << v_size;
              std::cout << "   v_globalindex = " << v_globalindex;
              std::cout << "   maximum_level = " << ni.maximum_level;
              std::cout << "   minimum_level = " << ni.minimum_level;
              std::cout << std::endl;
            }

          }

          // Now we still have to update minimum_touching_level for this
          // cell

          typedef typename IntersectionIterator::Intersection Intersection;

          // Loop over faces
          unsigned int intersection_index = 0;
          for(const auto& intersection : intersections(gv,cell)) {
            ++intersection_index;

            auto reference_face_element = referenceElement(intersection.geometry());

            const int eLocalIndex =  intersection.indexInInside();
            const int e_level = intersection.inside().level();

            // number of vertices on the face
            const int e_v_size = reference_element.size(eLocalIndex,1,dim);

            if(intersection.boundary()) {

              // loop over vertices on the face
              for(int i=0; i<e_v_size;++i){
                const int e_v_index = reference_element.subEntity(eLocalIndex,1,i,dim);
                const IndexType v_globalindex = indexSet.subIndex(cell,e_v_index,dim);

                const FacePoint facelocal_position = reference_face_element.position(i,dim-1);

                /*
                  typename BoundaryFunction::Traits::RangeType boundary_value;
                  boundaryFunction.evaluate(IntersectionGeometry<Intersection>(*fit,intersection_index),
                  facelocal_position,boundary_value);
                  if(boundary_value)
                  node_info[v_globalindex].is_boundary=true;
                  else
                  node_info[v_globalindex].is_boundary=false;
                */

                NodeInfo& ni = node_info[v_globalindex];
                ni.is_boundary = boundaryFunction.isDirichlet(IntersectionGeometry<Intersection>(intersection,intersection_index),facelocal_position);
                ni.addTouchingLevel(e_level);
              }

              // We are done here - the remaining tests are only required for neighbor intersections
              continue;
            }

            const int f_level = intersection.outside().level();

            // a conforming face has no hanging nodes
            if(intersection.conforming())
              continue;

            // so far no support for initially non conforming grids
            assert(e_level != f_level);

            // this check needs to be performed on the element containing
            // the hanging node only
            if(e_level < f_level)
              continue;

            // loop over vertices on the face
            for(int i=0; i<e_v_size;++i){
              const int e_v_index = reference_element.subEntity(eLocalIndex,1,i,dim);
              const IndexType v_globalindex = indexSet.subIndex(cell,e_v_index,dim);

              node_info[v_globalindex].addTouchingLevel(f_level);
            }

          } // end of loop over faces

        } // end loop over codim<0> leaf elements
      }

      HangingNodeManager(Grid & _grid, const BoundaryFunction & _boundaryFunction)
        : grid(_grid),
          boundaryFunction(_boundaryFunction),
          cell_mapper(grid.leafGridView())
      { analyzeView(); }

      const std::vector<NodeState> hangingNodes(const Cell& e) const
      {
        const typename GridView::IndexSet& indexSet = grid.leafGridView().indexSet();
        std::vector<NodeState> is_hanging;

        auto reference_element = referenceElement(e.geometry());

        // number of vertices in this element
        const IndexType v_size = reference_element.size(dim);

        // make sure the return array is big enough
        is_hanging.resize(v_size);

        // update minimum_level and maximum_level for vertices in this
        // cell
        // loop over vertices of the element
        for(IndexType i=0; i<v_size; ++i){
          const IndexType v_globalindex = indexSet.subIndex(e,i,dim);

          // here we make use of the fact that a node is hanging if and
          // only if it touches a cell of a level smaller than the
          // smallest level of all element containing the node
          const NodeInfo & v_info = node_info[v_globalindex];
          if(v_info.minimum_touching_level < v_info.minimum_level){
            is_hanging[i].is_hanging = true;
            is_hanging[i].is_boundary = v_info.is_boundary;
#ifndef NDEBUG
            if(verbosity>1){
              const Point & local  = reference_element.position(i,dim);
              const Point global = e.geometry().global(local);
              if(verbosity)
                std::cout << "Found hanging node with id " << v_globalindex << " at " << global << std::endl;
            }
#endif
          }
          else{
            is_hanging[i].is_hanging = false;
            is_hanging[i].is_boundary = v_info.is_boundary;
          }
        }

        return is_hanging;
      }

      void adaptToIsolatedHangingNodes()
      {
        if(verbosity)
          std::cout << "Begin isolation of hanging nodes" << std::endl;

        const typename GridView::IndexSet& indexSet = grid.leafGridView().indexSet();

        size_t iterations(0);

        bool reiterate;

        // Iterate until the isolation condition is achieved.
        do{
          size_t refinements(0);
          reiterate = false;

          GridView gv = grid.leafGridView();

          // loop over all codim<0> leaf elements of the partially refined grid
          for(const auto& cell : elements(gv)) {

            auto reference_element = referenceElement(cell.geometry());

            //std::cout << "cell center = " << it->geometry().center() << std::endl;

            // get the refinement level of the element
            const unsigned short level = cell.level();

            //std::cout << "level = " << level << std::endl;

            // number of vertices in this element
            const IndexType v_size = reference_element.size(dim);

            // update minimum_level and maximum_level for vertices in this
            // cell
            // loop over vertices of the element
            for(IndexType i=0; i<v_size; ++i){

              const IndexType v_globalindex = indexSet.subIndex(cell,i,dim);

              const NodeInfo & v_info = node_info[v_globalindex];

              //std::cout << "maximum_level = " << v_info.maximum_level << std::endl;

              const unsigned short level_diff = v_info.maximum_level - level;
              if( level_diff > 1){
                grid.mark(1, cell);   // Mark this element for an extra refinement if it has a hanging node belonging to a neighbouring element of a refinement level + 2 or more
                reiterate = true;    // Once an element has to be refined, the procedure needs to be repeated!
                refinements++;       // Count the number of refinements.

                if(verbosity>10){
                  // This will produce a lot of output on the screen!
                  std::cout << "   cell-id=" << cell_mapper.index(cell);
                  std::cout << "   level=" << level;
                  std::cout << "   v_size=" << v_size;
                  std::cout << "   v_globalindex = " << v_globalindex;
                  std::cout << std::endl;
                  std::cout << "Refining element nr " << cell_mapper.index(cell)
                            << " to isolate hanging nodes. Level diff = "
                            << v_info.maximum_level << " - " << level<< std::endl;
                }
                break;
              }
            } // end of loop over vertices


            if( cell.geometry().type().isSimplex() ){
              //
              // SPECIAL CASE for SIMPLICES:
              // Add extra check to find out "neighbouring" elements of level +2 or more
              // which share only a hanging node, but do not share an intersection
              // with the current element.
              //
              if( !reiterate ){

                //std::cout << "Extra check for SIMPLICES:" << std::endl;

                unsigned int intersection_index = 0;

                bool bJumpOut = false;

                // Loop over faces
                for(const auto& intersection : intersections(gv,cell)) {
                  ++intersection_index;

                  // only internal faces need to be taken care of
                  if(!intersection.boundary()) {

                    const int e_level = intersection.inside().level();
                    const int f_level = intersection.outside().level();

                    if( f_level > e_level ){

                      // We have to locate the potential hanging node
                      // for which we do the extra Check.

                      // get element-local index of the intersection
                      const int eLocalIndex =  intersection.indexInInside();

                      // Number of vertices on the face:
                      // A face(=edge) in a triangle has two vertices.
                      // A face(=triangle) in a tetrahedron has three vertices.
                      // const int e_v_size = reference_element.size(eLocalIndex,1,dim);

                      int nEdgeCenters = 0;
                      if( dim == 2 ){
                        // 2D-case: We need to check later for each vertex of the
                        // neigbouring element if it lies on the center of the element edge.
                        // Take care: fit->geometry().center() might return the face
                        // center of a refined neighbouring element!
                        // But we need the face center of the coarse face of the
                        // current element. Therefore loop over vertices on the face
                        // to calculate the proper face center for the coarse face!
                        nEdgeCenters = 1;
                      }
                      else{
                        // 3D-case: We need to check later for each vertex of the
                        // neigbouring element if it lies on the center of one of
                        // the 3 edges of the element face.
                        nEdgeCenters = 3;
                      }
                      std::vector<Dune::FieldVector<ctype,dim> >
                        edgecenter( nEdgeCenters, Dune::FieldVector<ctype,dim>(0) );
                      //std::cout << " edgecenter = " << edgecenter << std::endl;

                      // loop over center of the face (2d) or center of the edges of the face(3d)
                      for(int counter=0; counter<nEdgeCenters; ++counter){

                        int cornerIndex1 = counter % dim;
                        int cornerIndex2 = (counter+1) % dim;

                        const int e_v_index_1 =
                          reference_element.subEntity(eLocalIndex,1,cornerIndex1,dim);

                        const int e_v_index_2 =
                          reference_element.subEntity(eLocalIndex,1,cornerIndex2,dim);

                        auto vertex1 = cell.template subEntity<dim>(e_v_index_1);

                        auto vertex2 = cell.template subEntity<dim>(e_v_index_2);

                        edgecenter[counter] += vertex1.geometry().center();
                        edgecenter[counter] += vertex2.geometry().center();
                        edgecenter[counter] /= ctype( 2 );
                        //std::cout << " edgecenter = " << edgecenter << std::endl;


                        //
                        // check for the neighbouring element now...
                        //
                        auto nb_reference_element = referenceElement(intersection.outside().geometry());

                        // number of vertices in that neigbouring element
                        const IndexType nb_v_size = nb_reference_element.size(dim);

                        // loop over vertices of the neigbouring element
                        for(IndexType i=0; i<nb_v_size; ++i){

                          auto nb_vertex = intersection.outside().template subEntity<dim>(i);

                          bool doExtraCheck = false;

                          Dune::FieldVector<ctype,dim> center_diff ( edgecenter[counter] );

                          center_diff -= nb_vertex.geometry().center();

                          //std::cout << "nb_vertex = " << nb_vertex->geometry().center() << std::endl;

                          if( center_diff.two_norm2() < 5e-12 ){
                            doExtraCheck = true;
                          }


                          if( doExtraCheck ){

                            //std::cout << "doExtraCheck for node at "
                            // << nb_vertex->geometry().center() << std::endl;

                            const IndexType nb_v_globalindex = indexSet.index(nb_vertex);

                            const NodeInfo & nb_v_info = node_info[nb_v_globalindex];

                            const unsigned short level_diff = nb_v_info.maximum_level - level;

                            if( level_diff > 1){
                              bJumpOut = true;
                              grid.mark(1, cell);   // Mark this element for an extra refinement if it has a hanging node belonging to a neighbouring element of a refinement level + 2 or more
                              reiterate = true;    // Once an element has to be refined, the procedure needs to be repeated!
                              refinements++;       // Count the number of refinements.

                              if(verbosity>10){
                                // This will produce a lot of output on the screen!
                                std::cout << "   cell-id=" << cell_mapper.index(cell);
                                std::cout << "   level=" << level;
                                std::cout << "   v_size=" << v_size;
                                std::cout << std::endl;
                                std::cout << "   Extra refining for element nr " << cell_mapper.index(cell)
                                          << " to isolate hanging nodes. Level diff = "
                                          << nb_v_info.maximum_level << " - " << level<< std::endl;
                              }
                              break;

                            } // end if level_diff > 1

                          } // end if( doExtraCheck )
                          if( bJumpOut ) break;
                        } // end of loop over vertices of the neigbouring element
                        if( bJumpOut ) break;
                      } // end counter loop

                    } // end if( f_level > e_level )

                  } // end if not boundary
                  if( bJumpOut ) break;
                } // end of loop over faces of the element

              } // end if(!reiterate)

            } // end if geometry().type().isSimplex()

          } // end of loop over all codim<0> leaf elements


          if(reiterate){
            if(verbosity)
              std::cout << "Re-adapt for isolation of hanging nodes..." << std::endl;
            grid.preAdapt();
            grid.adapt();
            grid.postAdapt();
            analyzeView();
          }

          iterations++;
          if(verbosity)
            std::cout << "In iteration " << iterations << " there were "
                      << refinements << " grid cells to be refined additionally to isolate hanging nodes." << std::endl;
        }while(reiterate);
      }

    }; // end class HangingNodeManager

  } // end namespace PDELab
} // end namespace Dune
#endif // DUNE_PDELAB_CONSTRAINTS_HANGINGNODEMANAGER_HH
