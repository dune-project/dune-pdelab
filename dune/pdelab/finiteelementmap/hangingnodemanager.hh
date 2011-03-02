// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef HANGINGNODEMANAGER_HH
#define HANGINGNODEMANAGER_HH

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
      enum{ verbosity = 1 };
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
      typedef typename GridView::template Codim<0>::EntityPointer CellEntityPointer;
      typedef typename GridView::template Codim<dim>::EntityPointer VertexEntityPointer;
      typedef typename GridView::template Codim<0>::Iterator Iterator;
      typedef typename GridView::IntersectionIterator IntersectionIterator;
      typedef typename GridView::Grid::ctype ctype;
      typedef typename Dune::FieldVector<ctype,dim> Point;
      typedef typename Dune::FieldVector<ctype,dim-1> FacePoint;
  
      typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView,
        MCMGElementLayout> CellMapper;
      typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView,
        MCMGVertexLayout> VertexMapper;

      Grid & grid;
      const BoundaryFunction & boundaryFunction;
      CellMapper cell_mapper;
      VertexMapper vertex_mapper;

    public:
      size_t map(VertexEntityPointer & e){return vertex_mapper.map(*e);}

      void analyzeView()
      {
	cell_mapper.update();
	vertex_mapper.update();

	node_info = std::vector<NodeInfo>(vertex_mapper.size());

	const GridView & gv = grid.leafView();

	Iterator it = gv.template begin<0>();
	Iterator eit = gv.template end<0>();

	for(;it!=eit;++it){

	  const Dune::GenericReferenceElement<double,dim> & 
	    reference_element = 
	    Dune::GenericReferenceElements<double,dim>::general(it->geometry().type()); 

      
	  // level of this element
	  const unsigned short level = it->level();
      
	  // number of vertices in this element
	  const IndexType v_size = reference_element.size(dim);

	  // update minimum_level and maximum_level for vertices in this
	  // cell
	  for(IndexType i=0; i<v_size; ++i){
	    const VertexEntityPointer vertex = it->template subEntity<dim>(i);
	    const IndexType v_globalindex = vertex_mapper.map( *vertex );
	    unsigned short & min = node_info[v_globalindex].minimum_level;
	    unsigned short & max = node_info[v_globalindex].maximum_level;
	    if (level < min) min = level;
	    if (level > max) max = level;
	  }

	  // Now we still have to update minimum_touching_level for this
	  // cell

          unsigned int intersection_index = 0;
	  IntersectionIterator fit = gv.ibegin(*it);
	  IntersectionIterator efit = gv.iend(*it);
	  typedef typename IntersectionIterator::Intersection Intersection;

	  for(;fit!=efit;++fit,++intersection_index){

	    const Dune::GenericReferenceElement<double,dim-1> & 
	      reference_face_element = 
	      Dune::GenericReferenceElements<double,dim-1>::general(fit->geometry().type()); 

	    const int eLocalIndex =  fit->indexInInside();
	    const int e_level = fit->inside()->level();

	    // numbero of vertices in face
	    const int e_v_size = reference_element.size(eLocalIndex,1,dim);

	    if((*fit).boundary()) {
	      for(int i=0; i<e_v_size;++i){
		const int e_v_index = reference_element.subEntity(eLocalIndex,1,i,dim);
		const VertexEntityPointer vertex = it->template subEntity<dim>(e_v_index);
		const IndexType v_globalindex = vertex_mapper.map( *vertex );

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

                if( boundaryFunction.isDirichlet( IntersectionGeometry<Intersection>(*fit,intersection_index),
                                                  facelocal_position) )
		  node_info[v_globalindex].is_boundary=true;
                else
		  node_info[v_globalindex].is_boundary=false;

		unsigned short & min = node_info[v_globalindex].minimum_touching_level;
		if( e_level < min) min = e_level;
	      }
	      continue;
	    }

	    const int f_level = fit->outside()->level();
	  
	    // a conforming face has no hanging nodes
	    if(fit->conforming())
	      continue;

	    // so far no support for initially non conforming grids
	    assert(e_level != f_level);

	    // this check needs to be performed on the element containing
	    // the hanging node only
	    if(e_level < f_level)
	      continue;

	    for(int i=0; i<e_v_size;++i){
	      const int e_v_index = reference_element.subEntity(eLocalIndex,1,i,dim);
	      const VertexEntityPointer vertex = it->template subEntity<dim>(e_v_index);
	      const IndexType v_globalindex = vertex_mapper.map( *vertex );
	      unsigned short & min = node_info[v_globalindex].minimum_touching_level;
	      if( f_level < min) min = f_level;
	    }

	  }
	}
      }

      HangingNodeManager(Grid & _grid, const BoundaryFunction & _boundaryFunction)
	: grid(_grid),
	  boundaryFunction(_boundaryFunction),
	  cell_mapper(grid.leafView()),
	  vertex_mapper(grid.leafView())
      { analyzeView(); }

      const std::vector<NodeState> hangingNodes(const CellEntityPointer & e) const
      {
	std::vector<NodeState> is_hanging;
    
	const Dune::GenericReferenceElement<double,dim> & 
	  reference_element = 
	  Dune::GenericReferenceElements<double,dim>::general(e->geometry().type()); 

	// number of vertices in this element
	const IndexType v_size = reference_element.size(dim);

	// make sure the return array is big enough
	is_hanging.resize(v_size);

	// update minimum_level and maximum_level for vertices in this
	// cell
	for(IndexType i=0; i<v_size; ++i){
	  const VertexEntityPointer & vertex = e->template subEntity<dim>(i);
	  const IndexType v_globalindex = vertex_mapper.map( *vertex );

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
	      const Point global = e->geometry().global(local);
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

	size_t iterations(0);

	bool reiterate;
    
	// Iterate until the isolation condition is achieved. 
	do{
	  size_t refinements(0);
	  reiterate = false;
      
	  const GridView & gv = grid.leafView();

	  Iterator it = gv.template begin<0>();
	  Iterator eit = gv.template end<0>();

	  for(;it!=eit;++it){

	    const Dune::GenericReferenceElement<double,dim> & 
	      reference_element = 
	      Dune::GenericReferenceElements<double,dim>::general(it->geometry().type()); 

	    const unsigned short level = it->level();
	    // number of vertices in this element
	    const IndexType v_size = reference_element.size(dim);

	    // update minimum_level and maximum_level for vertices in this
	    // cell
	    for(IndexType i=0; i<v_size; ++i){
	      const VertexEntityPointer & vertex = it->template subEntity<dim>(i);
	      const IndexType v_globalindex = vertex_mapper.map( *vertex );
	      const NodeInfo & v_info = node_info[v_globalindex];

	      const unsigned short level_diff = v_info.maximum_level - level;
	      if( level_diff > 1){
		grid.mark(1, *it);
		reiterate = true;
		refinements++;
		if(verbosity>1){
		  std::cout << "   Refining element nr " << cell_mapper.map(*it) 
			    << " to isolate hanging nodes " 
			    << v_info.maximum_level << " - " << level<< std::endl;
		}
		break;
	      }
	    } // i

	  } // it

	  if(reiterate){
	    grid.preAdapt();
	    grid.adapt();
	    grid.postAdapt();
	    analyzeView();
	  }

	  iterations++;
	  if(verbosity)
	    std::cout << "In iteration " << iterations << " refined " 
		      << refinements << " grid cells" << std::endl;
	}while(reiterate);
      }

    };

  }
}
#endif 
