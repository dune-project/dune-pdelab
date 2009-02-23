// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifdef HAVE_CONFIG_H
#include "config.h"     
#endif
#include<iostream>
#include<vector>
#include<map>
#include<dune/common/mpihelper.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/static_assert.hh>
#include<dune/grid/yaspgrid.hh>
#include"../finiteelementmap/p0fem.hh"
#include"../finiteelementmap/p12dfem.hh"
#include"../finiteelementmap/pk2dfem.hh"
#include"../gridfunctionspace/gridfunctionspace.hh"
#include"../gridfunctionspace/gridfunctionspaceutilities.hh"
#include"../gridfunctionspace/interpolate.hh"
#include"../common/function.hh"
#include"../common/vtkexport.hh"
#include"gridexamples.hh"


// define some grid functions to interpolate from
template<typename GV, typename RF>
class F
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
                                                  F<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,F<GV,RF> > BaseT;

  F (const GV& gv) : BaseT(gv) {}
  inline void evaluateGlobal (const typename Traits::DomainType& x, 
							  typename Traits::RangeType& y) const
  {
    typename Traits::DomainType center;
    for (int i=0; i<GV::dimension; i++) center[i] = 0.5;
    center -= x;
	y = exp(-center.two_norm2());
  }
};

// define some grid functions to interpolate from
template<typename GV>
class B
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,int,1>,
                                                  B<GV> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,int,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,B<GV> > BaseT;

  B (const GV& gv) : BaseT(gv) {}
  inline void evaluateGlobal (const typename Traits::DomainType& x, 
							  typename Traits::RangeType& y) const
  {
    if (x[0]<1E-6)
      y = 1; // Dirichlet boundary
    else
      y = 0;
  }
};


// Parameter class for grid function space that assembles constraints
class P12DConstraintsAssembler
{
public:
  enum { assembleBoundary = true };
  enum { assembleIntersection = false };

  // boundary constraints
  // F : grid function returning boundary condition type
  // IG : intersection geometry
  // LFS : local function space
  // T : TransformationType
  template<typename F, typename IG, typename LFS, typename T>
  static void boundary (const F& f, const IG& ig, const LFS& lfs, T& trafo)
  {
    // 2D here, get midpoint of edge
    typedef typename IG::LocalGeometry::ctype ctype;
    typedef Dune::FieldVector<ctype,IG::LocalGeometry::mydimension> ILPoint;
    typedef Dune::FieldVector<ctype,IG::LocalGeometry::coorddimension> ELPoint;
    ILPoint ip(0.5);
    ELPoint ep = ig.intersectionSelfLocal().global(ip);

    // determine type of boundary condition
	typename F::Traits::RangeType bctype;
    f.evaluate(ep,bctype);

    // if dirichlet boundary, the two end nodes of the edge are constrained
    if (bctype>0)
      {
        // determine the constrained nodes (requires knowledge about reference element)
        typename T::value_type::second_type empty;
        int edge = ig.numberInSelf();
        trafo[(edge+1)%3] = empty; // first node
        trafo[(edge+2)%3] = empty; // second node
      }
  }
};

template<typename I>
class IntersectionGeometry
{
public:
  typedef typename I::Geometry Geometry;
  typedef typename I::LocalGeometry LocalGeometry;
  typedef typename I::Entity Entity;
  typedef typename I::EntityPointer EntityPointer;
  typedef typename Entity::ctype ctype;
  enum { dimension=Entity::dimension };
  enum { dimensionworld=Entity::dimensionworld };

  IntersectionGeometry (const I& i_)
    : i(i_)
  {}


  //! return true if intersection is with interior or exterior boundary (see the cases above)
    bool boundary () const
    {
      return i.boundary();
    }

  /**
     \brief Identifier for boundary segment from macro grid.

     One can attach a boundary Id to a boundary segment on the macro
     grid. This Id will also be used for all fragments of these
     boundary segments.

     The numbering is defined as:
     - Id==0 for all intersections without boundary()==false
     - Id>=0 for all intersections without boundary()==true
     
     The way the Identifiers are attached to the grid may differ
     between the different grid implementations.

   */
  int boundaryId () const
  {
    return i.boundaryId();
  } 

  //! @brief return true if intersection is shared with another element.
  bool neighbor () const 
    {
      return i.neighbor();
    }

  /*! @brief geometrical information about this intersection in local
    coordinates of the inside() entity.

    This method returns a Geometry object that provides a mapping from
    local coordinates of the intersection to local coordinates of the
    inside() entity.
  */
  const LocalGeometry& intersectionSelfLocal () const
    {
      return i.intersectionSelfLocal();
    }
  /*! @brief geometrical information about this intersection in local
    coordinates of the outside() entity.

    This method returns a Geometry object that provides a mapping from
    local coordinates of the intersection to local coordinates of the
    outside() entity.
  */
  const LocalGeometry& intersectionNeighborLocal () const
    {
      return i.intersectionNeighborLocal();
    }

  /*! @brief geometrical information about this intersection in global coordinates.

    This method returns a Geometry object that provides a mapping from
    local coordinates of the intersection to global (world) coordinates.
  */
  const Geometry& intersectionGlobal () const
    {
      return i.intersectionGlobal();
    }

  //! Local number of codim 1 entity in the inside() Entity where intersection is contained in
  int numberInSelf () const
    {
      return i.numberInSelf ();
    }

  //! Local number of codim 1 entity in outside() Entity where intersection is contained in
  int numberInNeighbor () const
    {
      return i.numberInNeighbor ();
    }

  /*! @brief Return an outer normal (length not necessarily 1)

    The returned vector may depend on local position within the intersection.
  */
  Dune::FieldVector<ctype, dimensionworld> outerNormal (const Dune::FieldVector<ctype, dimension-1>& local) const
    {
      return i.outerNormal(local);
    }

  /*! @brief return outer normal scaled with the integration element
	@copydoc outerNormal
    The normal is scaled with the integration element of the intersection. This
	method is redundant but it may be more efficent to use this function
	rather than computing the integration element via intersectionGlobal().
  */
  Dune::FieldVector<ctype, dimensionworld> integrationOuterNormal (const Dune::FieldVector<ctype, dimension-1>& local) const
    {
      return i.integrationOuterNormal(local);
    }

  /*! @brief Return unit outer normal (length == 1)

  The returned vector may depend on the local position within the intersection.
  It is scaled to have unit length.
  */
  Dune::FieldVector<ctype, dimensionworld> unitOuterNormal (const Dune::FieldVector<ctype, dimension-1>& local) const
    {
      return i.unitOuterNormal(local);
    }

  /*! @brief return EntityPointer to the Entity on the inside of this
    intersection. That is the Entity where we started this .
   */
  EntityPointer inside() const
    {
      return i.inside();
    }

  /*! @brief return EntityPointer to the Entity on the outside of this
    intersection. That is the neighboring Entity.

    @warning Don't call this method if there is no neighboring Entity
    (neighbor() returns false). In this case the result is undefined.
   */
  EntityPointer outside() const
    {
      return i.outside();
    }

private:
  const I& i;
};



// generate a P1 function and output it
template<class GV> 
void testpk (const GV& gv)
{
  typedef typename GV::Grid::ctype DF;
  const int dim = GV::dimension;

  // instantiate finite element maps
  typedef Dune::PDELab::P12DLocalFiniteElementMap<DF,double> P1FEM;
  P1FEM p1fem;
  
  // make a grid function space
  typedef Dune::PDELab::GridFunctionSpace<GV,P1FEM,P12DConstraintsAssembler> P1GFSU; 
  P1GFSU p1gfsu(gv,p1fem);

  // make constrained space
  typedef Dune::PDELab::ProxyGridFunctionSpace<P1GFSU> P1GFS; 
  P1GFS p1gfs(p1gfsu);

  // make coefficent Vectors
  typedef typename P1GFS::template VectorContainer<double>::Type P1V;
  P1V p1xg(p1gfs);
  p1xg = 0.0;

  // construct a grid function
  typedef F<GV,double> FType;
  FType f(gv);
  typedef B<GV> BType;
  BType b(gv);

  // do interpolation
  Dune::PDELab::interpolate(f,p1gfs,p1xg);

  // make local function space object
  typename P1GFSU::LocalFunctionSpace p1lfs(p1gfsu);

  // loop over elements
  typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
  typedef typename GV::IntersectionIterator IntersectionIterator;
  typedef typename IntersectionIterator::Intersection Intersection;
  for (ElementIterator it = gv.template begin<0>();
	   it!=gv.template end<0>(); ++it)
	{
      p1lfs.bind(*it);
      typename P1GFSU::Traits::TransformationType trafo;
      IntersectionIterator endit = gv.iend(*it);
      for (IntersectionIterator iit = gv.ibegin(*it); iit!=endit; ++iit)
        {
          if (iit->boundary())
            {
              P1GFSU::Traits::ConstraintsAssemblerType::boundary(Dune::PDELab::GridFunctionToLocalFunctionAdapter<BType>(b,*it),
                                                                 IntersectionGeometry<Intersection>(*iit),
                                                                 p1lfs,trafo);
            }
        }

      typedef typename P1GFSU::Traits::TransformationType::iterator mapiterator;

      if (trafo.size()>0)
        {
          for (mapiterator i=trafo.begin(); i!=trafo.end(); ++i)
            std::cout << p1lfs.globalIndex(i->first) << " ";
          std::cout << std::endl;
        }
    }

  // make discrete function object
  typedef Dune::PDELab::DiscreteGridFunction<P1GFS,P1V> P1DGF;
  P1DGF p1dgf(p1gfs,p1xg);

  // output grid function with VTKWriter
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<P1DGF>(p1dgf,"p1"));
  vtkwriter.write("testconstraints",Dune::VTKOptions::ascii);
}

int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

#if HAVE_ALBERTA
 	AlbertaLDomain albertagrid;
  	albertagrid.globalRefine(4);
    testpk(albertagrid.leafView());
#endif

	// test passed
	return 0;

  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
	return 1;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
	return 1;
  }
} 
