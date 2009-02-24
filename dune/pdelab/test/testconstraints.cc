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
#include"../gridfunctionspace/constraints.hh"
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

// define some boundary grid functions to define boundary conditions
template<typename GV>
class B
  : public Dune::PDELab::BoundaryGridFunctionBase<Dune::PDELab::BoundaryGridFunctionTraits<GV,int,1,
                                                                                           Dune::FieldVector<int,1> >,
                                                  B<GV> >
{
  const GV& gv;

public:
  typedef Dune::PDELab::BoundaryGridFunctionTraits<GV,int,1,Dune::FieldVector<int,1> > Traits;
  typedef Dune::PDELab::BoundaryGridFunctionBase<Traits,B<GV> > BaseT;

  B (const GV& gv_) : gv(gv_) {}

  template<typename I>
  inline void evaluate (const Dune::PDELab::IntersectionGeometry<I>& ig, 
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {  
    typedef Dune::PDELab::IntersectionGeometry<I> IG;

    // map from local coordinates in intersection to global coordinates
    Dune::FieldVector<typename IG::ctype,IG::Geometry::coorddimension> 
      xg = ig.intersectionGlobal().global(x);

    // set boundary condition according to coordinates
    // here we could also use boundaryid etc.
    if (xg[0]<1E-6)
      y = 1; // Dirichlet boundary
    else
      y = 0;
  }

  //! get a reference to the GridView
  inline const GV& getGridView ()
  {
    return gv;
  }
};


// Parameter class for grid function space that assembles constraints
class P12DConstraints
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
  void boundary (const F& f, const IG& ig, const LFS& lfs, T& trafo) const
  {
    // 2D here, get midpoint of edge
    typename F::Traits::DomainType ip(0.5);

    // determine type of boundary condition
	typename F::Traits::RangeType bctype;
    f.evaluate(ig,ip,bctype);

    // if dirichlet boundary, the two end nodes of the edge are constrained
    if (bctype>0)
      {
        // determine the constrained nodes (requires knowledge about reference element)
        typename T::RowType empty;
        int edge = ig.numberInSelf();
        trafo[(edge+1)%3] = empty; // first node
        trafo[(edge+2)%3] = empty; // second node
      }
  }
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
  
  // make constrained space
  typedef Dune::PDELab::GridFunctionSpace<GV,P1FEM,P12DConstraints> P1GFS; 
  P1GFS p1gfs(gv,p1fem);

  // make coefficent Vectors
  typedef typename P1GFS::template VectorContainer<double>::Type P1V;
  P1V p1xg(p1gfs);
  p1xg = 0.0;

  // make constraints map
  typedef typename P1GFS::template ConstraintsContainer<double>::Type P1C;
  P1C p1cg;
  p1cg.clear();

  // interpolate from grid function
  typedef F<GV,double> FType;
  FType f(gv);
  Dune::PDELab::interpolate(f,p1gfs,p1xg);

  // set up constraints from boundary condition function
  typedef B<GV> BType;
  BType b(gv);
  Dune::PDELab::constraints(b,p1gfs,p1cg);

  // set Dirichlet nodes to zero
  Dune::PDELab::set(p1cg,0.0,p1xg);

  // make discrete function object
  typedef Dune::PDELab::DiscreteGridFunction<P1GFS,P1V> P1DGF;
  P1DGF p1dgf(p1gfs,p1xg);

  // output grid function with VTKWriter
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<P1DGF>(p1dgf,"p1"));
  vtkwriter.write("testconstraintsp1",Dune::VTKOptions::ascii);
}

int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

#if HAVE_ALUGRID
//  	ALUUnitSquare alugrid;
//   	alugrid.globalRefine(3);
//     testpk(alugrid.leafView());
#endif

#if HAVE_UG
 	UGUnitSquare uggrid;
  	uggrid.globalRefine(3);
    testpk(uggrid.leafView());
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
