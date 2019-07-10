// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>

#include <dune/common/filledarray.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/pdelab.hh>

// test function trees
template<int dim>
struct test;

template<>
struct test<2> {
  template<class GV>
  static void testleafgridfunction(const GV& gv)
  {
    // instantiate finite element maps
    auto gt = Dune::GeometryTypes::quadrilateral;
    typedef Dune::PDELab::P0LocalFiniteElementMap<float,double,GV::dimension> P0FEM;
    P0FEM p0fem(gt);
    typedef Dune::PDELab::QkLocalFiniteElementMap<GV,float,double,1> Q12DFEM;
    Q12DFEM q12dfem(gv);
    typedef Dune::PDELab::QkLocalFiniteElementMap<GV,float,double,2> Q22DFEM;
    Q22DFEM q22dfem(gv);

    // make a grid function space
    typedef Dune::PDELab::GridFunctionSpace<GV,P0FEM> P0GFS;
    P0GFS p0gfs(gv,p0fem);
    typedef Dune::PDELab::GridFunctionSpace<GV,Q12DFEM> GFS1;
    GFS1 gfs1(gv,q12dfem);
    typedef Dune::PDELab::GridFunctionSpace<GV,Q22DFEM> GFS2;
    GFS2 gfs2(gv,q22dfem);
    typedef Dune::PDELab::GridFunctionSpace<GV,Q22DFEM,Dune::PDELab::NoConstraints,
      Dune::PDELab::ISTL::VectorBackend<> > GFS3;
    GFS3 gfs3(gv,q22dfem);

    // test power
    typedef Dune::PDELab::PowerGridFunctionSpace<GFS2,2,Dune::PDELab::ISTL::VectorBackend<> > PGFS2;
    PGFS2 pgfs2(gfs2,gfs2);
    typedef Dune::PDELab::PowerGridFunctionSpace<GFS2,3,Dune::PDELab::ISTL::VectorBackend<> > PGFS3;
    PGFS3 pgfs3(gfs2,gfs2,gfs2);
    typedef Dune::PDELab::PowerGridFunctionSpace<GFS2,4,Dune::PDELab::ISTL::VectorBackend<> > PGFS4;
    PGFS4 pgfs4(gfs2,gfs2,gfs2,gfs2);
    typedef Dune::PDELab::PowerGridFunctionSpace<GFS2,5,Dune::PDELab::ISTL::VectorBackend<> > PGFS5;
    PGFS5 pgfs5(gfs2,gfs2,gfs2,gfs2,gfs2);
    typedef Dune::PDELab::PowerGridFunctionSpace<GFS2,6,Dune::PDELab::ISTL::VectorBackend<> > PGFS6;
    PGFS6 pgfs6(gfs2,gfs2,gfs2,gfs2,gfs2,gfs2);
    typedef Dune::PDELab::PowerGridFunctionSpace<GFS2,7,Dune::PDELab::ISTL::VectorBackend<> > PGFS7;
    PGFS7 pgfs7(gfs2,gfs2,gfs2,gfs2,gfs2,gfs2,gfs2);
    typedef Dune::PDELab::PowerGridFunctionSpace<GFS2,8,Dune::PDELab::ISTL::VectorBackend<> > PGFS8;
    PGFS8 pgfs8(gfs2,gfs2,gfs2,gfs2,gfs2,gfs2,gfs2,gfs2);
    typedef Dune::PDELab::PowerGridFunctionSpace<GFS2,9,Dune::PDELab::ISTL::VectorBackend<> > PGFS9;
    PGFS9 pgfs9(gfs2,gfs2,gfs2,gfs2,gfs2,gfs2,gfs2,gfs2,gfs2);
    typedef Dune::PDELab::PowerGridFunctionSpace<GFS2,10,Dune::PDELab::ISTL::VectorBackend<> > PGFS10;
    PGFS10 pgfs10(gfs2,gfs2,gfs2,gfs2,gfs2,gfs2,gfs2,gfs2,gfs2,gfs2);
    typedef Dune::PDELab::PowerGridFunctionSpace<GFS2,17,Dune::PDELab::ISTL::VectorBackend<> > PGFS17;
    PGFS17 pgfs17(gfs2);
    typedef Dune::PDELab::PowerGridFunctionSpace<GFS2,17,Dune::PDELab::ISTL::VectorBackend<>,
        Dune::PDELab:: EntityBlockedOrderingTag> PGFS17B;
    PGFS17B pgfs17b(gfs2);

    // make coefficent Vectors - we need to use copies of the spaces because the original
    // spaces are now part of a larger hierarchy
    {
      using V1 = Dune::PDELab::Backend::Vector<GFS1,double>;
      GFS1 gfs1(gv,q12dfem);
      V1 x1(gfs1);
      x1 = 0.0;
      using V2 = Dune::PDELab::Backend::Vector<GFS2,double>;
      GFS2 gfs2(gv,q22dfem);
      V2 x2(gfs2);
      x2 = 0.0;
    }

    // test composite
    typedef Dune::PDELab::CompositeGridFunctionSpace<Dune::PDELab::ISTL::VectorBackend<>,
        Dune::PDELab::LexicographicOrderingTag,GFS1,PGFS2> CGFS2;
    CGFS2 cgfs2(gfs1,pgfs2);
    typedef Dune::PDELab::CompositeGridFunctionSpace<Dune::PDELab::ISTL::VectorBackend<>,
        Dune::PDELab::LexicographicOrderingTag,GFS1,PGFS2,CGFS2> CGFS3;
    CGFS3 cgfs3(gfs1,pgfs2,cgfs2);
    typedef Dune::PDELab::CompositeGridFunctionSpace<Dune::PDELab::ISTL::VectorBackend<>,
        Dune::PDELab::LexicographicOrderingTag,GFS1,PGFS2,CGFS2,CGFS3> CGFS4;
    CGFS4 cgfs4(gfs1,pgfs2,cgfs2,cgfs3);
    typedef Dune::PDELab::CompositeGridFunctionSpace<Dune::PDELab::ISTL::VectorBackend<>,
        Dune::PDELab::LexicographicOrderingTag,GFS1,PGFS2,CGFS2,CGFS3,CGFS4> CGFS5;
    CGFS5 cgfs5(gfs1,pgfs2,cgfs2,cgfs3,cgfs4);
    typedef Dune::PDELab::CompositeGridFunctionSpace<Dune::PDELab::ISTL::VectorBackend<>,
        Dune::PDELab::LexicographicOrderingTag,GFS1,PGFS2,CGFS2,CGFS3,CGFS4,CGFS5> CGFS6;
    CGFS6 cgfs6(gfs1,pgfs2,cgfs2,cgfs3,cgfs4,cgfs5);
    typedef Dune::PDELab::CompositeGridFunctionSpace<Dune::PDELab::ISTL::VectorBackend<>,
        Dune::PDELab::LexicographicOrderingTag,GFS1,PGFS2,CGFS2,CGFS3,CGFS4,CGFS5,CGFS6> CGFS7;
    CGFS7 cgfs7(gfs1,pgfs2,cgfs2,cgfs3,cgfs4,cgfs5,cgfs6);
    typedef Dune::PDELab::CompositeGridFunctionSpace<Dune::PDELab::ISTL::VectorBackend<>,
        Dune::PDELab::LexicographicOrderingTag,GFS1,PGFS2,CGFS2,CGFS3,CGFS4,CGFS5,CGFS6,CGFS7> CGFS8;
    CGFS8 cgfs8(gfs1,pgfs2,cgfs2,cgfs3,cgfs4,cgfs5,cgfs6,cgfs7);
    typedef Dune::PDELab::CompositeGridFunctionSpace<Dune::PDELab::ISTL::VectorBackend<>,
        Dune::PDELab::LexicographicOrderingTag,GFS1,PGFS2,CGFS2,CGFS3,CGFS4,CGFS5,CGFS6,CGFS7,CGFS8> CGFS9;
    CGFS9 cgfs9(gfs1,pgfs2,cgfs2,cgfs3,cgfs4,cgfs5,cgfs6,cgfs7,cgfs8);
  }
};

template<>
struct test<3> {
  template<class GV>
  static void testleafgridfunction(const GV& gv)
  {
    // instantiate finite element maps
    auto gt = Dune::GeometryTypes::hexahedron;
    typedef Dune::PDELab::P0LocalFiniteElementMap<float,double,GV::dimension> P0FEM;
    P0FEM p0fem(gt);
    typedef Dune::PDELab::PkLocalFiniteElementMap<GV,float,double,1> P1FEM;
    P1FEM p1fem(gv);
    typedef Dune::PDELab::QkLocalFiniteElementMap<GV,float,double,1> Q1FEM;
    Q1FEM q1fem(gv);

    // make a grid function space
    typedef Dune::PDELab::GridFunctionSpace<GV,P0FEM> P0GFS;
    P0GFS p0gfs(gv,p0fem);
// Doesn't work, we need a grid with triangular elemets for that
//  typedef Dune::PDELab::GridFunctionSpace<GV,P1FEM> P1GFS;
//  P1GFS p1gfs(gv,p1fem);
    typedef Dune::PDELab::GridFunctionSpace<GV,Q1FEM> Q1GFS;
    Q1GFS q1gfs(gv,q1fem);
  }
};

template<class GV>
void testleafgridfunction(const GV& gv)
{
  test<GV::dimension>::testleafgridfunction(gv);
}


int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper::instance(argc, argv);

    // 2D
    {
      std::cout << "2D tests" << std::endl;
      // need a grid in order to test grid functions
      Dune::FieldVector<double,2> L(1.0);
      std::array<int,2> N(Dune::filledArray<2,int>(1));
      Dune::YaspGrid<2> grid(L,N);
      grid.globalRefine(1);

      testleafgridfunction(grid.leafGridView());
    }

    // 3D
    {
      std::cout << "3D tests" << std::endl;
      // need a grid in order to test grid functions
      Dune::FieldVector<double,3> L(1.0);
      std::array<int,3> N(Dune::filledArray<3,int>(1));
      Dune::YaspGrid<3> grid(L,N);
      grid.globalRefine(1);

      testleafgridfunction(grid.leafGridView());
    }

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
