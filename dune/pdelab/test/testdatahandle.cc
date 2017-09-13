// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/pdelab/finiteelementmap/p0fem.hh>
#include <dune/pdelab/finiteelementmap/pkfem.hh>
#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/pdelab/finiteelementmap/qkdg.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/powergridfunctionspace.hh>
#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include <dune/pdelab/gridfunctionspace/vtk.hh>

template<int codim, typename GFS>
void print_info(const GFS& gfs)
{
  std::cout << "codim " << codim << ":" << std::endl;
  std::cout << "contains = " << std::boolalpha << gfs.dataHandleContains(codim) << "  "
            << "fixed size = " << gfs.dataHandleFixedSize(codim) << std::endl;
  std::vector<typename GFS::Ordering::Traits::ContainerIndex> indices;
  typedef typename GFS::Traits::GridView GV;
  GV gv = gfs.gridView();
  for (typename GV::template Codim<codim>::Iterator it = gv.template begin<codim>();
       it != gv.template end<codim>();
       ++it)
    {
      std::size_t size = gfs.dataHandleSize(*it);
      std::cout << "entity index = " << gv.indexSet().index(*it) << "  "
                << "size = " << size << std::endl;
      indices.resize(size);
      gfs.dataHandleContainerIndices(*it,indices);
      for (typename std::vector<typename GFS::Ordering::Traits::ContainerIndex>::const_iterator it = indices.begin();
           it != indices.end();
           ++it)
        std::cout << *it << std::endl;
    }
}

template<typename GFS>
void info(const GFS& gfs, std::string name)
{
  std::cout << name << std::endl
            << gfs.ordering().size() << std::endl;
  print_info<0>(gfs);
  //print_info<1>(gfs);
  print_info<2>(gfs);

  std::cout << std::endl << std::endl;
}

template<class GV>
static void testdatahandle(const GV& gv)
{
  // instantiate finite element maps
  auto gt = Dune::GeometryTypes::quadrilateral;
  typedef Dune::PDELab::P0LocalFiniteElementMap<float,double,GV::dimension> P0FEM;
  P0FEM p0fem(gt);
  typedef Dune::PDELab::QkLocalFiniteElementMap<GV,float,double,1> Q12DFEM;
  Q12DFEM q12dfem(gv);
  typedef Dune::PDELab::QkLocalFiniteElementMap<GV,float,double,2> Q22DFEM;
  Q22DFEM q22dfem(gv);
  typedef Dune::PDELab::QkDGLocalFiniteElementMap<double,double,2,2> DG22DFEM;
  DG22DFEM dg22dfem;

  typedef Dune::PDELab::NoConstraints NoConstraints;
  typedef Dune::PDELab::ISTL::VectorBackend<> VBE;

  // make a grid function space
  typedef Dune::PDELab::GridFunctionSpace<GV,P0FEM,NoConstraints,VBE> P0GFS;
  P0GFS p0gfs(gv,p0fem);
  p0gfs.name("p0gfs");
  typedef Dune::PDELab::GridFunctionSpace<GV,Q12DFEM,NoConstraints,VBE> GFS1;
  GFS1 gfs1(gv,q12dfem);
  gfs1.name("gfs1");
  typedef Dune::PDELab::GridFunctionSpace<GV,DG22DFEM,NoConstraints,VBE> GFS2;
  GFS2 gfs2(gv,dg22dfem);
  gfs2.name("gfs2");
  typedef Dune::PDELab::GridFunctionSpace<GV,Q22DFEM,NoConstraints,VBE> GFS3;
  GFS3 gfs3(gv,q22dfem);
  gfs3.name("gfs3");

  using V = Dune::PDELab::Backend::Vector<GFS1,int>;

  V v(gfs1);
  v = 1;

  Dune::PDELab::AddDataHandle<GFS1,V> adddh(gfs1,v);
  gfs1.gridView().communicate(adddh,Dune::All_All_Interface,Dune::ForwardCommunication);

  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
  Dune::PDELab::addSolutionToVTKWriter(vtkwriter,gfs1,v);
  vtkwriter.write("testdatahandle-gfs1",Dune::VTK::ascii);

#if 0
  info(p0gfs,"P0");
  info(gfs1,"P1");
  info(gfs2,"DG2");
  info(gfs3,"Q2");

  typedef Dune::PDELab::ISTL::VectorBackend<> VBE;

  typedef Dune::PDELab::PowerGridFunctionSpace<GFS1,3,VBE> PGFS1;
  PGFS1 pgfs1(gfs1);

  info(pgfs1,"P1^3");

  typedef Dune::PDELab::PowerGridFunctionSpace<GFS1,3,VBE,Dune::PDELab::EntityBlockedOrderingTag> PGFS2;
  PGFS2 pgfs2(gfs1);

  info(pgfs2,"P1^3 (entity-wise blocked)");

  typedef Dune::PDELab::CompositeGridFunctionSpace<
    VBE,
    Dune::PDELab::LexicographicOrderingTag,
    GFS1,
    GFS2
    > CGFS1;

  CGFS1 cgfs1(gfs1,gfs2);

  info(cgfs1,"P1 * DG2");

  typedef Dune::PDELab::CompositeGridFunctionSpace<
    VBE,
    Dune::PDELab::EntityBlockedOrderingTag,
    GFS1,
    GFS2,
    GFS1
    > CGFS2;

  CGFS2 cgfs2(gfs1,gfs2,gfs1);

  info(cgfs2,"P1 * DG2 * P1 (entity-wise blocked)");

  typedef Dune::PDELab::PowerGridFunctionSpace<
    GFS3,
    3,
    Dune::PDELab::ISTL::VectorBackend<
      Dune::PDELab::ISTL::Blocking::fixed
      >,
    Dune::PDELab::EntityBlockedOrderingTag
    > VGFS;
  VGFS vgfs(gfs3);

  info(vgfs,"P2^3 (entity-wise blocked, matrix blocks)");

  typedef Dune::PDELab::CompositeGridFunctionSpace<
    Dune::PDELab::ISTL::VectorBackend<
      Dune::PDELab::ISTL::Blocking::bcrs
      >,
    Dune::PDELab::LexicographicOrderingTag,
    GFS1,
    VGFS
    > THGFS;

  THGFS thgfs(gfs1,vgfs);

  info(thgfs,"P1 * (P2^3) (entity-wise blocked with matrix blocks at P2 level, macro blocks");

#if 0
  typedef Dune::PDELab::GridFunctionSpace<GV,Q22DFEM,Dune::PDELab::NoConstraints,
                                          Dune::PDELab::ISTL::VectorBackend<1>,
                                          Dune::PDELab::GridFunctionRestrictedMapper> GFS3;
  GFS3 gfs3(gv,q22dfem);

  // test power
  typedef Dune::PDELab::PowerGridFunctionSpace<GFS2,2> PGFS2;
  PGFS2 pgfs2(gfs2,gfs2);
  typedef Dune::PDELab::PowerGridFunctionSpace<GFS2,3> PGFS3;
  PGFS3 pgfs3(gfs2,gfs2,gfs2);
  typedef Dune::PDELab::PowerGridFunctionSpace<GFS2,4> PGFS4;
  PGFS4 pgfs4(gfs2,gfs2,gfs2,gfs2);
  typedef Dune::PDELab::PowerGridFunctionSpace<GFS2,5> PGFS5;
  PGFS5 pgfs5(gfs2,gfs2,gfs2,gfs2,gfs2);
  typedef Dune::PDELab::PowerGridFunctionSpace<GFS2,6> PGFS6;
  PGFS6 pgfs6(gfs2,gfs2,gfs2,gfs2,gfs2,gfs2);
  typedef Dune::PDELab::PowerGridFunctionSpace<GFS2,7> PGFS7;
  PGFS7 pgfs7(gfs2,gfs2,gfs2,gfs2,gfs2,gfs2,gfs2);
  typedef Dune::PDELab::PowerGridFunctionSpace<GFS2,8> PGFS8;
  PGFS8 pgfs8(gfs2,gfs2,gfs2,gfs2,gfs2,gfs2,gfs2,gfs2);
  typedef Dune::PDELab::PowerGridFunctionSpace<GFS2,9> PGFS9;
  PGFS9 pgfs9(gfs2,gfs2,gfs2,gfs2,gfs2,gfs2,gfs2,gfs2,gfs2);
  typedef Dune::PDELab::PowerGridFunctionSpace<GFS2,10> PGFS10;
  PGFS10 pgfs10(gfs2,gfs2,gfs2,gfs2,gfs2,gfs2,gfs2,gfs2,gfs2,gfs2);
  typedef Dune::PDELab::PowerGridFunctionSpace<GFS2,17> PGFS17;
  PGFS17 pgfs17(gfs2);
  typedef Dune::PDELab::PowerGridFunctionSpace<GFS2,17,Dune::PDELab::GridFunctionSpaceBlockwiseMapper> PGFS17B;
  PGFS17B pgfs17b(gfs2);

  // make coefficent Vectors
  using V1 = Dune::PDELab::Backend::Vector<GFS1,double>;
  V1 x1(gfs1);
  x1 = 0.0;
  using V2 = Dune::PDELab::Backend::Vector<GFS2,double>;
  V2 x2(gfs2);
  x2 = 0.0;

  // test composite
  typedef Dune::PDELab::CompositeGridFunctionSpace<Dune::PDELab::GridFunctionSpaceLexicographicMapper,
                                                   GFS1,PGFS2> CGFS2;
  CGFS2 cgfs2(gfs1,pgfs2);
  typedef Dune::PDELab::CompositeGridFunctionSpace<Dune::PDELab::GridFunctionSpaceLexicographicMapper,
                                                   GFS1,PGFS2,CGFS2> CGFS3;
  CGFS3 cgfs3(gfs1,pgfs2,cgfs2);
  typedef Dune::PDELab::CompositeGridFunctionSpace<Dune::PDELab::GridFunctionSpaceLexicographicMapper,
                                                   GFS1,PGFS2,CGFS2,CGFS3> CGFS4;
  CGFS4 cgfs4(gfs1,pgfs2,cgfs2,cgfs3);
  typedef Dune::PDELab::CompositeGridFunctionSpace<Dune::PDELab::GridFunctionSpaceLexicographicMapper,
                                                   GFS1,PGFS2,CGFS2,CGFS3,CGFS4> CGFS5;
  CGFS5 cgfs5(gfs1,pgfs2,cgfs2,cgfs3,cgfs4);
  typedef Dune::PDELab::CompositeGridFunctionSpace<Dune::PDELab::GridFunctionSpaceLexicographicMapper,
                                                   GFS1,PGFS2,CGFS2,CGFS3,CGFS4,CGFS5> CGFS6;
  CGFS6 cgfs6(gfs1,pgfs2,cgfs2,cgfs3,cgfs4,cgfs5);
  typedef Dune::PDELab::CompositeGridFunctionSpace<Dune::PDELab::GridFunctionSpaceLexicographicMapper,
                                                   GFS1,PGFS2,CGFS2,CGFS3,CGFS4,CGFS5,CGFS6> CGFS7;
  CGFS7 cgfs7(gfs1,pgfs2,cgfs2,cgfs3,cgfs4,cgfs5,cgfs6);
  typedef Dune::PDELab::CompositeGridFunctionSpace<Dune::PDELab::GridFunctionSpaceLexicographicMapper,
                                                   GFS1,PGFS2,CGFS2,CGFS3,CGFS4,CGFS5,CGFS6,CGFS7> CGFS8;
  CGFS8 cgfs8(gfs1,pgfs2,cgfs2,cgfs3,cgfs4,cgfs5,cgfs6,cgfs7);
  typedef Dune::PDELab::CompositeGridFunctionSpace<Dune::PDELab::GridFunctionSpaceLexicographicMapper,
                                                   GFS1,PGFS2,CGFS2,CGFS3,CGFS4,CGFS5,CGFS6,CGFS7,CGFS8> CGFS9;
  CGFS9 cgfs9(gfs1,pgfs2,cgfs2,cgfs3,cgfs4,cgfs5,cgfs6,cgfs7,cgfs8);

#endif
#endif
}


int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi

    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

    // 2D
    {
      if (helper.rank() == 0)
        std::cout << "2D tests" << std::endl;
      // need a grid in order to test grid functions
      Dune::FieldVector<double,2> L(1.0);
      L[0] *= 5;
      std::array<int,2> N(Dune::fill_array<int,2>(2));
      N[0] = 10;
      std::bitset<2> B(false);
      Dune::YaspGrid<2> grid(helper.getCommunicator(),L,N,B,1);
      // grid.globalRefine(1);

      testdatahandle(grid.leafGridView());
    }

#if 0
    // 3D
    {
      std::cout << "3D tests" << std::endl;
      // need a grid in order to test grid functions
      Dune::FieldVector<double,3> L(1.0);
      std::array<int,3> N(Dune::fill_array<int,3>(1));
      std::bitset<3> B(false);
      Dune::YaspGrid<3> grid(L,N,B,1);
      grid.globalRefine(1);

      testleafgridfunction(grid.leafGridView());
    }
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
