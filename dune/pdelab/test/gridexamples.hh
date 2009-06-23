// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_GRIDEXAMPLES_HH
#define DUNE_PDELAB_GRIDEXAMPLES_HH

#include<dune/grid/yaspgrid.hh>
#if HAVE_ALBERTA
#include<dune/grid/albertagrid.hh>
#include <dune/grid/albertagrid/dgfparser.hh>
#include <dune/grid/albertagrid/gridfactory.hh>
#endif
#if HAVE_UG 
#include <dune/grid/uggrid/uggridfactory.hh>
#endif
#if HAVE_ALUGRID
#include<dune/grid/alugrid.hh>
#endif

class YaspUnitSquare : public Dune::YaspGrid<2>
{
public:
  YaspUnitSquare () : Dune::YaspGrid<2>(Dune::FieldVector<double,2>(1.0),
					  Dune::FieldVector<int,2>(1),
					  Dune::FieldVector<bool,2>(false),0)
  {}
};

#if HAVE_ALUGRID
class ALUUnitSquare : public Dune::ALUSimplexGrid<2,2> 
{
public:
  ALUUnitSquare () : Dune::ALUSimplexGrid<2,2>("grids/2dsimplex.alu") {}
};

#endif //HAVE_ALUGRID


#if HAVE_ALBERTA
#  if ALBERTA_DIM == 2
class AlbertaLDomain : public Dune::AlbertaGrid<2,2>
{
public:
  AlbertaLDomain () : Dune::AlbertaGrid<2,2>("grids/ldomain.al") {}
};

class AlbertaUnitSquare : public Dune::AlbertaGrid<2,2>
{
public:
  AlbertaUnitSquare () : Dune::AlbertaGrid<2,2>("grids/2dgrid.al") {}
};

class AlbertaReentrantCorner : public Dune::GridPtr<Dune::AlbertaGrid<2,2> >
{
public:
  AlbertaReentrantCorner()
    : Dune::GridPtr<Dune::AlbertaGrid<2,2> >("grids/2dreentrantcorner.dgf")
  { }
};
#  endif //ALBERTA_DIM == 2
#endif //HAVE_ALBERTA

template<typename Grid>
class TriangulatedLDomainMaker {
  dune_static_assert(Grid::dimension == 2, "Dimension of grid must be 2");
  dune_static_assert(Grid::dimensionworld == 2, "Dimension of world must be 2");
public:
  static Dune::SmartPointer<Grid> create() {
    Dune::GridFactory<Grid> gf;

    Dune::FieldVector<typename Grid::ctype, 2> pos;
	pos[0] =-1.0;  pos[1] =-1.0; gf.insertVertex(pos);
	pos[0] = 0.0;  pos[1] =-1.0; gf.insertVertex(pos);
	pos[0] =-1.0;  pos[1] = 0.0; gf.insertVertex(pos);
	pos[0] = 0.0;  pos[1] = 0.0; gf.insertVertex(pos);
	pos[0] = 1.0;  pos[1] = 0.0; gf.insertVertex(pos);
	pos[0] =-1.0;  pos[1] = 1.0; gf.insertVertex(pos);
	pos[0] = 0.0;  pos[1] = 1.0; gf.insertVertex(pos);
	pos[0] = 1.0;  pos[1] = 1.0; gf.insertVertex(pos);

    Dune::GeometryType type; type.makeTriangle();
	std::vector<unsigned int> vid(3);
	vid[0] = 0;  vid[1] = 1;  vid[2] = 2; gf.insertElement(type, vid);
	vid[0] = 2;  vid[1] = 1;  vid[2] = 3; gf.insertElement(type, vid);
	vid[0] = 2;  vid[1] = 3;  vid[2] = 5; gf.insertElement(type, vid);
	vid[0] = 5;  vid[1] = 3;  vid[2] = 6; gf.insertElement(type, vid);
	vid[0] = 3;  vid[1] = 4;  vid[2] = 6; gf.insertElement(type, vid);
	vid[0] = 6;  vid[1] = 4;  vid[2] = 7; gf.insertElement(type, vid);

    return gf.createGrid();
  }
};

// Kuhn triangulation with 6 tets, does not contain unit tet, all tets have
// (0,7) as a common edge
template<typename Grid>
class KuhnTriangulatedUnitCubeMaker {
  dune_static_assert(Grid::dimension == 3, "Dimension of grid must be 3");
  dune_static_assert(Grid::dimensionworld == 3, "Dimension of world must be 3");
public:
  static Dune::SmartPointer<Grid> create() {
    Dune::GridFactory<Grid> gf;
    Dune::FieldVector<typename Grid::ctype, 3> pos;

    pos[0] = 0; pos[1] = 0; pos[2] = 0; gf.insertVertex(pos);
    pos[0] = 1; pos[1] = 0; pos[2] = 0; gf.insertVertex(pos);
    pos[0] = 0; pos[1] = 1; pos[2] = 0; gf.insertVertex(pos);
    pos[0] = 1; pos[1] = 1; pos[2] = 0; gf.insertVertex(pos);
    pos[0] = 0; pos[1] = 0; pos[2] = 1; gf.insertVertex(pos);
    pos[0] = 1; pos[1] = 0; pos[2] = 1; gf.insertVertex(pos);
    pos[0] = 0; pos[1] = 1; pos[2] = 1; gf.insertVertex(pos);
    pos[0] = 1; pos[1] = 1; pos[2] = 1; gf.insertVertex(pos);

    Dune::GeometryType type;
    type.makeTetrahedron();
    std::vector<unsigned int> vid(4);

    vid[0] = 0; vid[1] = 1; vid[2] = 3; vid[3] = 7; gf.insertElement(type, vid);
    vid[0] = 0; vid[1] = 1; vid[2] = 5; vid[3] = 7; gf.insertElement(type, vid);
    vid[0] = 0; vid[1] = 4; vid[2] = 5; vid[3] = 7; gf.insertElement(type, vid);
    vid[0] = 0; vid[1] = 4; vid[2] = 6; vid[3] = 7; gf.insertElement(type, vid);
    vid[0] = 0; vid[1] = 2; vid[2] = 6; vid[3] = 7; gf.insertElement(type, vid);
    vid[0] = 0; vid[1] = 2; vid[2] = 3; vid[3] = 7; gf.insertElement(type, vid);

    return gf.createGrid();
  }
};

#endif
