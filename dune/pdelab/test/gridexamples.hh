// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_TEST_GRIDEXAMPLES_HH
#define DUNE_PDELAB_TEST_GRIDEXAMPLES_HH

#include <memory>
#include <vector>

#include <dune/common/fvector.hh>

#include <dune/geometry/type.hh>

#include<dune/grid/yaspgrid.hh>
#include<dune/grid/common/gridfactory.hh>
#if HAVE_ALBERTA
#include<dune/grid/albertagrid.hh>
#include <dune/grid/albertagrid/dgfparser.hh>
#include <dune/grid/albertagrid/gridfactory.hh>
#endif
#if HAVE_UG
#include <dune/grid/uggrid/uggridfactory.hh>
#endif
#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif

class YaspUnitSquare : public Dune::YaspGrid<2>
{
public:

  YaspUnitSquare () : Dune::YaspGrid<2>(Dune::FieldVector<double,2>(1.0),
                                        {{1,1}},
                                        std::bitset<2>(false),
                                        0)
  {}
};

#if HAVE_DUNE_ALUGRID
class ALUUnitSquare : public Dune::ALUGrid<2,2,Dune::simplex,Dune::nonconforming>
{
public:
  ALUUnitSquare () : Dune::ALUGrid<2,2,Dune::simplex,Dune::nonconforming>(GRIDSDIR "/2dsimplex.alu") {}
};

#endif //HAVE_DUNE_ALUGRID


#if HAVE_ALBERTA
#  if ALBERTA_DIM == 2
class AlbertaLDomain : public Dune::AlbertaGrid<2,2>
{
public:
  AlbertaLDomain () : Dune::AlbertaGrid<2,2>(GRIDSDIR "/ldomain.al") {}
};

class AlbertaUnitSquare : public Dune::AlbertaGrid<2,2>
{
public:
  AlbertaUnitSquare () : Dune::AlbertaGrid<2,2>(GRIDSDIR "/2dgrid.al") {}
};

class AlbertaReentrantCorner : public Dune::GridPtr<Dune::AlbertaGrid<2,2> >
{
public:
  AlbertaReentrantCorner()
    : Dune::GridPtr<Dune::AlbertaGrid<2,2> >(GRIDSDIR "/2dreentrantcorner.dgf")
  { }
};
#  endif //ALBERTA_DIM == 2
#endif //HAVE_ALBERTA

template<typename Grid>
class TriangulatedLDomainMaker {
  static_assert(Grid::dimension == 2, "Dimension of grid must be 2");
  static_assert(Grid::dimensionworld == 2, "Dimension of world must be 2");
public:
  static std::unique_ptr<Grid> create() {
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

    auto type = Dune::GeometryTypes::triangle;
	std::vector<unsigned int> vid(3);
	vid[0] = 0;  vid[1] = 1;  vid[2] = 2; gf.insertElement(type, vid);
	vid[0] = 2;  vid[1] = 1;  vid[2] = 3; gf.insertElement(type, vid);
	vid[0] = 2;  vid[1] = 3;  vid[2] = 5; gf.insertElement(type, vid);
	vid[0] = 5;  vid[1] = 3;  vid[2] = 6; gf.insertElement(type, vid);
	vid[0] = 3;  vid[1] = 4;  vid[2] = 6; gf.insertElement(type, vid);
	vid[0] = 6;  vid[1] = 4;  vid[2] = 7; gf.insertElement(type, vid);

    return std::unique_ptr<Grid>(gf.createGrid());
  }
};

//////////////////////////////////////////////////////////////////////
//
// UnitTriangle
//

template<typename Grid>
class UnitTriangleMaker {
  static_assert(Grid::dimension == 2, "Dimension of grid must be 2");
  static_assert(Grid::dimensionworld == 2, "Dimension of world must be 2");
public:
  static std::unique_ptr<Grid> create() {
    Dune::GridFactory<Grid> gf;
    Dune::FieldVector<typename Grid::ctype, 2> pos;

    pos[0] = 0; pos[1] = 0; gf.insertVertex(pos);
    pos[0] = 1; pos[1] = 0; gf.insertVertex(pos);
    pos[0] = 0; pos[1] = 1; gf.insertVertex(pos);

    auto type = Dune::GeometryTypes::triangle;
    std::vector<unsigned int> vid(3);

    vid[0] = 0; vid[1] = 1; vid[2] = 2; gf.insertElement(type, vid);

    return std::unique_ptr<Grid>(gf.createGrid());
  }
};

#if HAVE_DUNE_ALUGRID
template<>
class UnitTriangleMaker<Dune::ALUGrid<2,2,Dune::simplex,Dune::nonconforming> > {
  typedef Dune::ALUGrid<2,2,Dune::simplex,Dune::nonconforming> Grid;
public:
  static std::unique_ptr<Grid> create() {
    return std::unique_ptr<Grid>(new Grid(GRIDSDIR "/2dtriangle.alu"));
  }
};
#endif // HAVE_DUNE_ALUGRID

//////////////////////////////////////////////////////////////////////
//
// TriangulatedUnitSquare
//

template<typename Grid>
class TriangulatedUnitSquareMaker {
  static_assert(Grid::dimension == 2, "Dimension of grid must be 2");
  static_assert(Grid::dimensionworld == 2, "Dimension of world must be 2");
public:
  static std::unique_ptr<Grid> create() {
    Dune::GridFactory<Grid> gf;
    Dune::FieldVector<typename Grid::ctype, 2> pos;

    pos[0] = 0; pos[1] = 0; gf.insertVertex(pos);
    pos[0] = 1; pos[1] = 0; gf.insertVertex(pos);
    pos[0] = 0; pos[1] = 1; gf.insertVertex(pos);
    pos[0] = 1; pos[1] = 1; gf.insertVertex(pos);

    auto type = Dune::GeometryTypes::triangle;
    std::vector<unsigned int> vid(3);

    vid[0] = 0; vid[1] = 1; vid[2] = 2; gf.insertElement(type, vid);
    vid[0] = 1; vid[1] = 2; vid[2] = 3; gf.insertElement(type, vid);

    return std::unique_ptr<Grid>(gf.createGrid());
  }
};

#if HAVE_DUNE_ALUGRID
template<>
class TriangulatedUnitSquareMaker<Dune::ALUGrid<2,2,Dune::simplex,Dune::nonconforming> > {
  typedef Dune::ALUGrid<2,2,Dune::simplex,Dune::nonconforming> Grid;
public:
  static std::unique_ptr<Grid> create() {
    return std::unique_ptr<Grid>(new Grid(GRIDSDIR "/2dsimplex.alu"));
  }
};
#endif //HAVE_DUNE_ALUGRID

//////////////////////////////////////////////////////////////////////
//
// UnitTetrahedron
//

template<typename Grid>
class UnitTetrahedronMaker {
  static_assert(Grid::dimension == 3, "Dimension of grid must be 3");
  static_assert(Grid::dimensionworld == 3, "Dimension of world must be 3");
public:
  static std::unique_ptr<Grid> create() {
    Dune::GridFactory<Grid> gf;
    Dune::FieldVector<typename Grid::ctype, 3> pos;

    pos[0] = 0; pos[1] = 0; pos[2] = 0; gf.insertVertex(pos);
    pos[0] = 1; pos[1] = 0; pos[2] = 0; gf.insertVertex(pos);
    pos[0] = 0; pos[1] = 1; pos[2] = 0; gf.insertVertex(pos);
    pos[0] = 0; pos[1] = 0; pos[2] = 1; gf.insertVertex(pos);

    auto type = Dune::GeometryTypes::tetrahedron;
    std::vector<unsigned int> vid(4);

    vid[0] = 0; vid[1] = 1; vid[2] = 2; vid[3] = 3; gf.insertElement(type, vid);

    return std::unique_ptr<Grid>(gf.createGrid());
  }
};

//////////////////////////////////////////////////////////////////////
//
// TriangulatedUnitCube
//

// Minimal triangulation with 5 tets, does contain unit tet
// AlbertaSimplexGrid<3,3> cannot refine this, see Flyspry#569
template<typename Grid>
class TriangulatedUnitCubeMaker {
  static_assert(Grid::dimension == 3, "Dimension of grid must be 3");
  static_assert(Grid::dimensionworld == 3, "Dimension of world must be 3");
public:
  static std::unique_ptr<Grid> create() {
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

    auto type = Dune::GeometryTypes::tetrahedron;
    std::vector<unsigned int> vid(4);

    // tet at vertex 0
    vid[0] = 0; vid[1] = 1; vid[2] = 2; vid[3] = 4; gf.insertElement(type, vid);
    // tet at vertex 3
    vid[0] = 1; vid[1] = 2; vid[2] = 3; vid[3] = 7; gf.insertElement(type, vid);
    // central tet
    vid[0] = 1; vid[1] = 2; vid[2] = 4; vid[3] = 7; gf.insertElement(type, vid);
    // tet at vertex 5
    vid[0] = 1; vid[1] = 4; vid[2] = 5; vid[3] = 7; gf.insertElement(type, vid);
    // tet at vertex 6
    vid[0] = 2; vid[1] = 4; vid[2] = 6; vid[3] = 7; gf.insertElement(type, vid);

    return std::unique_ptr<Grid>(gf.createGrid());
  }
};

#if HAVE_ALBERTA
# if ALBERTA_DIM == 3
#  ifndef ALLOW_ALBERTA_MINIMAL_TRIANGULATED_CUBE
// AlbertaSimplexGrid<3,3> cannot refine the minimal triangulated cube, see
// Flyspry#569.  If you want to use it nevertheless, define
// ALLOW_ALBERTA_MINIMAL_TRIANGULATED_CUBE before including gridexamples.hh.
// specialize the template to make any attempt to use the create() method fail.
template<>
class TriangulatedUnitCubeMaker<Dune::AlbertaGrid<3,3> >
{};
#  endif //ALLOW_ALBERTA_MINIMAL_TRIANGULATED_CUBE
# endif //ALBERTA_DIM == 3
#endif //HAVE_ALBERTA

//////////////////////////////////////////////////////////////////////
//
// KuhnTriangulatedUnitCubeMaker
//

// Kuhn triangulation with 6 tets, does not contain unit tet, all tets have
// (0,7) as a common edge
template<typename Grid>
class KuhnTriangulatedUnitCubeMaker {
  static_assert(Grid::dimension == 3, "Dimension of grid must be 3");
  static_assert(Grid::dimensionworld == 3, "Dimension of world must be 3");
public:
  static std::unique_ptr<Grid> create() {
    Dune::GridFactory<Grid> gf;

    int fake_argc = 0;
    char **fake_argv = NULL;

    if(Dune::MPIHelper::instance(fake_argc, fake_argv).rank() == 0) {
      Dune::FieldVector<typename Grid::ctype, 3> pos;

      pos[0] = 0; pos[1] = 0; pos[2] = 0; gf.insertVertex(pos);
      pos[0] = 1; pos[1] = 0; pos[2] = 0; gf.insertVertex(pos);
      pos[0] = 0; pos[1] = 1; pos[2] = 0; gf.insertVertex(pos);
      pos[0] = 1; pos[1] = 1; pos[2] = 0; gf.insertVertex(pos);
      pos[0] = 0; pos[1] = 0; pos[2] = 1; gf.insertVertex(pos);
      pos[0] = 1; pos[1] = 0; pos[2] = 1; gf.insertVertex(pos);
      pos[0] = 0; pos[1] = 1; pos[2] = 1; gf.insertVertex(pos);
      pos[0] = 1; pos[1] = 1; pos[2] = 1; gf.insertVertex(pos);

      auto type = Dune::GeometryTypes::tetrahedron;
      std::vector<unsigned int> vid(4);

      vid[0] = 0; vid[1] = 1; vid[2] = 3; vid[3] = 7; gf.insertElement(type, vid);
      vid[0] = 0; vid[1] = 1; vid[2] = 5; vid[3] = 7; gf.insertElement(type, vid);
      vid[0] = 0; vid[1] = 4; vid[2] = 5; vid[3] = 7; gf.insertElement(type, vid);
      vid[0] = 0; vid[1] = 4; vid[2] = 6; vid[3] = 7; gf.insertElement(type, vid);
      vid[0] = 0; vid[1] = 2; vid[2] = 6; vid[3] = 7; gf.insertElement(type, vid);
      vid[0] = 0; vid[1] = 2; vid[2] = 3; vid[3] = 7; gf.insertElement(type, vid);
    }

    std::unique_ptr<Grid> gp(gf.createGrid());
    gp->loadBalance();
    return gp;
  }
};

#endif // DUNE_PDELAB_TEST_GRIDEXAMPLES_HH
