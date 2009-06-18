// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_STDDOMAINS_HH
#define DUNE_PDELAB_STDDOMAINS_HH

#include <string>
#include <sstream>

#include <dune/common/fvector.hh>
#include <dune/common/geometrytype.hh>
#include <dune/common/smartpointer.hh>
#include <dune/common/static_assert.hh>

#ifdef HAVE_ALBERTA
#include <dune/grid/albertagrid/gridfactory.hh>
#endif
#ifdef HAVE_ALUGRID
#include <dune/grid/alugrid.hh>
#endif
#ifdef HAVE_UG
#include <dune/grid/uggrid/uggridfactory.hh>
#endif

namespace Dune {
  namespace PDELab {

    // UnitTriangle

    template<typename Grid>
    class UnitTriangleMaker {
      dune_static_assert(Grid::dimension == 2, "Dimension of grid must be 2");
      dune_static_assert(Grid::dimensionworld == 2, "Dimension of world must be 2");
    public:
      static Dune::SmartPointer<Grid> create() {
        GridFactory<Grid> gf;
        FieldVector<typename Grid::ctype, 2> pos;

        pos[0] = 0; pos[1] = 0; gf.insertVertex(pos);
        pos[0] = 1; pos[1] = 0; gf.insertVertex(pos);
        pos[0] = 0; pos[1] = 1; gf.insertVertex(pos);

        GeometryType type;
        type.makeTriangle();
        std::vector<unsigned int> vid(3);

        vid[0] = 0; vid[1] = 1; vid[2] = 2; gf.insertElement(type, vid);

        return gf.createGrid();
      }
    };

#ifdef HAVE_ALUGRID
    template<>
    class UnitTriangleMaker<Dune::ALUSimplexGrid<2,2> > {
      typedef Dune::ALUSimplexGrid<2,2> Grid;
    public:
      static Dune::SmartPointer<Grid> create() {
        return new Grid("grids/2dtriangle.alu");
      }
    };
#endif // HAVE_ALUGRID

    // TriangulatedUnitSquare

    template<typename Grid>
    class TriangulatedUnitSquareMaker {
      dune_static_assert(Grid::dimension == 2, "Dimension of grid must be 2");
      dune_static_assert(Grid::dimensionworld == 2, "Dimension of world must be 2");
    public:
      static Dune::SmartPointer<Grid> create() {
        Dune::GridFactory<Grid> gf;
        Dune::FieldVector<typename Grid::ctype, 2> pos;

        pos[0] = 0; pos[1] = 0; gf.insertVertex(pos);
        pos[0] = 1; pos[1] = 0; gf.insertVertex(pos);
        pos[0] = 0; pos[1] = 1; gf.insertVertex(pos);
        pos[0] = 1; pos[1] = 1; gf.insertVertex(pos);

        Dune::GeometryType type;
        type.makeTriangle();
        std::vector<unsigned int> vid(3);

        vid[0] = 0; vid[1] = 1; vid[2] = 2; gf.insertElement(type, vid);
        vid[0] = 1; vid[1] = 2; vid[2] = 3; gf.insertElement(type, vid);

        return gf.createGrid();
      }
    };

#ifdef HAVE_ALUGRID
    template<>
    class TriangulatedUnitSquareMaker<Dune::ALUSimplexGrid<2,2> > {
      typedef Dune::ALUSimplexGrid<2,2> Grid;
    public:
      static Dune::SmartPointer<Grid> create() {
        return new Grid("grids/2dsimplex.alu");
      }
    };
#endif // HAVE_ALUGRID

    // UnitTetrahedron

    template<typename Grid>
    class UnitTetrahedronMaker {
      dune_static_assert(Grid::dimension == 3, "Dimension of grid must be 3");
      dune_static_assert(Grid::dimensionworld == 3, "Dimension of world must be 3");
    public:
      static Dune::SmartPointer<Grid> create() {
        Dune::GridFactory<Grid> gf;
        Dune::FieldVector<typename Grid::ctype, 3> pos;

        pos[0] = 0; pos[1] = 0; pos[2] = 0; gf.insertVertex(pos);
        pos[0] = 1; pos[1] = 0; pos[2] = 0; gf.insertVertex(pos);
        pos[0] = 0; pos[1] = 1; pos[2] = 0; gf.insertVertex(pos);
        pos[0] = 0; pos[1] = 0; pos[2] = 1; gf.insertVertex(pos);

        Dune::GeometryType type;
        type.makeTetrahedron();
        std::vector<unsigned int> vid(4);

        vid[0] = 0; vid[1] = 1; vid[2] = 2; vid[3] = 3; gf.insertElement(type, vid);

        return gf.createGrid();
      }
    };

    // TriangulatedUnitCube

    template<typename Grid>
    class TriangulatedUnitCubeMaker {
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

        return gf.createGrid();
      }
    };

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_STDDOMAINS_HH
