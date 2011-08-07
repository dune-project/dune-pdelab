// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=8 sw=4 sts=4:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include<iostream>
#include<map>
#include<string>
#include<vector>

#include <dune/common/array.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/mpihelper.hh>

#include<dune/grid/alugrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include<dune/pdelab/backend/istlvectorbackend.hh>
#include<dune/pdelab/finiteelementmap/mfdconstraints.hh>
#include<dune/pdelab/finiteelementmap/mimeticfem.hh>
#include<dune/pdelab/finiteelementmap/p0fem.hh>
#include<dune/pdelab/constraints/constraints.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/intersectionindexset.hh>

int main(int argc, char** argv)
{
    try
    {
        Dune::MPIHelper::instance(argc, argv);

        typedef Dune::ALUCubeGrid<3, 3> GridType;
        Dune::FieldVector<double, 3> lower_left(0.0);
        Dune::FieldVector<double, 3> upper_right(1.0);
        Dune::array<unsigned, 3> n;
        n.fill(1);
        Dune::shared_ptr<GridType> grid_ptr =
            Dune::StructuredGridFactory<GridType>::createCubeGrid(lower_left, upper_right, n);
        GridType& grid = *grid_ptr;
        typedef GridType::LeafGridView GV;
        const GV& gridview = grid.leafView();

        grid.globalRefine(1);
        grid.mark(1, *gridview.begin<0>());
        grid.preAdapt();
        grid.adapt();
        grid.postAdapt();

        // set up index set for intersections
        typedef Dune::PDELab::IntersectionIndexSet<GV> IIS;
        IIS iis(gridview);

        // make finite element maps
        typedef Dune::PDELab::P0LocalFiniteElementMap<double,double,3> CellFEM;
        Dune::GeometryType gt;
        gt.makeHexahedron();
        CellFEM cell_fem(gt);
        typedef Dune::PDELab::MimeticLocalFiniteElementMap<IIS,double,double,3> FaceFEM;
        FaceFEM face_fem(iis, Dune::GeometryType::cube);

        // make function spaces
        typedef Dune::PDELab::GridFunctionSpace<GV,CellFEM,
            Dune::PDELab::NoConstraints,Dune::PDELab::ISTLVectorBackend<1> > CellGFS;
        CellGFS cell_gfs(gridview, cell_fem);
        typedef Dune::PDELab::GridFunctionSpace<GV,FaceFEM,
            Dune::PDELab::MimeticConstraints,Dune::PDELab::ISTLVectorBackend<1>,
            Dune::PDELab::GridFunctionStaticSize<IIS> > FaceGFS;
        FaceGFS face_gfs(gridview, face_fem, iis);
        typedef Dune::PDELab::CompositeGridFunctionSpace
            <Dune::PDELab::GridFunctionSpaceLexicographicMapper,CellGFS,FaceGFS> GFS;
        GFS gfs(cell_gfs, face_gfs);

        typedef Dune::PDELab::LocalFunctionSpace<GFS> LFS;
        LFS lfs(gfs);
        typedef GV::Codim<0>::Iterator ElementIterator;
        for(ElementIterator it = gridview.begin<0>(), itend = gridview.end<0>(); it != itend; ++it)
            lfs.bind(*it);
    }
    catch (const Dune::Exception &e)
    {
        std::cerr << "Dune reported error: " << e << std::endl;
        throw;
    }
    catch (...)
    {
        std::cerr << "Unknown exception thrown." << std::endl;
        throw;
    }
}
