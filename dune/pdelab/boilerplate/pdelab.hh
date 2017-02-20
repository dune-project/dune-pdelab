// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
#ifndef DUNE_PDELAB_BOILERPLATE_PDELAB_HH
#define DUNE_PDELAB_BOILERPLATE_PDELAB_HH

/** \file
    \brief Provide some classes to reduce boiler plate code in pdelab applications

    These classes are experimental !

    To see examples how they might simplify your life, have a look at the
    dune-pdelab-howto module, in particular at the examples in
    dune-pdelab-howto/src/boilerplatetutorial/.
*/

// first of all we include a lot of dune grids and pdelab files
#include <iostream>
#include <memory>

#include <dune/common/parallel/mpihelper.hh> // include mpi helper class
#include <dune/common/parametertreeparser.hh>
#include <dune/common/classname.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/onedgrid.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/yaspgrid.hh>
#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif
#if HAVE_ALBERTA
#include<dune/grid/albertagrid.hh>
#include <dune/grid/albertagrid/dgfparser.hh>
#endif
#if HAVE_UG
#include<dune/grid/uggrid.hh>
#endif
#if HAVE_DUNE_ALUGRID
#include<dune/alugrid/grid.hh>
#include <dune/alugrid/dgf.hh>
#endif
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/io/file/gmshreader.hh>

#include <dune/istl/bvector.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/solvercategory.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/io.hh>

#include <dune/istl/paamg/amg.hh>
#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/common/functionutilities.hh>
#include <dune/pdelab/common/vtkexport.hh>
#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/constraints/hangingnode.hh>
#include <dune/pdelab/constraints/p0.hh>
#include <dune/pdelab/constraints/p0ghost.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/gridoperator/onestep.hh>
#include <dune/pdelab/stationary/linearproblem.hh>
#include <dune/pdelab/finiteelementmap/pkfem.hh>
#include <dune/pdelab/finiteelementmap/p0fem.hh>
#include <dune/pdelab/finiteelementmap/opbfem.hh>
#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/pdelab/finiteelementmap/qkdg.hh>
#include <dune/pdelab/adaptivity/adaptivity.hh>
#include <dune/pdelab/instationary/onestep.hh>
#include <dune/pdelab/common/instationaryfilenamehelper.hh>
#include <dune/pdelab/newton/newton.hh>

namespace Dune {
    namespace PDELab {

        // make grids
        template<typename T>
        class StructuredGrid
        {
        public:
            // export types
            typedef T Grid;
            typedef typename T::ctype ctype;
            static const int dim = T::dimension;
            static const int dimworld = T::dimensionworld;

            // constructors
            StructuredGrid (Dune::GeometryType::BasicType meshtype, unsigned int cells)
            {
                FieldVector<ctype,dimworld> lowerLeft(0.0);
                FieldVector<ctype,dimworld> upperRight(1.0);
                array<unsigned int,dim> elements; elements.fill(cells);

                StructuredGridFactory<T> factory;

                if (meshtype==Dune::GeometryType::cube)
                    gridp = factory.createCubeGrid(lowerLeft,upperRight,elements);
                else if (meshtype==Dune::GeometryType::simplex)
                    gridp = factory.createSimplexGrid(lowerLeft,upperRight,elements);
                else
                    {
                        DUNE_THROW(GridError, className<StructuredGrid>()
                                   << "::StructuredGrid(): grid type must be simplex or cube ");
                    }
            }


            StructuredGrid (Dune::GeometryType::BasicType meshtype,
                            array<double,dimworld> lower_left, array<double,dimworld> upper_right,
                            array<unsigned int,dim> cells)
            {
                FieldVector<ctype,dimworld> lowerLeft;
                FieldVector<ctype,dimworld> upperRight;
                array<unsigned int,dim> elements;

                // copy data to correct types for StructuredGridFactory
                for (size_t i=0; i<dimworld; i++)
                    {
                        lowerLeft[i] = lower_left[i];
                        upperRight[i] = upper_right[i];
                    }
                for (size_t i=0; i<dim; i++)
                    {
                        elements[i] = cells[i];
                    }

                StructuredGridFactory<T> factory;

                if (meshtype==Dune::GeometryType::cube)
                    gridp = factory.createCubeGrid(lowerLeft,upperRight,elements);
                else if (meshtype==Dune::GeometryType::simplex)
                    gridp = factory.createSimplexGrid(lowerLeft,upperRight,elements);
                else
                    {
                        DUNE_THROW(GridError, className<StructuredGrid>()
                                   << "::StructuredGrid(): grid type must be simplex or cube ");
                    }
            }

            // return shared pointer
            std::shared_ptr<T> getSharedPtr ()
            {
                return gridp;
            }

            // return grid reference
            T& getGrid ()
            {
                return *gridp;
            }

            // return grid reference const version
            const T& getGrid () const
            {
                return *gridp;
            }

            T& operator*()
            {
                return *gridp;
            }

            T* operator->()
            {
                return gridp.operator->();
            }

            const T& operator*() const
            {
                return *gridp;
            }

            const T* operator->() const
            {
                return gridp.operator->();
            }


        private:
            std::shared_ptr<T> gridp; // hold a shared pointer to a grid
        };

        // specialization for yaspgrid; treats paralle case right
        template<int dim>
        class StructuredGrid<YaspGrid<dim> >
        {
        public:

            // export types
            typedef YaspGrid<dim> Grid;
            typedef typename Grid::ctype ctype;
            static const int dimworld = Grid::dimensionworld;

            // simple constructor for the unit cube
            StructuredGrid (Dune::GeometryType::BasicType meshtype, unsigned int cells, int overlap=1)
            {
                // check element type
                if (meshtype!=Dune::GeometryType::cube)
                    std::cout << "StructuredGrid(): element type " << meshtype << " is ignored" << std::endl;

                // copy data to correct types for YaspGrid
                Dune::FieldVector<double,dimworld> L(1.0);
                std::array<int,dimworld> N(Dune::fill_array<int,dimworld>(cells));
                std::bitset<dimworld> B(false);

                // instantiate the grid
                gridp = std::shared_ptr<Grid>(new Grid(L,N,B,overlap,Dune::MPIHelper::getCollectiveCommunication()));
            }

            // constructor with sizes given
            StructuredGrid (Dune::GeometryType::BasicType meshtype,
                            array<double,dimworld> lower_left, array<double,dimworld> upper_right,
                            array<unsigned int,dim> cells, int overlap=1)
            {
                // check that lower right corner is the origin
                for(int d = 0; d < dimworld; ++d)
                    if(std::abs(lower_left[d]) > std::abs(upper_right[d])*1e-10)
                        DUNE_THROW(GridError, className<StructuredGrid>()
                                   << "::createCubeGrid(): The lower coordinates "
                                   "must be at the origin for YaspGrid.");

                // check element type
                if (meshtype!=Dune::GeometryType::cube)
                    std::cout << "StructuredGrid(): element type " << meshtype << " is ignored" << std::endl;

                // copy data to correct types for YaspGrid
                Dune::FieldVector<double,dimworld> L;
                std::array<int,dimworld> N;
                std::bitset<dimworld> B(false);
                for (size_t i=0; i<dimworld; i++)
                    {
                        L[i] = upper_right[i];
                        N[i] = cells[i];
                    }

                // instantiate the grid
                gridp = std::shared_ptr<Grid>(new Grid(L,N,B,overlap,Dune::MPIHelper::getCollectiveCommunication()));
            }

            // constructor with periodicity argument
            StructuredGrid (Dune::GeometryType::BasicType meshtype,
                            array<double,dimworld> lower_left, array<double,dimworld> upper_right,
                            array<unsigned int,dim> cells, array<bool,dim> periodic, int overlap=1)
            {
                // check that lower right corner is the origin
                for(int d = 0; d < dimworld; ++d)
                    if(std::abs(lower_left[d]) > std::abs(upper_right[d])*1e-10)
                        DUNE_THROW(GridError, className<StructuredGrid>()
                                   << "::createCubeGrid(): The lower coordinates "
                                   "must be at the origin for YaspGrid.");

                // check element type
                if (meshtype!=Dune::GeometryType::cube)
                    std::cout << "StructuredGrid(): element type " << meshtype << " is ignored" << std::endl;

                // copy data to correct types for YaspGrid
                Dune::FieldVector<double,dimworld> L;
                std::array<int,dimworld> N;
                std::bitset<dimworld> B(false);
                for (size_t i=0; i<dimworld; i++)
                    {
                        L[i] = upper_right[i];
                        N[i] = cells[i];
                        B[i] = periodic[i];
                    }

                // instantiate the grid
                gridp = std::shared_ptr<Grid>(new Grid(L,N,B,overlap,Dune::MPIHelper::getCollectiveCommunication()));
            }

            // return shared pointer
            std::shared_ptr<Grid> getSharedPtr ()
            {
                return gridp;
            }

            // return grid reference
            Grid& getGrid ()
            {
                return *gridp;
            }

            // return grid reference const version
            const Grid& getGrid () const
            {
                return *gridp;
            }

            Grid& operator*()
            {
                return *gridp;
            }

            Grid* operator->()
            {
                return gridp.operator->();
            }

            const Grid& operator*() const
            {
                return *gridp;
            }

            const Grid* operator->() const
            {
                return gridp.operator->();
            }

        private:
            std::shared_ptr<Grid> gridp; // hold a shared pointer to a grid
        };

        // unstructured grid read from gmsh file
        template<typename T>
        class UnstructuredGrid
        {
        public:
            // export types
            typedef T Grid;
            typedef typename T::ctype ctype;
            static const int dim = T::dimension;
            static const int dimworld = T::dimensionworld;

            // constructors
            UnstructuredGrid (std::string filename, bool verbose = true, bool insert_boundary_segments=true)
            {
                Dune::GridFactory<T> factory;
                Dune::GmshReader<T>::read(factory,filename,verbose,insert_boundary_segments);
                gridp = std::shared_ptr<T>(factory.createGrid());
            }

            // return shared pointer
            std::shared_ptr<T> getSharedPtr ()
            {
                return gridp;
            }

            // return grid reference
            T& getGrid ()
            {
                return *gridp;
            }

            // return grid reference const version
            const T& getGrid () const
            {
                return *gridp;
            }

            T& operator*()
            {
                return *gridp;
            }

            T* operator->()
            {
                return gridp.operator->();
            }

            const T& operator*() const
            {
                return *gridp;
            }

            const T* operator->() const
            {
                return gridp.operator->();
            }

        private:
            std::shared_ptr<T> gridp; // hold a shared pointer to a grid
        };


        //============================================================================
        // Continuous Lagrange Finite Element Space
        //============================================================================

        // finite element map base template
        template<typename GV, typename C, typename R, unsigned int degree, unsigned int dim, Dune::GeometryType::BasicType gt>
        class CGFEMBase
        {};

        template<typename GV, typename C, typename R, unsigned int degree, unsigned int dim>
        class CGFEMBase<GV,C,R,degree,dim,Dune::GeometryType::simplex>
        {
        public:
            typedef PkLocalFiniteElementMap<GV,C,R,degree> FEM;

            CGFEMBase (const GV& gridview)
            {
                femp = std::shared_ptr<FEM>(new FEM(gridview));
            }

            FEM& getFEM() {return *femp;}
            const FEM& getFEM() const {return *femp;}

        private:
            std::shared_ptr<FEM> femp;
        };

        template<typename GV, typename C, typename R, unsigned int degree, unsigned int dim>
        class CGFEMBase<GV,C,R,degree,dim,Dune::GeometryType::cube>
        {
        public:
            typedef QkLocalFiniteElementMap<GV,C,R,degree> FEM;

            CGFEMBase (const GV& gridview)
            {
                femp = std::shared_ptr<FEM>(new FEM(gridview));
            }

            FEM& getFEM() {return *femp;}
            const FEM& getFEM() const {return *femp;}

        private:
            std::shared_ptr<FEM> femp;
        };

        //============================================================================

        // define enumeration type that differentiate conforming and nonconforming meshes
        enum MeshType {
            conforming,
            nonconforming
        };

        // constraints base template
        template<typename Grid, unsigned int degree, Dune::GeometryType::BasicType gt, MeshType mt, SolverCategory::Category st, typename BCType, typename GV = typename Grid::LeafGridView>
        class CGCONBase
        {};

        template<typename Grid, typename BCType, typename GV>
        class CGCONBase<Grid,1,Dune::GeometryType::simplex,MeshType::nonconforming,SolverCategory::sequential,BCType,GV>
        {
        public:
            typedef HangingNodesDirichletConstraints<Grid,HangingNodesConstraintsAssemblers::SimplexGridP1Assembler,BCType> CON;

            CGCONBase (Grid& grid, const BCType& bctype, const GV& gv)
            {
                conp = std::shared_ptr<CON>(new CON(grid,true,bctype));
            }

            CGCONBase (Grid& grid, const BCType& bctype)
            {
                conp = std::shared_ptr<CON>(new CON(grid,true,bctype));
            }

            template<typename GFS>
            void postGFSHook (const GFS& gfs) {}
            CON& getCON() {return *conp;}
            const CON& getCON() const {return *conp;}
            template<typename GFS, typename DOF>
            void make_consistent (const GFS& gfs, DOF& x) const {}
        private:
            std::shared_ptr<CON> conp;
        };

        template<typename Grid, typename BCType, typename GV>
        class CGCONBase<Grid,1,Dune::GeometryType::cube,MeshType::nonconforming,SolverCategory::sequential,BCType,GV>
        {
        public:
            typedef HangingNodesDirichletConstraints<Grid,HangingNodesConstraintsAssemblers::CubeGridQ1Assembler,BCType> CON;

            CGCONBase (Grid& grid, const BCType& bctype, const GV& gv)
            {
                conp = std::shared_ptr<CON>(new CON(grid,true,bctype));
            }

            CGCONBase (Grid& grid, const BCType& bctype)
            {
                conp = std::shared_ptr<CON>(new CON(grid,true,bctype));
            }

            template<typename GFS>
            void postGFSHook (const GFS& gfs) {}
            CON& getCON() {return *conp;}
            const CON& getCON() const {return *conp;}
            template<typename GFS, typename DOF>
            void make_consistent (const GFS& gfs, DOF& x) const {}
        private:
            std::shared_ptr<CON> conp;
        };

        template<typename Grid, unsigned int degree, Dune::GeometryType::BasicType gt,typename BCType, typename GV>
        class CGCONBase<Grid,degree,gt,MeshType::conforming,SolverCategory::sequential,BCType,GV>
        {
        public:
            typedef ConformingDirichletConstraints CON;

            CGCONBase (Grid& grid, const BCType& bctype, const GV& gv)
            {
                conp = std::shared_ptr<CON>(new CON());
            }

            CGCONBase (Grid& grid, const BCType& bctype)
            {
                conp = std::shared_ptr<CON>(new CON());
            }

            template<typename GFS>
            void postGFSHook (const GFS& gfs) {}
            CON& getCON() {return *conp;}
            const CON& getCON() const {return *conp;}
            template<typename GFS, typename DOF>
            void make_consistent (const GFS& gfs, DOF& x) const {}
        private:
            std::shared_ptr<CON> conp;
        };

        template<typename Grid, unsigned int degree, Dune::GeometryType::BasicType gt,typename BCType, typename GV>
        class CGCONBase<Grid,degree,gt,MeshType::conforming,SolverCategory::overlapping,BCType,GV>
        {
        public:
            typedef OverlappingConformingDirichletConstraints CON;

            CGCONBase (Grid& grid, const BCType& bctype, const GV& gv)
            {
                conp = std::shared_ptr<CON>(new CON());
            }

            CGCONBase (Grid& grid, const BCType& bctype)
            {
                conp = std::shared_ptr<CON>(new CON());
            }

            template<typename GFS>
            void postGFSHook (const GFS& gfs) {}
            CON& getCON() {return *conp;}
            const CON& getCON() const {return *conp;}
            template<typename GFS, typename DOF>
            void make_consistent (const GFS& gfs, DOF& x) const
            {
                // make vector consistent; this is needed for all overlapping solvers
                ISTL::ParallelHelper<GFS> helper(gfs);
                helper.maskForeignDOFs(Backend::native(x));
                Dune::PDELab::AddDataHandle<GFS,DOF> adddh(gfs,x);
                if (gfs.gridView().comm().size()>1)
                    gfs.gridView().communicate(adddh,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);
            }
        private:
            std::shared_ptr<CON> conp;
        };

        template<typename Grid, unsigned int degree, Dune::GeometryType::BasicType gt,typename BCType, typename GV>
        class CGCONBase<Grid,degree,gt,MeshType::conforming,SolverCategory::nonoverlapping,BCType,GV>
        {
        public:
            using CON = Dune::PDELab::ConformingDirichletConstraints;
            CGCONBase (Grid& grid, const BCType& bctype)
            {
                conp = std::shared_ptr<CON>(new CON);
            }

            template<typename GFS>
            CON& getCON() {return *conp;}
            const CON& getCON() const {return *conp;}
            template<typename GFS, typename DOF>
            void make_consistent (const GFS& gfs, DOF& x) const {}
        private:
            std::shared_ptr<CON> conp;
        };


        // continuous Lagrange finite elements
        template<typename T, typename N, unsigned int degree, typename BCType,
                 Dune::GeometryType::BasicType gt, MeshType mt, SolverCategory::Category st = SolverCategory::sequential,
                 typename VBET=ISTL::VectorBackend<> >
        class CGSpace {
        public:

            // export types
            typedef T Grid;
            typedef typename T::LeafGridView GV;
            typedef typename T::ctype ctype;
            static const int dim = T::dimension;
            static const int dimworld = T::dimensionworld;

            typedef CGFEMBase<GV,ctype,N,degree,dim,gt> FEMB;
            typedef CGCONBase<Grid,degree,gt,mt,st,BCType> CONB;

            typedef typename FEMB::FEM FEM;
            typedef typename CONB::CON CON;

            typedef VBET VBE;
            typedef GridFunctionSpace<GV,FEM,CON,VBE> GFS;

            typedef N NT;
            using DOF = Backend::Vector<GFS,N>;
            typedef Dune::PDELab::DiscreteGridFunction<GFS,DOF> DGF;
            typedef typename GFS::template ConstraintsContainer<N>::Type CC;
            typedef VTKGridFunctionAdapter<DGF> VTKF;

            // constructor making the grid function space an all that is needed
            CGSpace (Grid& grid, const BCType& bctype)
                : gv(grid.leafGridView()), femb(gv), conb(grid,bctype)
            {
                gfsp = std::shared_ptr<GFS>(new GFS(gv,femb.getFEM(),conb.getCON()));
                gfsp->name("cgspace");
                // initialize ordering
                gfsp->update();
                conb.postGFSHook(*gfsp);
                ccp = std::shared_ptr<CC>(new CC());
            }

            FEM& getFEM()
            {
                return femb.getFEM();
            }

            const FEM& getFEM() const
            {
                return femb.getFEM();
            }

            // return gfs reference
            GFS& getGFS ()
            {
                return *gfsp;
            }

            // return gfs reference const version
            const GFS& getGFS () const
            {
                return *gfsp;
            }

            // return gfs reference
            CC& getCC ()
            {
                return *ccp;
            }

            // return gfs reference const version
            const CC& getCC () const
            {
                return *ccp;
            }

            void assembleConstraints (const BCType& bctype)
            {
                ccp->clear();
                constraints(bctype,*gfsp,*ccp);
            }

            void clearConstraints ()
            {
                ccp->clear();
            }

            void setConstrainedDOFS (DOF& x, NT nt) const
            {
                set_constrained_dofs(*ccp,nt,x);
                conb.make_consistent(*gfsp,x);
            }

            void setNonConstrainedDOFS (DOF& x, NT nt) const
            {
                set_nonconstrained_dofs(*ccp,nt,x);
                conb.make_consistent(*gfsp,x);
            }

            void copyConstrainedDOFS (const DOF& xin, DOF& xout) const
            {
                copy_constrained_dofs(*ccp,xin,xout);
                conb.make_consistent(*gfsp,xout);
            }

            void copyNonConstrainedDOFS (const DOF& xin, DOF& xout) const
            {
                copy_nonconstrained_dofs(*ccp,xin,xout);
                conb.make_consistent(*gfsp,xout);
            }

        private:
            GV gv; // need this object here because FEM and GFS store a const reference !!
            FEMB femb;
            CONB conb;
            std::shared_ptr<GFS> gfsp;
            std::shared_ptr<CC> ccp;
        };

        // template specialization for nonoverlapping case
        template<typename T, typename N, unsigned int degree, typename BCType,
                 Dune::GeometryType::BasicType gt, MeshType mt,
                 typename VBET>
        class CGSpace<T, N, degree, BCType, gt, mt, SolverCategory::nonoverlapping, VBET> {
        public:

            // export types
            typedef T Grid;
            typedef typename T::LeafGridView GV;
            typedef typename T::ctype ctype;
            typedef typename Dune::PDELab::NonOverlappingEntitySet<GV> ES;
            static const int dim = T::dimension;
            static const int dimworld = T::dimensionworld;

            typedef CGFEMBase<ES,ctype,N,degree,dim,gt> FEMB;
            typedef CGCONBase<Grid,degree,gt,mt,SolverCategory::nonoverlapping,BCType> CONB;

            typedef typename FEMB::FEM FEM;
            typedef typename CONB::CON CON;

            typedef VBET VBE;
            typedef GridFunctionSpace<ES,FEM,CON,VBE> GFS;

            typedef N NT;
            using DOF = Backend::Vector<GFS,N>;
            typedef Dune::PDELab::DiscreteGridFunction<GFS,DOF> DGF;
            typedef typename GFS::template ConstraintsContainer<N>::Type CC;
            typedef VTKGridFunctionAdapter<DGF> VTKF;

            // constructor making the grid function space an all that is needed
            CGSpace (Grid& grid, const BCType& bctype)
                : gv(grid.leafGridView()), es(gv), femb(es), conb(grid,bctype)
            {
                gfsp = std::shared_ptr<GFS>(new GFS(es,femb.getFEM(),conb.getCON()));
                gfsp->name("cgspace");
                // initialize ordering
                gfsp->update();
                // conb.postGFSHook(*gfsp);
                ccp = std::shared_ptr<CC>(new CC());
            }

            FEM& getFEM()
            {
                return femb.getFEM();
            }

            const FEM& getFEM() const
            {
                return femb.getFEM();
            }

            // return gfs reference
            GFS& getGFS ()
            {
                return *gfsp;
            }

            // return gfs reference const version
            const GFS& getGFS () const
            {
                return *gfsp;
            }

            // return gfs reference
            CC& getCC ()
            {
                return *ccp;
            }

            // return gfs reference const version
            const CC& getCC () const
            {
                return *ccp;
            }

            void assembleConstraints (const BCType& bctype)
            {
                ccp->clear();
                constraints(bctype,*gfsp,*ccp);
            }

            void clearConstraints ()
            {
                ccp->clear();
            }

            void setConstrainedDOFS (DOF& x, NT nt) const
            {
                set_constrained_dofs(*ccp,nt,x);
                conb.make_consistent(*gfsp,x);
            }

            void setNonConstrainedDOFS (DOF& x, NT nt) const
            {
                set_nonconstrained_dofs(*ccp,nt,x);
                conb.make_consistent(*gfsp,x);
            }

            void copyConstrainedDOFS (const DOF& xin, DOF& xout) const
            {
                copy_constrained_dofs(*ccp,xin,xout);
                conb.make_consistent(*gfsp,xout);
            }

            void copyNonConstrainedDOFS (const DOF& xin, DOF& xout) const
            {
                copy_nonconstrained_dofs(*ccp,xin,xout);
                conb.make_consistent(*gfsp,xout);
            }

        private:
            GV gv; // need this object here because FEM and GFS store a const reference !!
            ES es;
            FEMB femb;
            CONB conb;
            std::shared_ptr<GFS> gfsp;
            std::shared_ptr<CC> ccp;
        };



        //============================================================================
        // Discontinuous Finite Element Space
        //============================================================================

        // constraints base template
        template<SolverCategory::Category st>
        class DGCONBase
        {};

        template<>
        class DGCONBase<SolverCategory::sequential>
        {
        public:
            typedef NoConstraints CON;
            DGCONBase ()
            {
                conp = std::shared_ptr<CON>(new CON());
            }
            CON& getCON() {return *conp;}
            const CON& getCON() const {return *conp;}
            template<typename GFS, typename DOF>
            void make_consistent (const GFS& gfs, DOF& x) const {}
        private:
            std::shared_ptr<CON> conp;
        };

        template<>
        class DGCONBase<SolverCategory::nonoverlapping>
        {
        public:
            typedef P0ParallelGhostConstraints CON;
            DGCONBase ()
            {
                conp = std::shared_ptr<CON>(new CON());
            }
            CON& getCON() {return *conp;}
            const CON& getCON() const {return *conp;}
            template<typename GFS, typename DOF>
            void make_consistent (const GFS& gfs, DOF& x) const {}
        private:
            std::shared_ptr<CON> conp;
        };

        template<>
        class DGCONBase<SolverCategory::overlapping>
        {
        public:
            typedef P0ParallelConstraints CON;
            DGCONBase ()
            {
                conp = std::shared_ptr<CON>(new CON());
            }
            CON& getCON() {return *conp;}
            const CON& getCON() const {return *conp;}
            template<typename GFS, typename DOF>
            void make_consistent (const GFS& gfs, DOF& x) const
            {
                // make vector consistent; this is needed for all overlapping solvers
                ISTL::ParallelHelper<GFS> helper(gfs);
                helper.maskForeignDOFs(Backend::native(x));
                Dune::PDELab::AddDataHandle<GFS,DOF> adddh(gfs,x);
                if (gfs.gridView().comm().size()>1)
                    gfs.gridView().communicate(adddh,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);
            }
        private:
            std::shared_ptr<CON> conp;
        };

        // Discontinuous space
        // default implementation, use only specializations below
        template<typename T, typename N, unsigned int degree,
                 Dune::GeometryType::BasicType gt, SolverCategory::Category st = SolverCategory::sequential,
                 typename VBET=ISTL::VectorBackend<ISTL::Blocking::fixed,Dune::PB::PkSize<degree,T::dimension>::value> >
        class DGPkSpace
        {
        public:

            // export types
            typedef T Grid;
            typedef typename T::LeafGridView GV;
            typedef typename T::ctype ctype;
            static const int dim = T::dimension;
            static const int dimworld = T::dimensionworld;
            typedef N NT;
#if HAVE_GMP
            typedef OPBLocalFiniteElementMap<ctype,NT,degree,dim,gt,Dune::GMPField<512>,Dune::PB::BasisType::Pk> FEM;
#else
            typedef OPBLocalFiniteElementMap<ctype,NT,degree,dim,gt> FEM;
#endif
            typedef DGCONBase<st> CONB;
            typedef typename CONB::CON CON;
            typedef VBET VBE;
            typedef GridFunctionSpace<GV,FEM,CON,VBE> GFS;
            using DOF = Backend::Vector<GFS,N>;
            typedef Dune::PDELab::DiscreteGridFunction<GFS,DOF> DGF;
            typedef typename GFS::template ConstraintsContainer<N>::Type CC;
            typedef VTKGridFunctionAdapter<DGF> VTKF;

            // constructor making the grid function space an all that is needed
            DGPkSpace (const GV& gridview) : gv(gridview), conb()
            {
                femp = std::shared_ptr<FEM>(new FEM());
                gfsp = std::shared_ptr<GFS>(new GFS(gv,*femp));
                // initialize ordering
                gfsp->update();
                ccp = std::shared_ptr<CC>(new CC());
            }

            FEM& getFEM() { return *femp; }
            const FEM& getFEM() const { return *femp; }

            // return gfs reference
            GFS& getGFS () { return *gfsp; }

            // return gfs reference const version
            const GFS& getGFS () const {return *gfsp;}

            // return gfs reference
            CC& getCC () { return *ccp;}

            // return gfs reference const version
            const CC& getCC () const { return *ccp;}

            template<class BCTYPE>
            void assembleConstraints (const BCTYPE& bctype)
            {
                ccp->clear();
                constraints(bctype,*gfsp,*ccp);
            }

            void clearConstraints ()
            {
                ccp->clear();
            }

            void setConstrainedDOFS (DOF& x, NT nt) const
            {
                set_constrained_dofs(*ccp,nt,x);
                conb.make_consistent(*gfsp,x);
            }

            void setNonConstrainedDOFS (DOF& x, NT nt) const
            {
                set_nonconstrained_dofs(*ccp,nt,x);
                conb.make_consistent(*gfsp,x);
            }

            void copyConstrainedDOFS (const DOF& xin, DOF& xout) const
            {
                copy_constrained_dofs(*ccp,xin,xout);
                conb.make_consistent(*gfsp,xout);
            }

            void copyNonConstrainedDOFS (const DOF& xin, DOF& xout) const
            {
                copy_nonconstrained_dofs(*ccp,xin,xout);
                conb.make_consistent(*gfsp,xout);
            }

        private:
            GV gv; // need this object here because FEM and GFS store a const reference !!
            CONB conb;
            std::shared_ptr<FEM> femp;
            std::shared_ptr<GFS> gfsp;
            std::shared_ptr<CC> ccp;
        };

        // Discontinuous space
        // default implementation, use only specializations below
        template<typename T, typename N, unsigned int degree,
                 Dune::GeometryType::BasicType gt, SolverCategory::Category st = SolverCategory::sequential,
                 //typename VBET=ISTL::VectorBackend<ISTL::Blocking::fixed,Dune::PB::PkSize<degree,T::dimension>::value> >
                 typename VBET=ISTL::VectorBackend<> >
        class DGQkOPBSpace
        {
        public:

            // export types
            typedef T Grid;
            typedef typename T::LeafGridView GV;
            typedef typename T::ctype ctype;
            static const int dim = T::dimension;
            static const int dimworld = T::dimensionworld;
            typedef N NT;
#if HAVE_GMP
            typedef OPBLocalFiniteElementMap<ctype,NT,degree,dim,gt,Dune::GMPField<512>,Dune::PB::BasisType::Qk> FEM;
#else
            typedef OPBLocalFiniteElementMap<ctype,NT,degree,dim,gt,N,Dune::PB::BasisType::Qk> FEM;
#endif
            typedef DGCONBase<st> CONB;
            typedef typename CONB::CON CON;
            typedef VBET VBE;
            typedef GridFunctionSpace<GV,FEM,CON,VBE> GFS;
            using DOF = Backend::Vector<GFS,N>;
            typedef Dune::PDELab::DiscreteGridFunction<GFS,DOF> DGF;
            typedef typename GFS::template ConstraintsContainer<N>::Type CC;
            typedef VTKGridFunctionAdapter<DGF> VTKF;

            // constructor making the grid function space an all that is needed
            DGQkOPBSpace (const GV& gridview) : gv(gridview), conb()
            {
                femp = std::shared_ptr<FEM>(new FEM());
                gfsp = std::shared_ptr<GFS>(new GFS(gv,*femp));
                // initialize ordering
                gfsp->update();
                ccp = std::shared_ptr<CC>(new CC());
            }

            FEM& getFEM() { return *femp; }
            const FEM& getFEM() const { return *femp; }

            // return gfs reference
            GFS& getGFS () { return *gfsp; }

            // return gfs reference const version
            const GFS& getGFS () const {return *gfsp;}

            // return gfs reference
            CC& getCC () { return *ccp;}

            // return gfs reference const version
            const CC& getCC () const { return *ccp;}

            template<class BCTYPE>
            void assembleConstraints (const BCTYPE& bctype)
            {
                ccp->clear();
                constraints(bctype,*gfsp,*ccp);
            }

            void clearConstraints ()
            {
                ccp->clear();
            }

            void setConstrainedDOFS (DOF& x, NT nt) const
            {
                set_constrained_dofs(*ccp,nt,x);
                conb.make_consistent(*gfsp,x);
            }

            void setNonConstrainedDOFS (DOF& x, NT nt) const
            {
                set_nonconstrained_dofs(*ccp,nt,x);
                conb.make_consistent(*gfsp,x);
            }

            void copyConstrainedDOFS (const DOF& xin, DOF& xout) const
            {
                copy_constrained_dofs(*ccp,xin,xout);
                conb.make_consistent(*gfsp,xout);
            }

            void copyNonConstrainedDOFS (const DOF& xin, DOF& xout) const
            {
                copy_nonconstrained_dofs(*ccp,xin,xout);
                conb.make_consistent(*gfsp,xout);
            }

        private:
            GV gv; // need this object here because FEM and GFS store a const reference !!
            CONB conb;
            std::shared_ptr<FEM> femp;
            std::shared_ptr<GFS> gfsp;
            std::shared_ptr<CC> ccp;
        };

        // Discontinuous space
        // default implementation, use only specializations below
        template<typename T, typename N, unsigned int degree,
                 Dune::GeometryType::BasicType gt, SolverCategory::Category st = SolverCategory::sequential,
                 typename VBET=ISTL::VectorBackend<ISTL::Blocking::fixed,Dune::QkStuff::QkSize<degree,T::dimension>::value> >
        class DGQkSpace
        {
        public:

            // export types
            typedef T Grid;
            typedef typename T::LeafGridView GV;
            typedef typename T::ctype ctype;
            static const int dim = T::dimension;
            static const int dimworld = T::dimensionworld;
            typedef N NT;
            typedef QkDGLocalFiniteElementMap<ctype,NT,degree,dim> FEM;
            typedef DGCONBase<st> CONB;
            typedef typename CONB::CON CON;
            typedef VBET VBE;
            typedef GridFunctionSpace<GV,FEM,CON,VBE> GFS;
            using DOF = Backend::Vector<GFS,N>;
            typedef Dune::PDELab::DiscreteGridFunction<GFS,DOF> DGF;
            typedef typename GFS::template ConstraintsContainer<N>::Type CC;
            typedef VTKGridFunctionAdapter<DGF> VTKF;

            // constructor making the grid function space an all that is needed
            DGQkSpace (const GV& gridview) : gv(gridview), conb()
            {
                femp = std::shared_ptr<FEM>(new FEM());
                gfsp = std::shared_ptr<GFS>(new GFS(gv,*femp));
                // initialize ordering
                gfsp->update();
                ccp = std::shared_ptr<CC>(new CC());
            }

            FEM& getFEM() { return *femp; }
            const FEM& getFEM() const { return *femp; }

            // return gfs reference
            GFS& getGFS () { return *gfsp; }

            // return gfs reference const version
            const GFS& getGFS () const {return *gfsp;}

            // return gfs reference
            CC& getCC () { return *ccp;}

            // return gfs reference const version
            const CC& getCC () const { return *ccp;}

            template<class BCTYPE>
            void assembleConstraints (const BCTYPE& bctype)
            {
                ccp->clear();
                constraints(bctype,*gfsp,*ccp);
            }

            void clearConstraints ()
            {
                ccp->clear();
            }

            void setConstrainedDOFS (DOF& x, NT nt) const
            {
                set_constrained_dofs(*ccp,nt,x);
                conb.make_consistent(*gfsp,x);
            }

            void setNonConstrainedDOFS (DOF& x, NT nt) const
            {
                set_nonconstrained_dofs(*ccp,nt,x);
                conb.make_consistent(*gfsp,x);
            }

            void copyConstrainedDOFS (const DOF& xin, DOF& xout) const
            {
                copy_constrained_dofs(*ccp,xin,xout);
                conb.make_consistent(*gfsp,xout);
            }

            void copyNonConstrainedDOFS (const DOF& xin, DOF& xout) const
            {
                copy_nonconstrained_dofs(*ccp,xin,xout);
                conb.make_consistent(*gfsp,xout);
            }

        private:
            GV gv; // need this object here because FEM and GFS store a const reference !!
            CONB conb;
            std::shared_ptr<FEM> femp;
            std::shared_ptr<GFS> gfsp;
            std::shared_ptr<CC> ccp;
        };


        // Discontinuous space using QK with Gauss Lobatto points (use only for cube elements)
        template<typename T, typename N, unsigned int degree,
                 Dune::GeometryType::BasicType gt, SolverCategory::Category st = SolverCategory::sequential,
                 //typename VBET=ISTL::VectorBackend<ISTL::Blocking::fixed,Dune::QkStuff::QkSize<degree,T::dimension>::value> >
                 typename VBET=ISTL::VectorBackend<> >
        class DGQkGLSpace
        {
        public:

            // export types
            typedef T Grid;
            typedef typename T::LeafGridView GV;
            typedef typename T::ctype ctype;
            static const int dim = T::dimension;
            static const int dimworld = T::dimensionworld;
            typedef N NT;
            typedef QkDGLocalFiniteElementMap<ctype,NT,degree,dim,QkDGBasisPolynomial::lobatto> FEM;
            typedef DGCONBase<st> CONB;
            typedef typename CONB::CON CON;
            typedef VBET VBE;
            typedef GridFunctionSpace<GV,FEM,CON,VBE> GFS;
            using DOF = Backend::Vector<GFS,N>;
            typedef Dune::PDELab::DiscreteGridFunction<GFS,DOF> DGF;
            typedef typename GFS::template ConstraintsContainer<N>::Type CC;
            typedef VTKGridFunctionAdapter<DGF> VTKF;

            // constructor making the grid function space an all that is needed
            DGQkGLSpace (const GV& gridview) : gv(gridview), conb()
            {
                femp = std::shared_ptr<FEM>(new FEM());
                gfsp = std::shared_ptr<GFS>(new GFS(gv,*femp));
                // initialize ordering
                gfsp->update();
                ccp = std::shared_ptr<CC>(new CC());
            }

            FEM& getFEM() { return *femp; }
            const FEM& getFEM() const { return *femp; }

            // return gfs reference
            GFS& getGFS () { return *gfsp; }

            // return gfs reference const version
            const GFS& getGFS () const {return *gfsp;}

            // return gfs reference
            CC& getCC () { return *ccp;}

            // return gfs reference const version
            const CC& getCC () const { return *ccp;}

            template<class BCTYPE>
            void assembleConstraints (const BCTYPE& bctype)
            {
                ccp->clear();
                constraints(bctype,*gfsp,*ccp);
            }

            void clearConstraints ()
            {
                ccp->clear();
            }

            void setConstrainedDOFS (DOF& x, NT nt) const
            {
                set_constrained_dofs(*ccp,nt,x);
                conb.make_consistent(*gfsp,x);
            }

            void setNonConstrainedDOFS (DOF& x, NT nt) const
            {
                set_nonconstrained_dofs(*ccp,nt,x);
                conb.make_consistent(*gfsp,x);
            }

            void copyConstrainedDOFS (const DOF& xin, DOF& xout) const
            {
                copy_constrained_dofs(*ccp,xin,xout);
                conb.make_consistent(*gfsp,xout);
            }

            void copyNonConstrainedDOFS (const DOF& xin, DOF& xout) const
            {
                copy_nonconstrained_dofs(*ccp,xin,xout);
                conb.make_consistent(*gfsp,xout);
            }

        private:
            GV gv; // need this object here because FEM and GFS store a const reference !!
            CONB conb;
            std::shared_ptr<FEM> femp;
            std::shared_ptr<GFS> gfsp;
            std::shared_ptr<CC> ccp;
        };


        // Discontinuous space using Legendre polynomials (use only for cube elements)
        template<typename T, typename N, unsigned int degree,
                 Dune::GeometryType::BasicType gt, SolverCategory::Category st = SolverCategory::sequential,
                 //typename VBET=ISTL::VectorBackend<ISTL::Blocking::fixed,Dune::QkStuff::QkSize<degree,T::dimension>::value> >
                 typename VBET=ISTL::VectorBackend<> >
        class DGLegendreSpace
        {
        public:

            // export types
            typedef T Grid;
            typedef typename T::LeafGridView GV;
            typedef typename T::ctype ctype;
            static const int dim = T::dimension;
            static const int dimworld = T::dimensionworld;
            typedef N NT;
            typedef QkDGLocalFiniteElementMap<ctype,NT,degree,dim,QkDGBasisPolynomial::legendre> FEM;
            typedef DGCONBase<st> CONB;
            typedef typename CONB::CON CON;
            typedef VBET VBE;
            typedef GridFunctionSpace<GV,FEM,CON,VBE> GFS;
            using DOF = Backend::Vector<GFS,N>;
            typedef Dune::PDELab::DiscreteGridFunction<GFS,DOF> DGF;
            typedef typename GFS::template ConstraintsContainer<N>::Type CC;
            typedef VTKGridFunctionAdapter<DGF> VTKF;

            // constructor making the grid function space an all that is needed
            DGLegendreSpace (const GV& gridview) : gv(gridview), conb()
            {
                femp = std::shared_ptr<FEM>(new FEM());
                gfsp = std::shared_ptr<GFS>(new GFS(gv,*femp));
                // initialize ordering
                gfsp->update();
                ccp = std::shared_ptr<CC>(new CC());
            }

            FEM& getFEM() { return *femp; }
            const FEM& getFEM() const { return *femp; }

            // return gfs reference
            GFS& getGFS () { return *gfsp; }

            // return gfs reference const version
            const GFS& getGFS () const {return *gfsp;}

            // return gfs reference
            CC& getCC () { return *ccp;}

            // return gfs reference const version
            const CC& getCC () const { return *ccp;}

            template<class BCTYPE>
            void assembleConstraints (const BCTYPE& bctype)
            {
                ccp->clear();
                constraints(bctype,*gfsp,*ccp);
            }

            void clearConstraints ()
            {
                ccp->clear();
            }

            void setConstrainedDOFS (DOF& x, NT nt) const
            {
                set_constrained_dofs(*ccp,nt,x);
                conb.make_consistent(*gfsp,x);
            }

            void setNonConstrainedDOFS (DOF& x, NT nt) const
            {
                set_nonconstrained_dofs(*ccp,nt,x);
                conb.make_consistent(*gfsp,x);
            }

            void copyConstrainedDOFS (const DOF& xin, DOF& xout) const
            {
                copy_constrained_dofs(*ccp,xin,xout);
                conb.make_consistent(*gfsp,xout);
            }

            void copyNonConstrainedDOFS (const DOF& xin, DOF& xout) const
            {
                copy_nonconstrained_dofs(*ccp,xin,xout);
                conb.make_consistent(*gfsp,xout);
            }

        private:
            GV gv; // need this object here because FEM and GFS store a const reference !!
            CONB conb;
            std::shared_ptr<FEM> femp;
            std::shared_ptr<GFS> gfsp;
            std::shared_ptr<CC> ccp;
        };


        // Discontinuous P0 space
        template<typename T, typename N,
                 Dune::GeometryType::BasicType gt, SolverCategory::Category st = SolverCategory::sequential,
                 typename VBET=ISTL::VectorBackend<> >
        class P0Space
        {
        public:

            // export types
            typedef T Grid;
            typedef typename T::LeafGridView GV;
            typedef typename T::ctype ctype;
            static const int dim = T::dimension;
            static const int dimworld = T::dimensionworld;
            typedef N NT;
            typedef Dune::PDELab::P0LocalFiniteElementMap<ctype,NT,dim> FEM;
            typedef DGCONBase<st> CONB;
            typedef typename CONB::CON CON;
            typedef VBET VBE;
            typedef GridFunctionSpace<GV,FEM,CON,VBE> GFS;
            using DOF = Backend::Vector<GFS,N>;
            typedef Dune::PDELab::DiscreteGridFunction<GFS,DOF> DGF;
            typedef typename GFS::template ConstraintsContainer<N>::Type CC;
            typedef VTKGridFunctionAdapter<DGF> VTKF;

            // constructor making the grid function space an all that is needed
            P0Space (const GV& gridview) : gv(gridview), conb()
            {
                femp = std::shared_ptr<FEM>(new FEM(Dune::GeometryType(gt,dim)));
                gfsp = std::shared_ptr<GFS>(new GFS(gv,*femp));
                // initialize ordering
                gfsp->update();
                ccp = std::shared_ptr<CC>(new CC());
            }

            FEM& getFEM() { return *femp; }
            const FEM& getFEM() const { return *femp; }

            // return gfs reference
            GFS& getGFS () { return *gfsp; }

            // return gfs reference const version
            const GFS& getGFS () const {return *gfsp;}

            // return gfs reference
            CC& getCC () { return *ccp;}

            // return gfs reference const version
            const CC& getCC () const { return *ccp;}

            template<class BCTYPE>
            void assembleConstraints (const BCTYPE& bctype)
            {
                ccp->clear();
                constraints(bctype,*gfsp,*ccp);
            }

            void clearConstraints ()
            {
                ccp->clear();
            }

            void setConstrainedDOFS (DOF& x, NT nt) const
            {
                set_constrained_dofs(*ccp,nt,x);
                conb.make_consistent(*gfsp,x);
            }

            void setNonConstrainedDOFS (DOF& x, NT nt) const
            {
                set_nonconstrained_dofs(*ccp,nt,x);
                conb.make_consistent(*gfsp,x);
            }

            void copyConstrainedDOFS (const DOF& xin, DOF& xout) const
            {
                copy_constrained_dofs(*ccp,xin,xout);
                conb.make_consistent(*gfsp,xout);
            }

            void copyNonConstrainedDOFS (const DOF& xin, DOF& xout) const
            {
                copy_nonconstrained_dofs(*ccp,xin,xout);
                conb.make_consistent(*gfsp,xout);
            }

        private:
            GV gv; // need this object here because FEM and GFS store a const reference !!
            CONB conb;
            std::shared_ptr<FEM> femp;
            std::shared_ptr<GFS> gfsp;
            std::shared_ptr<CC> ccp;
        };


        // how can we most easily specify a grid function
        // pass a function space as parameter
        template<typename FS, typename Functor>
        class UserFunction
            : public GridFunctionBase<GridFunctionTraits<typename FS::GV, typename FS::NT,
                                                         1,FieldVector<typename FS::NT,1> >
                                      ,UserFunction<FS,Functor> >
        {
        public:
            typedef GridFunctionTraits<typename FS::GV, typename FS::NT,
                                       1,FieldVector<typename FS::NT,1> > Traits;

            //! constructor
            UserFunction (const FS& fs_, const Functor& f_)
                : fs(fs_), f(f_)
            {}

            //! \copydoc GridFunctionBase::evaluate()
            inline void evaluate (const typename Traits::ElementType& e,
                                  const typename Traits::DomainType& x,
                                  typename Traits::RangeType& y) const
            {
                typename Traits::DomainType x_ = e.geometry().global(x);
                std::vector<double> x__(x.size());
                for (size_t i=0; i<x.size(); ++i) x__[i]=x_[i];
                y = f(x__);
            }

            inline const typename FS::GV& getGridView () const
            {
                return fs.getGFS().gridView();
            }

        private:
            const FS fs; // store a copy of the function space
            const Functor f;
        };


        template<typename FS, typename LOP, SolverCategory::Category st = SolverCategory::sequential>
        class GalerkinGlobalAssembler
        {
        public:
            // export types
            typedef ISTL::BCRSMatrixBackend<> MBE;
            typedef Dune::PDELab::GridOperator<typename FS::GFS,typename FS::GFS,LOP,MBE,
                                               typename FS::NT,typename FS::NT,typename FS::NT,
                                               typename FS::CC,typename FS::CC> GO;
            typedef typename GO::Jacobian MAT;

            GalerkinGlobalAssembler (const FS& fs, LOP& lop, const std::size_t nonzeros)
            {
                gop = std::shared_ptr<GO>(new GO(fs.getGFS(),fs.getCC(),fs.getGFS(),fs.getCC(),lop,MBE(nonzeros)));
            }

            // return grid reference
            GO& getGO ()
            {
                return *gop;
            }

            // return grid reference const version
            const GO& getGO () const
            {
                return *gop;
            }

            GO& operator*()
            {
                return *gop;
            }

            GO* operator->()
            {
                return gop.operator->();
            }

            const GO& operator*() const
            {
                return *gop;
            }

            const GO* operator->() const
            {
                return gop.operator->();
            }

        private:
            std::shared_ptr<GO> gop;
        };


        template<typename FS, typename LOP, SolverCategory::Category st = SolverCategory::sequential>
        class GalerkinGlobalAssemblerNewBackend
        {
        public:
            // export types
            typedef Dune::PDELab::ISTL::BCRSMatrixBackend<> MBE;
            typedef Dune::PDELab::GridOperator<typename FS::GFS,typename FS::GFS,LOP,MBE,
                                               typename FS::NT,typename FS::NT,typename FS::NT,
                                               typename FS::CC,typename FS::CC> GO;
            typedef typename GO::Jacobian MAT;

            GalerkinGlobalAssemblerNewBackend (const FS& fs, LOP& lop, const MBE& mbe)
            {
                gop = std::shared_ptr<GO>(new GO(fs.getGFS(),fs.getCC(),fs.getGFS(),fs.getCC(),lop,mbe));
            }

            // return grid reference
            GO& getGO ()
            {
                return *gop;
            }

            // return grid reference const version
            const GO& getGO () const
            {
                return *gop;
            }

            GO& operator*()
            {
                return *gop;
            }

            GO* operator->()
            {
                return gop.operator->();
            }

            const GO& operator*() const
            {
                return *gop;
            }

            const GO* operator->() const
            {
                return gop.operator->();
            }

        private:
            std::shared_ptr<GO> gop;
        };


        // variant with two different function spaces
        template<typename FSU, typename FSV, typename LOP, SolverCategory::Category st>
        class GlobalAssembler
        {
        public:
            // export types
            typedef ISTL::BCRSMatrixBackend<> MBE;
            typedef Dune::PDELab::GridOperator<typename FSU::GFS,typename FSV::GFS,LOP,MBE,
                                               typename FSU::NT,typename FSU::NT,typename FSU::NT,
                                               typename FSU::CC,typename FSV::CC> GO;
            typedef typename GO::Jacobian MAT;

            GlobalAssembler (const FSU& fsu, const FSV& fsv, LOP& lop, const std::size_t nonzeros)
            {
                gop = std::shared_ptr<GO>(new GO(fsu.getGFS(),fsu.getCC(),fsv.getGFS(),fsv.getCC(),lop,MBE(nonzeros)));
            }

            // return grid reference
            GO& getGO ()
            {
                return *gop;
            }

            // return grid reference const version
            const GO& getGO () const
            {
                return *gop;
            }

            GO& operator*()
            {
                return *gop;
            }

            GO* operator->()
            {
                return gop.operator->();
            }

            const GO& operator*() const
            {
                return *gop;
            }

            const GO* operator->() const
            {
                return gop.operator->();
            }

        private:
            std::shared_ptr<GO> gop;
        };


        template<typename GO1, typename GO2, bool implicit = true>
        class OneStepGlobalAssembler
        {
        public:
            // export types
            typedef ISTL::BCRSMatrixBackend<> MBE;
            typedef Dune::PDELab::OneStepGridOperator<typename GO1::GO,typename GO2::GO,implicit> GO;
            typedef typename GO::Jacobian MAT;

            OneStepGlobalAssembler (GO1& go1, GO2& go2)
            {
                gop = std::shared_ptr<GO>(new GO(*go1,*go2));
            }

            // return grid reference
            GO& getGO ()
            {
                return *gop;
            }

            // return grid reference const version
            const GO& getGO () const
            {
                return *gop;
            }

            GO& operator*()
            {
                return *gop;
            }

            GO* operator->()
            {
                return gop.operator->();
            }

            const GO& operator*() const
            {
                return *gop;
            }

            const GO* operator->() const
            {
                return gop.operator->();
            }

        private:
            std::shared_ptr<GO> gop;
        };


        // packaging of the CG_AMG_SSOR solver: default version is sequential
        template<typename FS, typename ASS, SolverCategory::Category st = SolverCategory::sequential>
        class ISTLSolverBackend_CG_AMG_SSOR
        {
        public:
            // types exported
            typedef Dune::PDELab::ISTLBackend_SEQ_CG_AMG_SSOR<typename ASS::GO> LS;

            ISTLSolverBackend_CG_AMG_SSOR (const FS& fs, const ASS& ass, unsigned maxiter_=5000,
                                           int verbose_=1, bool reuse_=false, bool usesuperlu_=true)
            {
                lsp = std::shared_ptr<LS>(new LS(maxiter_,verbose_,reuse_,usesuperlu_));
            }

            LS& getLS () {return *lsp;}
            const LS& getLS () const { return *lsp;}
            LS& operator*(){return *lsp;}
            LS* operator->() { return lsp.operator->(); }
            const LS& operator*() const{return *lsp;}
            const LS* operator->() const{ return lsp.operator->();}

       private:
            std::shared_ptr<LS> lsp;
        };

        // packaging of the CG_AMG_SSOR solver: nonoverlapping version
        template<typename FS, typename ASS>
        class ISTLSolverBackend_CG_AMG_SSOR<FS,ASS, SolverCategory::nonoverlapping>
        {
        public:
            // types exported
            typedef Dune::PDELab::ISTLBackend_NOVLP_CG_AMG_SSOR<typename ASS::GO> LS;

            ISTLSolverBackend_CG_AMG_SSOR (const FS& fs, const ASS& ass, unsigned maxiter_=5000,
                                           int verbose_=1, bool reuse_=false, bool usesuperlu_=true)
            {
                lsp = std::shared_ptr<LS>(new LS(fs.getGFS(),maxiter_,verbose_,reuse_,usesuperlu_));
            }

            LS& getLS () {return *lsp;}
            const LS& getLS () const { return *lsp;}
            LS& operator*(){return *lsp;}
            LS* operator->() { return lsp.operator->(); }
            const LS& operator*() const{return *lsp;}
            const LS* operator->() const{ return lsp.operator->();}

       private:
            std::shared_ptr<LS> lsp;
        };

        // packaging of the CG_AMG_SSOR solver: overlapping version
        template<typename FS, typename ASS>
        class ISTLSolverBackend_CG_AMG_SSOR<FS,ASS, SolverCategory::overlapping>
        {
        public:
            // types exported
            typedef Dune::PDELab::ISTLBackend_CG_AMG_SSOR<typename ASS::GO> LS;

            ISTLSolverBackend_CG_AMG_SSOR (const FS& fs, const ASS& ass, unsigned maxiter_=5000,
                                           int verbose_=1, bool reuse_=false, bool usesuperlu_=true)
            {
                lsp = std::shared_ptr<LS>(new LS(fs.getGFS(),maxiter_,verbose_,reuse_,usesuperlu_));
            }

            LS& getLS () {return *lsp;}
            const LS& getLS () const { return *lsp;}
            LS& operator*(){return *lsp;}
            LS* operator->() { return lsp.operator->(); }
            const LS& operator*() const{return *lsp;}
            const LS* operator->() const{ return lsp.operator->();}

        private:
            std::shared_ptr<LS> lsp;
        };

        // packaging of the CG_SSOR solver: default version is sequential
        template<typename FS, typename ASS, SolverCategory::Category st = SolverCategory::sequential>
        class ISTLSolverBackend_CG_SSOR
        {
        public:
            // types exported
            typedef ISTLBackend_SEQ_CG_SSOR LS;

            ISTLSolverBackend_CG_SSOR (const FS& fs, const ASS& ass, unsigned maxiter_=5000,
                                       int steps_=5, int verbose_=1)
            {
                lsp = std::shared_ptr<LS>(new LS(maxiter_,verbose_));
            }

            LS& getLS () {return *lsp;}
            const LS& getLS () const { return *lsp;}
            LS& operator*(){return *lsp;}
            LS* operator->() { return lsp.operator->(); }
            const LS& operator*() const{return *lsp;}
            const LS* operator->() const{ return lsp.operator->();}

       private:
            std::shared_ptr<LS> lsp;
        };

        // packaging of the CG_SSOR solver: nonoverlapping version
        template<typename FS, typename ASS>
        class ISTLSolverBackend_CG_SSOR<FS,ASS,SolverCategory::nonoverlapping>
        {
        public:
            // types exported
            typedef ISTLBackend_NOVLP_CG_SSORk<typename ASS::GO> LS;

            ISTLSolverBackend_CG_SSOR (const FS& fs, const ASS& ass, unsigned maxiter_=5000,
                                       int steps_=5, int verbose_=1)
            {
                lsp = std::shared_ptr<LS>(new LS(fs.getGFS(),maxiter_,steps_,verbose_));
            }

            LS& getLS () {return *lsp;}
            const LS& getLS () const { return *lsp;}
            LS& operator*(){return *lsp;}
            LS* operator->() { return lsp.operator->(); }
            const LS& operator*() const{return *lsp;}
            const LS* operator->() const{ return lsp.operator->();}

       private:
            std::shared_ptr<LS> lsp;
        };

        // packaging of the CG_SSOR solver: overlapping version
        template<typename FS, typename ASS>
        class ISTLSolverBackend_CG_SSOR<FS,ASS,SolverCategory::overlapping>
        {
        public:
            // types exported
            typedef ISTLBackend_OVLP_CG_SSORk<typename FS::GFS, typename FS::CC> LS;

            ISTLSolverBackend_CG_SSOR (const FS& fs, const ASS& ass, unsigned maxiter_=5000,
                                       int steps_=5, int verbose_=1)
            {
                lsp = std::shared_ptr<LS>(new LS(fs.getGFS(),fs.getCC(),maxiter_,steps_,verbose_));
            }

            LS& getLS () {return *lsp;}
            const LS& getLS () const { return *lsp;}
            LS& operator*(){return *lsp;}
            LS* operator->() { return lsp.operator->(); }
            const LS& operator*() const{return *lsp;}
            const LS* operator->() const{ return lsp.operator->();}

       private:
            std::shared_ptr<LS> lsp;
        };


        // packaging of a default solver that should always work
        // in the sequential case : BCGS SSOR
        template<typename FS, typename ASS, SolverCategory::Category st = SolverCategory::sequential>
        class ISTLSolverBackend_IterativeDefault
        {
        public:
            // types exported
            typedef ISTLBackend_SEQ_BCGS_SSOR LS;

            ISTLSolverBackend_IterativeDefault (const FS& fs, const ASS& ass, unsigned maxiter_=5000, int verbose_=1)
            {
                lsp = std::shared_ptr<LS>(new LS(maxiter_,verbose_));
            }

            LS& getLS () {return *lsp;}
            const LS& getLS () const { return *lsp;}
            LS& operator*(){return *lsp;}
            LS* operator->() { return lsp.operator->(); }
            const LS& operator*() const{return *lsp;}
            const LS* operator->() const{ return lsp.operator->();}

       private:
            std::shared_ptr<LS> lsp;
        };

        // in the nonoverlapping case : BCGS SSORk
        template<typename FS, typename ASS>
        class ISTLSolverBackend_IterativeDefault<FS,ASS,SolverCategory::nonoverlapping>
        {
        public:
            // types exported
            typedef Dune::PDELab::ISTLBackend_NOVLP_BCGS_SSORk<typename ASS::GO> LS;

            ISTLSolverBackend_IterativeDefault (const FS& fs, const ASS& ass, unsigned maxiter_=5000, int verbose_=1)
            {
                lsp = std::shared_ptr<LS>(new LS(ass.getGO(),maxiter_,3,verbose_));
            }

            LS& getLS () {return *lsp;}
            const LS& getLS () const { return *lsp;}
            LS& operator*(){return *lsp;}
            LS* operator->() { return lsp.operator->(); }
            const LS& operator*() const{return *lsp;}
            const LS* operator->() const{ return lsp.operator->();}

       private:
            std::shared_ptr<LS> lsp;
        };

        // in the overlapping case : BCGS SSORk
        template<typename FS, typename ASS>
        class ISTLSolverBackend_IterativeDefault<FS,ASS,SolverCategory::overlapping>
        {
        public:
            // types exported
            typedef ISTLBackend_OVLP_BCGS_SSORk<typename FS::GFS, typename FS::CC> LS;

            ISTLSolverBackend_IterativeDefault (const FS& fs, const ASS& ass, unsigned maxiter_=5000, int verbose_=1)
            {
                lsp = std::shared_ptr<LS>(new LS(fs.getGFS(),fs.getCC(),maxiter_,3,verbose_));
            }

            LS& getLS () {return *lsp;}
            const LS& getLS () const { return *lsp;}
            LS& operator*(){return *lsp;}
            LS* operator->() { return lsp.operator->(); }
            const LS& operator*() const{return *lsp;}
            const LS* operator->() const{ return lsp.operator->();}

       private:
            std::shared_ptr<LS> lsp;
        };

        // packaging of a default solver that should always work
        // in the sequential case : BCGS SSOR
        template<typename FS, typename ASS, SolverCategory::Category st = SolverCategory::sequential>
        class ISTLSolverBackend_ExplicitDiagonal
        {
        public:
            // types exported
            typedef Dune::PDELab::ISTLBackend_SEQ_ExplicitDiagonal LS;

            ISTLSolverBackend_ExplicitDiagonal (const FS& fs, const ASS& ass, unsigned maxiter_=5000, int verbose_=1)
            {
                lsp = std::shared_ptr<LS>(new LS());
            }

            LS& getLS () {return *lsp;}
            const LS& getLS () const { return *lsp;}
            LS& operator*(){return *lsp;}
            LS* operator->() { return lsp.operator->(); }
            const LS& operator*() const{return *lsp;}
            const LS* operator->() const{ return lsp.operator->();}

       private:
            std::shared_ptr<LS> lsp;
        };

        // packaging of a default solver that should always work
        // in the sequential case : BCGS SSOR
        template<typename FS, typename ASS>
        class ISTLSolverBackend_ExplicitDiagonal<FS,ASS,SolverCategory::overlapping>
        {
        public:
            // types exported
            typedef Dune::PDELab::ISTLBackend_OVLP_ExplicitDiagonal<typename FS::GFS> LS;

            ISTLSolverBackend_ExplicitDiagonal (const FS& fs, const ASS& ass, unsigned maxiter_=5000, int verbose_=1)
            {
                lsp = std::shared_ptr<LS>(new LS(fs.getGFS()));
            }

            LS& getLS () {return *lsp;}
            const LS& getLS () const { return *lsp;}
            LS& operator*(){return *lsp;}
            LS* operator->() { return lsp.operator->(); }
            const LS& operator*() const{return *lsp;}
            const LS* operator->() const{ return lsp.operator->();}

       private:
            std::shared_ptr<LS> lsp;
        };

        // packaging of a default solver that should always work
        // in the sequential case : BCGS SSOR
        template<typename FS, typename ASS>
        class ISTLSolverBackend_ExplicitDiagonal<FS,ASS,SolverCategory::nonoverlapping>
        {
        public:
            // types exported
            typedef Dune::PDELab::ISTLBackend_NOVLP_ExplicitDiagonal<typename FS::GFS> LS;

            ISTLSolverBackend_ExplicitDiagonal (const FS& fs, const ASS& ass, unsigned maxiter_=5000, int verbose_=1)
            {
                lsp = std::shared_ptr<LS>(new LS(fs.getGFS()));
            }

            LS& getLS () {return *lsp;}
            const LS& getLS () const { return *lsp;}
            LS& operator*(){return *lsp;}
            LS* operator->() { return lsp.operator->(); }
            const LS& operator*() const{return *lsp;}
            const LS* operator->() const{ return lsp.operator->();}

       private:
            std::shared_ptr<LS> lsp;
        };


} // end namespace PDELab
    } // end namespace Dune

#endif // DUNE_PDELAB_BOILERPLATE_PDELAB_HH
