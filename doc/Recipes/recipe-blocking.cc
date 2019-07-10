#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>

#include <dune/common/filledarray.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/pdelab.hh>


/**
 * \page recipe-blocking Setting up blocked data structures
 *
 * When solving vector-valued PDEs the result at each grid point must be a vector.
 * When doing this the ordering of the global indices can be chosen, e.g.
 * lexicographic or entity blocked. For many solvers, e.g. AMG, the ordering makes a
 * difference in the efficiency of the solver. Further options, such as choosing to aggregate
 * over the blocks and specifying the norm to be used can produce faster results.
 *
 * First, we have to define one or more scalar grid function spaces to combine
 * \snippet recipe-blocking.cc Scalar grid function space
 *
 * There are two ways to define a vector-valued grid function space. If all scalar grid function
 * spaces are the same a power grid function space can be used. In contrast a composite grid function
 * space allows for different scalar grid function spaces. Both can be used with different blocking.
 *
 * Let's first define a power grid function space with lexiographic ordering
 * \snippet recipe-blocking.cc Lexiographic blocked type
 *
 * And then a composite grid function space with entity blocked ordering
 * \snippet recipe-blocking.cc Entity blocked type
 *
 * Full example code: @ref recipe-blocking.cc
 * \example recipe-blocking.cc
 * See explanation at @ref recipe-blocking
 */

int main(int argc, char** argv)
{
    try{
        // Initialize Mpi
        Dune::MPIHelper::instance(argc, argv);

        // need a grid in order to test grid functions
        constexpr unsigned int dim = 2;
        Dune::FieldVector<double,dim> L(5.0);
        std::array<int,dim> N(Dune::filledArray<dim,int>(64));

        typedef Dune::YaspGrid<dim> Grid;
        Grid grid(L,N);

        using DF = double;
        using GV = decltype(grid.leafGridView());

        // instantiate finite element maps
        typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,double,1> FEM;
        FEM fem(grid.leafGridView());

        //Set up constraints
        typedef Dune::PDELab::ConformingDirichletConstraints Constraints;

        // [Scalar grid function space]
        // Set up scalar grid function space
        typedef Dune::PDELab::GridFunctionSpace<GV, FEM, Constraints, Dune::PDELab::ISTL::VectorBackend<>> SCALAR_GFS;
        SCALAR_GFS U1(grid.leafGridView(),fem); U1.name("U1");
        SCALAR_GFS U2(grid.leafGridView(),fem); U2.name("U2");
        //! [Scalar grid function space]

        {
        // [Lexiographic blocked type]
        // Use lexiographical blocked ordering
        typedef Dune::PDELab::LexicographicOrderingTag LexiographicOrderingTag;

        // Set up power grid function space
        using VBE = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed>; // blocking vector backend
        typedef Dune::PDELab::PowerGridFunctionSpace<SCALAR_GFS,
                                                     dim, // block size
                                                     VBE, // blocked vector backend
                                                     LexiographicOrderingTag> GFS;
        GFS gfs(U1,U2);
        //! [Lexiographic blocked type]
        }

        {
        // [Entity blocked type]
        // Use entity blocked ordering
        typedef Dune::PDELab::EntityBlockedOrderingTag EntityOrderingTag;

        // Setting up a composite grid function space with the same scalar grid function space in both components
        using VBE = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed>; // blocking vector backend
        typedef Dune::PDELab::CompositeGridFunctionSpace<VBE, // blocked vector backend
                                                         EntityOrderingTag,
                                                         SCALAR_GFS, SCALAR_GFS> GFS;
        GFS gfs(U1,U2);
        //! [Entity blocked type]
        }

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
