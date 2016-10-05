#include "config.h"

#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/finiteelementmap/pkfem.hh>
#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/localoperator/l2.hh>

#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/gridfunctionspace/vtk.hh>

#include <dune/istl/io.hh>




class AssigningLOP
  : public Dune::PDELab::FullVolumePattern
  , public Dune::PDELab::LocalOperatorDefaultFlags
{

public:

  // pattern assembly flags
  enum { doPatternVolume = true };

  // residual assembly flags
  enum { doAlphaVolume = true };

  // jacobian of volume term
  template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
  void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                        M & mat) const
  {

    for (std::size_t i = 0; i < lfsv.size(); ++i)
      for (std::size_t j = 0; j < lfsu.size(); ++j)
        mat.accumulate(lfsv,i,lfsu,j,1.0);
  }

};


template<typename GV, typename FEM>
void run_test(const GV& gv, const FEM& fem, std::string suffix)
{

  using EntitySet = Dune::PDELab::NonOverlappingEntitySet<GV>;
  auto entity_set = EntitySet(gv);

  using RF = double;

  using VBE = Dune::PDELab::ISTL::VectorBackend<>;
  using NoConstraints = Dune::PDELab::NoConstraints;

  using GFS = Dune::PDELab::GridFunctionSpace<EntitySet,FEM,NoConstraints,VBE>;
  GFS gfs(entity_set,fem);
  gfs.name("x");

  using LOP = AssigningLOP;
  LOP lop;

  using CC = Dune::PDELab::EmptyTransformation;

  using MBE = Dune::PDELab::ISTL::BCRSMatrixBackend<>;

  using GridOperator = Dune::PDELab::GridOperator<
    GFS,GFS,
    LOP,
    MBE,
    RF,RF,RF,CC,CC
    >;

  GridOperator grid_operator(gfs,gfs,lop,MBE(9));

  typename GridOperator::Traits::Range x(gfs,0.0);

  typename GridOperator::Traits::Jacobian A(grid_operator,0.0);

  grid_operator.jacobian(x,A);
  grid_operator.make_consistent(A);

  auto& comm = gv.grid().comm();

  std::cout << comm.rank() << ": Starting VTK output" << std::endl;

  Dune::VTKWriter<GV> vtk_writer(gv);
  vtk_writer.write("nononverlapping_" + suffix);

  for (int r = 0; r < comm.size(); ++r)
    {
      if (r == comm.rank())
        {
          std::cout << "rank: " << r << std::endl;
          for (const auto& vertex : vertices(entity_set))
            {
              std::cout << std::setw(3) << entity_set.indexSet().index(vertex) << " "
                        << "(" << vertex.geometry().center() << ") "
                        << gv.grid().globalIdSet().id(vertex)
                        << std::endl;
            }
          std::cout << std::endl;

          Dune::printmatrix(std::cout,Dune::PDELab::Backend::native(A),"","");

        }

      comm.barrier();
    }
}


int main(int argc, char** argv)
{
  try {
    Dune::MPIHelper& mpi_helper = Dune::MPIHelper::instance(argc,argv);

    {

      using Grid = Dune::UGGrid<2>;
      auto grid_ptr = Dune::StructuredGridFactory<Grid>::createSimplexGrid({{0.0,0.0}},{{1.0,1.0}},{{4,4}});
      Grid& grid = *grid_ptr;

      grid.loadBalance();

      using GV = Grid::LeafGridView;
      auto gv = grid.leafGridView();

      using DF = Grid::ctype;
      using RF = double;

      using FEM = Dune::PDELab::PkLocalFiniteElementMap<GV,DF,RF,1>;
      FEM fem(gv);

      run_test(gv,fem,"structured_simplex");

    }

    {

      using Grid = Dune::UGGrid<2>;

      Dune::GridFactory<Grid> factory;

      // The following code builds a 2D grid with 5x3 cubes. On 2 MPI ranks, UG will by default
      // partition that grid in a way that introduces a concave corner on each rank, allowing us
      // to test the pattern extender algorithm in the GridOperator.
      if (mpi_helper.getCollectiveCommunication().rank() == 0)
        {

          const std::size_t nx = 5;
          const std::size_t ny = 3;

          std::size_t idx = 0;

          for (std::size_t j = 0; j < ny + 1; ++j)
            for (std::size_t i = 0; i < nx + 1; ++i, ++idx)
              {
                factory.insertVertex({static_cast<double>(i),static_cast<double>(j)});
                std::cout << idx << ": (" << i << "," << j << ")" << std::endl;
              }

          std::vector<unsigned int> corners(4);

          for (std::size_t j = 0; j < ny; ++j)
            for (std::size_t i = 0; i < nx; ++i)
              {
                corners[0] = (nx+1)*j       + i;
                corners[1] = (nx+1)*j       + i + 1;
                corners[2] = (nx+1)*(j + 1) + i;
                corners[3] = (nx+1)*(j + 1) + i + 1;
                factory.insertElement(Dune::GeometryType(Dune::GeometryType::cube,2),corners);
              }
        }

      // need to explicitly wrap in a shared_ptr here
      auto grid_ptr = std::shared_ptr<Grid>(factory.createGrid());

      Grid& grid = *grid_ptr;

      grid.loadBalance();

      using GV = Grid::LeafGridView;
      auto gv = grid.leafGridView();

      using DF = Grid::ctype;
      using RF = double;

      using FEM = Dune::PDELab::QkLocalFiniteElementMap<GV,DF,RF,1>;
      FEM fem(gv);

      run_test(gv,fem,"structured_cube_concave_corner");

    }

    return 0;
  }
  catch (Dune::Exception& e) {
    std::cerr << "Dune exception: " << e << std::endl;
  }
  catch (std::exception& e) {
    std::cerr << "std::exception: " << e.what() << std::endl;
  }
  catch (...) {
    std::cerr << "unknown exception" << std::endl;
  }
  return 1;
}
