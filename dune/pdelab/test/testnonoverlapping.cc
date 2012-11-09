#include "config.h"

#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/dgfparser/gridptr.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>

#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/finiteelementmap/q1fem.hh>
#include <dune/pdelab/constraints/constraints.hh>
#include <dune/pdelab/backend/istlvectorbackend.hh>
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

int main(int argc, char** argv)
{
  try {
    Dune::MPIHelper& mpi_helper = Dune::MPIHelper::instance(argc,argv);

    typedef Dune::UGGrid<2> Grid;

    Dune::GridPtr<Grid> grid_ptr("nonoverlapping1.dgf");
    Grid& grid = *grid_ptr;

    grid.loadBalance();

    typedef Grid::Partition<Dune::InteriorBorder_Partition>::LeafGridView GV;
    GV gv = grid.leafView<Dune::InteriorBorder_Partition>();

    typedef GV::ctype DF;
    typedef double RF;

    typedef Dune::PDELab::Q1LocalFiniteElementMap<DF,RF,2> FEM;
    FEM fem;

    typedef Dune::PDELab::ISTLVectorBackend<> VBE;
    typedef Dune::PDELab::NoConstraints NoConstraints;

    typedef Dune::PDELab::GridFunctionSpace<GV,FEM,NoConstraints,VBE,Dune::PDELab::NonOverlappingLeafOrderingTag> GFS;
    GFS gfs(gv,fem);
    gfs.name("x");

    typedef AssigningLOP LOP;
    LOP lop;

    typedef Dune::PDELab::EmptyTransformation CC;

    typedef Dune::PDELab::GridOperator<
      GFS,GFS,
      LOP,Dune::PDELab::ISTLMatrixBackend,
      RF,RF,RF,CC,CC,true
      > GridOperator;

    GridOperator grid_operator(gfs,gfs,lop);

    GridOperator::Traits::Range x(gfs,0.0);

    GridOperator::Traits::Jacobian A(grid_operator,0.0);

    grid_operator.jacobian(x,A);
    grid_operator.make_consistent(A);

    Dune::VTKWriter<GV> vtk_writer(gv);
    vtk_writer.write("nononverlapping1");

    for (int r = 0; r < mpi_helper.getCollectiveCommunication().size(); ++r)
      {
        if (r == mpi_helper.getCollectiveCommunication().rank())
          {
            std::cout << "rank: " << r << std::endl;
            for (GV::Codim<2>::Iterator it = gv.begin<2>(); it != gv.end<2>(); ++it)
              {
                std::cout << std::setw(3) << gv.indexSet().index(*it) << " "
                          << "(" << it->geometry().center() << ") "
                          << gv.grid().globalIdSet().id(*it)
                          << std::endl;
              }
            std::cout << std::endl;

            Dune::printmatrix(std::cout,A.base(),"","");

          }

        mpi_helper.getCollectiveCommunication().barrier();
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
