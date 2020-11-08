#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include <cmath>
#include <algorithm>
#include <memory>
#include <bitset>
#include <array>
#include <string>
#include <sstream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/parametertreeparser.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include <dune/pdelab.hh>

#include <dune/pdelab/test/geneo_unittests/utility.hh>


// symbolic constants for the dimension of the domain
constexpr int TWODIM { 2 };
constexpr int THREEDIM { 3 };


std::string generate_filename(const Dune::ParameterTree& ptree)
{
  const int dim { ptree.get<int>("grid.dim") };
  const std::string puType { ptree.get("partition.type", "") };

  if (!(puType == "standard" || puType == "sarkis"))
    DUNE_THROW(Dune::Exception, "Unknown partition of unity type.");

  std::stringstream ss;

  ss << "pou_" << puType << "_" << dim << "d_"
     << ptree.get<int>("grid.overlap") << "o_"
     << ptree.get<int>("partition.baseDegree") << "b";

  return ss.str();
}


// sanity check: do the partition of unity functions add up to 1 at every
// point? Or more precisely:
// \sum\limits_{j=1}^N R_j^T \Theta_j(v|_{\Omega_j}) = v, \forall v \in V_h
// ? For this, it suffices to verify above relation for the basis functions
template<class GFS, class V>
void perform_sanity_check(const GFS& gfs, V& partUnity)
{
  Dune::PDELab::AddDataHandle<GFS, V> partUnityHandle(gfs, partUnity);
  gfs.gridView().communicate(partUnityHandle, Dune::All_All_Interface,
                             Dune::ForwardCommunication);

  // in theory we could also write the partition of unity after communication
  // to vtk file and check against a reference. However this would
  // unnecessarily clutter the directory with reference vtks containing only
  // ones, so we iterate through the vector and check at runtime.
  for (const auto& entry : partUnity)
  {
    if (!Utility::is_numeric_one(entry))
    {
      DUNE_THROW(Dune::Exception,"Sanity check failed! "
                 "Partition of unity does not add up to 1.0.");
    }
  }

  std::cout << "Sanity check succeeded." << std::endl;
  return;
}


/*
 * The sarkisPartitionOfUnity method from geneo/partitionofunity.hh only works
 * in 2 dimensions. In 3 dimensions compilation fails due to invalid type
 * conversion. So, to avoid compiling the sarkis partition of unity in 3
 * dimensions, we split the test routine into two cases: The standard partition
 * of unity in standardTestdriver that works for any dimension, and the
 * sarkisTestdriver that is only called in the 2 dimensional case.
 */


template<class GV, class FEM>
void standardTestdriver(const GV& gv, const FEM& fem, const int ord,
                        const std::string& filename)
{
  using Dune::PDELab::Backend::native;
  using RangeField = double;

  using VBE = Dune::PDELab::ISTL::VectorBackend<>;

  // for partitions of unity, we need constraints on the subdomain boundaries
  using CON = Dune::PDELab::OverlappingConformingDirichletConstraints;

  using GFS = Dune::PDELab::GridFunctionSpace<GV, FEM, CON, VBE>;
  GFS gfs(gv, fem);

  using PartUnityCC =
    typename GFS::template ConstraintsContainer<RangeField>::Type;
  auto cc { PartUnityCC() };

  // exclude constraints at the domain boundary
  Dune::PDELab::NoDirichletConstraintsParameters bctype;

  // enforce constraints
  Dune::PDELab::constraints(bctype, gfs, cc);

  // create partition of unity
  using V = Dune::PDELab::Backend::Vector<GFS, RangeField>;
  auto partUnity { std::make_shared<V>(standardPartitionOfUnity<V>(gfs, cc)) };

  Utility::write_gridfunction_vtk(filename, gv, gfs, *partUnity, ord);
  perform_sanity_check(gfs, *partUnity);

  return;
}


template<class GV, class FEM, class GP>
void sarkisTestdriver(const GV& gv, const FEM& fem, const GP& gp,
                      const int ord, const std::string& filename)
{
  using Dune::PDELab::Backend::native;
  using RangeField = double;

  using VBE = Dune::PDELab::ISTL::VectorBackend<>;

  // for partitions of unity, we need constraints on the subdomain boundaries
  using CON = Dune::PDELab::OverlappingConformingDirichletConstraints;

  using GFS = Dune::PDELab::GridFunctionSpace<GV, FEM, CON, VBE>;
  GFS gfs(gv, fem);

  using PartUnityCC =
    typename GFS::template ConstraintsContainer<RangeField>::Type;
  auto cc { PartUnityCC() };

  // exclude constraints at the domain boundary
  Dune::PDELab::NoDirichletConstraintsParameters bctype;

  // enforce constraints
  Dune::PDELab::constraints(bctype, gfs, cc);

  // create an auxiliary local function space for the construction of the
  // Sarkis partition of unity
  using LFS =
    Dune::PDELab::LocalFunctionSpace<GFS, Dune::PDELab::AnySpaceTag>;
  LFS lfs(gfs);

  using V = Dune::PDELab::Backend::Vector<GFS, RangeField>;
  auto partUnity { std::make_shared<V>(sarkisPartitionOfUnity<V>(
                                         gfs, lfs, cc, gp.nCells.at(0),
                                         gp.nCells.at(1), gp.ovlp,
                                         gp.subdomLayout.at(0),
                                         gp.subdomLayout.at(1))) };

  Utility::write_gridfunction_vtk(filename, gv, gfs, *partUnity, ord);
  perform_sanity_check(gfs, *partUnity);

  return;
}


int main(int argc, char** argv)
{
  Dune::MPIHelper& mpihelper { Dune::MPIHelper::instance(argc, argv) };

  // open ini file
  Dune::ParameterTree ptree;

  // if single .ini file is passed as command line argument, use it
  if (argc == 2)
    Dune::ParameterTreeParser::readINITree(argv[1], ptree);
  // if no .ini file is passed as command line argument, use default
  else if (argc == 1)
    Dune::ParameterTreeParser::readINITree("test-geneo-partitionofunity.ini",
                                           ptree);
  else
    DUNE_THROW(Dune::InvalidStateException, "Error parsing .ini file. "
               "Either provide a single .ini file or no .ini file (default "
               "used).");

  // read static parameters and partition of unity type from ini file
  const int dim { ptree.get<int>("grid.dim") };
  const int ord { ptree.get<int>("partition.baseDegree") };
  const std::string partUnityType { ptree.get("partition.type", "") };

  std::string filename { generate_filename(ptree) };

  if ((partUnityType == "sarkis") && (dim == THREEDIM))
    DUNE_THROW(Dune::Exception, "Sarkis partition of unity is not yet "
               "supported in 3 dimensions.");

  using Real = double;
  using GVCoord = Real;

  if (dim == TWODIM)
  {
    Utility::GridParameters<GVCoord, TWODIM> gp(ptree);

    // grid partitioning
    using Partitioner = Dune::YaspFixedSizePartitioner<TWODIM>;
    auto partitioner { std::make_unique<Partitioner>(gp.subdomLayout) };

    // instantiate grid
    using Grid = Dune::YaspGrid<TWODIM>;
    auto grid { std::make_shared<Grid>(
                  gp.upperRight, gp.nCells, gp.isPeriodic, gp.ovlp,
                  Dune::MPIHelper::getCollectiveCommunication(),
                  partitioner.get()) };

    // create GridView to finest level
    using GV = Grid::LevelGridView;
    auto gv { grid->levelGridView(grid->maxLevel()) };

    if (ord == 1)
    {
      Dune::PDELab::QkLocalFiniteElementMap<GV, GVCoord, Real, 1> feMap(gv);

      if (partUnityType == "standard")
        standardTestdriver(gv, feMap, ord, filename);
      else if (partUnityType == "sarkis")
        sarkisTestdriver(gv, feMap, gp, ord, filename);
      else
        DUNE_THROW(Dune::Exception, "Unknown partition of unity type.");
    }
    else if (ord == 2)
    {
      Dune::PDELab::QkLocalFiniteElementMap<GV, GVCoord, Real, 2>  feMap(gv);

      if (partUnityType == "standard")
        standardTestdriver(gv, feMap, ord, filename);
      else if (partUnityType == "sarkis")
        sarkisTestdriver(gv, feMap, gp, ord, filename);
      else
        DUNE_THROW(Dune::Exception, "Unknown partition of unity type.");
    }
    else
      DUNE_THROW(Dune::Exception, "Degree higher than 3 not yet supported.");
  }
  else if (dim == THREEDIM)
  {
    Utility::GridParameters<GVCoord, THREEDIM> gp(ptree);

    // grid partitioning
    using Partitioner = Dune::YaspFixedSizePartitioner<THREEDIM>;
    auto partitioner { std::make_unique<Partitioner>(gp.subdomLayout) };

    // instantiate grid
    using Grid = Dune::YaspGrid<THREEDIM>;
    auto grid { std::make_shared<Grid>(
                  gp.upperRight, gp.nCells, gp.isPeriodic, gp.ovlp,
                  Dune::MPIHelper::getCollectiveCommunication(),
                  partitioner.get()) };

    // create GridView to finest level
    using GV = Grid::LevelGridView;
    auto gv { grid->levelGridView(grid->maxLevel()) };

    if (ord == 1)
    {
      Dune::PDELab::QkLocalFiniteElementMap<GV, GVCoord, Real, 1> feMap(gv);

      standardTestdriver(gv, feMap, ord, filename);
    }
    else if (ord == 2)
    {
      Dune::PDELab::QkLocalFiniteElementMap<GV, GVCoord, Real, 2> feMap(gv);

      standardTestdriver(gv, feMap, ord, filename);
    }
    else
      DUNE_THROW(Dune::Exception, "Degree higher than 2 not yet supported.");
  }

  return 0;
}
