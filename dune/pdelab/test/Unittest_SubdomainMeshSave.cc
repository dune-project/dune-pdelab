#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/pdelab.hh>

#include <dune/grid/utility/globalindexset.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/common/rangegenerators.hh>
#include <dune/grid/common/gridinfo.hh> // visualize grid information -> could be remove at the end


#include <dune/grid/utility/parmetisgridpartitioner.hh>
#include <dune/pdelab/backend/istl/geneo/OfflineOnline/SubDomainGmshWriter.hh>
#include <dune/pdelab/backend/istl/geneo/OfflineOnline/partitioner.hh>
#include <dune/pdelab/backend/istl/geneo/OfflineOnline/GenericEllipticProblem.hh>

void driver(Dune::MPIHelper& helper) {

  // ~~~~~~~~~~~~~~~~~~
//  Grid set up
  // ~~~~~~~~~~~~~~~~~~
  // define parameters
  const unsigned int dim = 2;
  const unsigned int degree = 1;
  const std::size_t nonzeros = std::pow(2*degree+1,dim);
  typedef double NumberType;

  typedef Dune::UGGrid<dim> GRID;
  Dune::GridFactory<GRID> factory;
  Dune::GmshReader<GRID>::read(factory, "grids/24x24.msh", true, true);
  std::unique_ptr<GRID> grid (factory.createGrid());

  typedef typename GRID::LeafGridView GV;
  auto gv = grid->leafGridView();

  // std::vector<unsigned> partition(grid->leafGridView().size(0));
  // std::cout << grid->leafGridView().size(0) << std::endl;
  // for (int i=0; i< partition.size(); i++) {
  //   if (i < partition.size()/2 ) {partition[i]=0;}
  //   else {partition[i]=1;}
  // }

  unsigned subdomains = 3; // grid->comm().size();
  unsigned uoverlap = 1u;
  std::string extensionmethod = "vertex";

  auto partition = parmetis_partitioning(gv,helper,subdomains);

  // ~~~~~~~~~~~~~~~~~~
  // Save domain decomposition
  // ~~~~~~~~~~~~~~~~~~
  // First recreate the overlapped subdomain
  auto overlappingsubdomains = grow_subdomains(gv,subdomains,partition,uoverlap,extensionmethod);

  // test overlap generation
  auto& indexset = gv.indexSet(); // to attach data to elements
  std::cout << indexset.size(0) << " elements" << std::endl;
  std::cout << indexset.size(dim) << " vertices" << std::endl;
  std::cout << subdomains << " subdomains" << std::endl;
  int vp = 0;
  std::vector<unsigned> adomain(indexset.size(0));
  std::vector<unsigned> k0(indexset.size(0));
  size_t k_0 = 0;
  unsigned count_elements=0;
  for (const auto& e : elements(gv))
  {
    adomain[indexset.index(e)] = 0;
    if (overlappingsubdomains[indexset.index(e)].count(vp)>0) adomain[indexset.index(e)] = 1;
    if (partition[indexset.index(e)]==vp) adomain[indexset.index(e)] = 2;
    k0[indexset.index(e)] = overlappingsubdomains[indexset.index(e)].size();
    k_0 = std::max(k_0,overlappingsubdomains[indexset.index(e)].size());
    count_elements += overlappingsubdomains[indexset.index(e)].size();
  }
  std::cout << k_0 << " k_0" << std::endl;
  std::cout << indexset.size(0) << " elements, " << count_elements << " in all subdomains, ratio=" << ((double)count_elements)/indexset.size(0) << std::endl;

  // Write subdomain as a msh file
  Dune::SubDomainGmshWriter<GV> sdwriter(gv, overlappingsubdomains, subdomains);
  sdwriter.write("", "subdomain.msh");
}

int main(int argc, char **argv)
{
  using Dune::PDELab::Backend::native;

  // initialize MPI, finalize is done automatically on exit
  Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc,argv);

  driver(helper);

  return 0;
}
