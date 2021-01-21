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

#include <dune/pdelab/backend/istl/geneo/OfflineOnline/OfflineTools.hh>
#include <dune/pdelab/backend/istl/geneo/OfflineOnline/GenericEllipticProblem.hh>

void driver(Dune::MPIHelper& helper) {

  // ~~~~~~~~~~~~~~~~~~
  // Create offline hierarchy
  // ~~~~~~~~~~~~~~~~~~
  std::string path_to_storage = "Offline/";

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

  const int algebraic_overlap = 2;

  // Getting map global ID to global Index
  std::vector<long unsigned int> NodalGID = extractNodalGID<GRID>(grid);

  // Transfer partitioning from ParMETIS to our grid
#if PARMETIS_MAJOR_VERSION
  // ~~~~~~~~~~~~~~~~~~
  // Save domain decomposition (meshes)
  // ~~~~~~~~~~~~~~~~~~
  // First recreate the overlapped subdomain
  unsigned subdomains = grid->comm().size();
  auto partition = parmetis_partitioning(gv,helper,subdomains);
  unsigned uoverlap = algebraic_overlap + 0u;
  std::string extensionmethod = "vertex";
  auto overlappingsubdomains = grow_subdomains(gv,subdomains,partition,uoverlap,extensionmethod);

  // Write subdomain as a msh file
  Dune::SubDomainGmshWriter<GV> sdwriter(gv, overlappingsubdomains, subdomains, 16);
  sdwriter.write(path_to_storage, "subdomain.msh");

  grid->loadBalance(partition, 0);
#else
  grid->loadBalance();
#endif

  // ~~~~~~~~~~~~~~~~~~
//  Type definitions
  // ~~~~~~~~~~~~~~~~~~
  const int components = 1;
  using K = double;
  using Vector = Dune::BlockVector<Dune::FieldVector<K,components>>;
  using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<K,components,components>>;

  typedef Dune::BlockVector<Dune::FieldVector<K, 1>> CoarseVector;
  typedef Dune::BCRSMatrix<Dune::FieldMatrix<K, 1, 1>> CoarseMatrix;

  typedef Dune::BlockVector<Dune::FieldVector<int, 1>> vector1i;

  using ESExcluder = Dune::PDELab::EntitySetExcluder<Vector, GV>;
  auto ghost_excluder = std::make_shared<Dune::PDELab::EntitySetGhostExcluder<Vector, GV>>();

  using ES = Dune::PDELab::ExcluderEntitySet<GV,Dune::Partitions::All, ESExcluder>;
  ES es(gv, ghost_excluder);

  // make problem parameters
  typedef GenericEllipticProblem<ES,NumberType> Problem;
  Problem problem;
  typedef Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<Problem> BCType;
  BCType bctype(es,problem);

  using Dune::PDELab::Backend::native;

  // ~~~~~~~~~~~~~~~~~~
//  FE fine Space definition
  // ~~~~~~~~~~~~~~~~~~
  typedef typename ES::Grid::ctype DF;
  // instantiate finite element maps
  typedef Dune::PDELab::QkLocalFiniteElementMap<ES,DF,double,1> FEM;
  FEM fem(es);
  // function space with no constraints on processor boundaries, needed for the GenEO eigenproblem
  typedef Dune::PDELab::GridFunctionSpace<ES,FEM,
                                          Dune::PDELab::ConformingDirichletConstraints,
                                          Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed,1>
                                          > GFS;
  GFS gfs(es,fem);
  int verbose = 0;
  if (gfs.gridView().comm().rank()==0) verbose = 2;
  // make a degree of freedom vector on fine grid and initialize it with interpolation of Dirichlet condition
  typedef Dune::PDELab::Backend::Vector<GFS,NumberType> V;
  V x(gfs,0.0);
  // Extract domain boundary constraints from problem definition, apply trace to solution vector
  typedef Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<Problem> G;
  G g(es,problem);
  Dune::PDELab::interpolate(g,gfs,x);
  // Set up constraints containers with boundary constraints, but without processor constraints
  typedef typename GFS::template ConstraintsContainer<NumberType>::Type CC;
  auto cc = CC();
  // assemble constraints
  Dune::PDELab::constraints(bctype,gfs,cc);
  // set initial guess
  V x0(gfs,0.0);
  Dune::PDELab::copy_nonconstrained_dofs(cc,x0,x);
  // LocalOperator for given problem
  typedef Dune::PDELab::ConvectionDiffusionFEM<Problem,FEM> LOP;
  LOP lop(problem);
  // LocalOperator wrapper zeroing out subdomains' interiors in order to set up overlap matrix
  typedef Dune::PDELab::ISTL::BCRSMatrixBackend<> MBE;
  // Construct GridOperators from LocalOperators
  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,NumberType,NumberType,NumberType,CC,CC> GO;
  auto go = GO(gfs,cc,gfs,cc,lop,MBE(nonzeros));
  // Assemble fine grid matrix defined without processor constraints
  typedef typename GO::Jacobian M;
  M A(go);
  go.jacobian(x,A);
  // set up and assemble right hand side w.r.t. l(v)-a(u_g,v)
  V d(gfs,0.0);
  go.residual(x,d);

  // ~~~~~~~~~~~~~~~~~~
//  Solving process begin here: First some parameters
  // ~~~~~~~~~~~~~~~~~~
  double eigenvalue_threshold = -1;
  int nev = 4;
  int nev_arpack = nev;
  // double shift = 0.001;

  // ~~~~~~~~~~~~~~~~~~
//  ADAPTER
  // ~~~~~~~~~~~~~~~~~~
  Dune::NonoverlappingOverlapAdapter<GV, Vector, Matrix> adapter(gv, native(A), nonzeros, algebraic_overlap);
  using Attribute = Dune::EPISAttribute;
  Dune::AllSet<Attribute> allAttribute;
  auto allinterface = std::shared_ptr<Dune::Interface>(new Dune::Interface());
  allinterface->build(*adapter.getRemoteIndices(), allAttribute, allAttribute); // all to all communication
  // build up buffered communicator allowing communication over a dof vector
  auto communicator = std::shared_ptr<Dune::BufferedCommunicator>(new Dune::BufferedCommunicator());
  communicator->build<Vector>(*allinterface);

  auto geneo_matrices = setupGenEOMatrices(go, adapter, A);
  std::shared_ptr<Matrix> A_extended = std::get<0>(geneo_matrices);
  std::shared_ptr<Matrix> A_overlap_extended = std::get<1>(geneo_matrices);
  std::shared_ptr<Vector> part_unity = std::get<2>(geneo_matrices);

  auto rank = adapter.gridView().comm().rank();

  /* Plot partition of unity for comparison with online runs */
  V vect(gfs, 0.0);
  adapter.restrictVector(*part_unity, native(vect));
  Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,Dune::refinementLevels(0));
  Dune::PDELab::vtk::DefaultFunctionNameGenerator fieldname;
  fieldname.prefix("PoU");
  Dune::PDELab::addSolutionToVTKWriter(vtkwriter,gfs,vect,fieldname);
  vtkwriter.write("PoU",Dune::VTK::ascii);

  // ~~~~~~~~~~~~~~~~~~
// Store new2old / old2new index -> potentially useful to post-process online calculus
  // ~~~~~~~~~~~~~~~~~~

  auto new2old_localindex = adapter.get_new2old_localindex();
  auto old2new_localindex = adapter.get_old2new_localindex();
  vector1i n2o(new2old_localindex.size()), o2n(old2new_localindex.size());
  for (int i=0; i<new2old_localindex.size();i++) { n2o[i]=new2old_localindex[i]; }
  for (int i=0; i<old2new_localindex.size();i++) { o2n[i]=old2new_localindex[i]; }

  std::string filename_n2o = path_to_storage + std::to_string(rank) + "_restrict.mm";
  Dune::storeMatrixMarket(n2o, filename_n2o);
  std::string filename_o2n = path_to_storage + std::to_string(rank) + "_extend.mm";
  Dune::storeMatrixMarket(o2n, filename_o2n);

  // ~~~~~~~~~~~~~~~~~~
// Get and store local index to global id
  // ~~~~~~~~~~~~~~~~~~
  vector1i GlobalIndices = extractGlobalIndices<vector1i, GV, Matrix, Vector>(adapter, NodalGID);
  std::string filename_GI = path_to_storage + std::to_string(rank) + "_GI.mm";
  Dune::storeMatrixMarket(GlobalIndices, filename_GI);

  // ~~~~~~~~~~~~~~~~~~
// Store PoU
  // ~~~~~~~~~~~~~~~~~~
  Vector PoU((*part_unity).size());
  for (int i=0; i<(*part_unity).size();i++)
    PoU[i]=(*part_unity)[i];
  std::string filename_PoU = path_to_storage + std::to_string(rank) + "_PoU.mm";
  Dune::storeMatrixMarket(PoU, filename_PoU, 16);

  // ~~~~~~~~~~~~~~~~~~
// Save neighbor ranks vector
  // ~~~~~~~~~~~~~~~~~~
  std::vector<int> neighbor_ranks = adapter.findNeighboringRanks();
  vector1i NR(neighbor_ranks.size());
  for (int i=0; i<neighbor_ranks.size();i++)
    NR[i]=neighbor_ranks[i];
  std::string filename_NR = path_to_storage + std::to_string(rank) + "_neighborRanks.mm";
  Dune::storeMatrixMarket(NR, filename_NR);

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  Subdomain basis computations
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  std::shared_ptr<Dune::PDELab::SubdomainBasis<Vector>> subdomainbasis;
  subdomainbasis = std::make_shared<Dune::PDELab::NonoverlappingGenEOBasis<GO, Matrix, Vector>>(adapter, A_extended, A_extended, part_unity, eigenvalue_threshold, nev, nev_arpack);

  //  Particular solution
  // ParticularSolution PartSol = ParticularSolution(adapter, A_extended, *part_unity);
  // PartSol.exactRHS(native(d));
  // PartSol.solveAndAppend(*subdomainbasis);
  // TODO :: compute the particular solution just on certain subdomains not all

  // ~~~~~~~~~~~~~~~~~~
  // Save local basis sizes
  // ~~~~~~~~~~~~~~~~~~
  // Communicate local coarse space dimensions
  vector1i local_basis_sizes(gv.comm().size());
  int local_size = subdomainbasis->basis_size();
  gv.comm().allgather(&local_size, 1, local_basis_sizes.data());
  if (rank==0) {
    std::string filename_ls = path_to_storage + "localBasisSizes.mm";
    Dune::storeMatrixMarket(local_basis_sizes, filename_ls);
  }

  // ~~~~~~~~~~~~~~~~~~
  // Save local basis
  // ~~~~~~~~~~~~~~~~~~
  for (int basis_index = 0; basis_index < subdomainbasis->basis_size(); basis_index++) {
    std::string filename_EV = path_to_storage + std::to_string(rank) + "_EV_" + std::to_string(basis_index) + ".mm";
    Dune::storeMatrixMarket(*subdomainbasis->get_basis_vector(basis_index), filename_EV, 16);
  }

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  Coarse space construction and saving
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  auto coarse_space = std::make_shared<Dune::PDELab::NonoverlappingSubdomainProjectedCoarseSpace<GV, Matrix, Vector>>(adapter, gv, *A_extended, subdomainbasis, verbose);
  // Get the coarse matrix from the coarse space
  std::shared_ptr<CoarseMatrix> AH = coarse_space->get_coarse_system();
  // Save the coarse space
  std::string filename = path_to_storage + "OfflineAH.mm";
  Dune::storeMatrixMarket(*AH, filename, 16);

  // ~~~~~~~~~~~~~~~~~~
//  Coarse space right hand side and saving
  // ~~~~~~~~~~~~~~~~~~
  // Use the correct rhs to solve the final system
  Vector b(adapter.getExtendedSize());
  adapter.extendVector(native(d), b);
  communicator->forward<AddGatherScatter<Vector>>(b,b); // make function known in other subdomains

  Vector b_offline(b.N());
  for (int i=0; i<b.N();i++)
    b_offline[i]=b[i];
  std::string filename_fineb = path_to_storage + std::to_string(gv.comm().rank()) + "_fineb.mm";
  Dune::storeMatrixMarket(b_offline, filename_fineb, 16);

  // if(gv.comm().rank()==1){
  //   std::cout << "Offline b: " << std::endl;
  //   std::cout << "[" << b[0];
  //   for (int i = 1; i < b.N(); i++)
  //     std::cout << "  " <<  b[i];
  //   std::cout << "]" << std::endl;
  // }

  // Initializate a load vector (right hand side) in the coarse space : coarse_d = RH * b
  CoarseVector coarse_d(coarse_space->basis_size(), coarse_space->basis_size());
  coarse_space->restrict(b,coarse_d);
  // Save the coarse rhs
  if (rank==0){
    std::string filename_cb = path_to_storage + "OfflinebH.mm";
    Dune::storeMatrixMarket(coarse_d, filename_cb, 16);
  }

  // ~~~~~~~~~~~~~~~~~~
// Solve the coarse space system  ::TODO:: MAKE THAT OPTIONAL for a normal offline go
  // ~~~~~~~~~~~~~~~~~~
  // Objective :  find x = RH^T * AH^-1 * RH * b
  // Use of UMFPack to solve the problem [AH * coarse_v = RH * b] instead of inversing AH :: objective is to have [coarse_v = AH^-1 *  RH * b]
  Dune::UMFPack<CoarseMatrix> coarse_solver(*AH, false); // Is there something better than UMFPACK?
  CoarseVector coarse_v(coarse_space->basis_size(), coarse_space->basis_size());
  Dune::InverseOperatorResult result;
  coarse_solver.apply(coarse_v,coarse_d,result);

  // Save the coarse solution to further comparaison online offline
  std::string filename_coarse_sol = path_to_storage + "pristineSol.mm";
  Dune::storeMatrixMarket(coarse_v, filename_coarse_sol, 16);

  // ~~~~~~~~~~~~~~~~~~
// Prolongate the solution in order to have v_fine_vovlp on the fine space : v_fine_vovlp = RH^T * coarse_v = RH^T * AH^-1 * RH * b
  // ~~~~~~~~~~~~~~~~~~
  Vector v_fine_ovlp(adapter.getExtendedSize());
  coarse_space->prolongate(coarse_v, v_fine_ovlp);

  Vector fineSolOffline(v_fine_ovlp.N());
  for (int i=0; i<v_fine_ovlp.N();i++)
    fineSolOffline[i]=v_fine_ovlp[i];
  std::string filename_fineSolSubdomainOnly = path_to_storage + std::to_string(gv.comm().rank()) + "_fineSolSubdomainOnly.mm";
  Dune::storeMatrixMarket(fineSolOffline, filename_fineSolSubdomainOnly, 15);

  communicator->forward<AddGatherScatter<Vector>>(v_fine_ovlp, v_fine_ovlp);
  V v(gfs, 0.0);
  adapter.restrictVector(v_fine_ovlp, v);

  // native(x) -= v;

  typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
  typedef Dune::PDELab::VTKGridFunctionAdapter<DGF> ADAPT;

  /* Write solution to VTK */
  Dune::VTKWriter<GV> vtkwriterSol(gv);
  DGF xdgfSol(gfs,v);
  auto adaptSol = std::make_shared<ADAPT>(xdgfSol,"solution");
  vtkwriterSol.addVertexData(adaptSol);
  vtkwriterSol.write("Offline_solution");

  // ~~~~~~~~~~~~~~~~~~
// Visualise all the basis in vtk format
  // ~~~~~~~~~~~~~~~~~~
  for (int basis_index = 0; basis_index < subdomainbasis->basis_size(); basis_index++) {
    Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,Dune::refinementLevels(0));
    V vect(gfs, 0.0);
    adapter.restrictVector(native(*subdomainbasis->get_basis_vector(basis_index)), native(vect));

    int rank = adapter.gridView().comm().rank();
    std::string filename = "BasisVector_"+std::to_string(basis_index);

    Dune::PDELab::vtk::DefaultFunctionNameGenerator fieldname;
    fieldname.prefix("EV");

    Dune::PDELab::addSolutionToVTKWriter(vtkwriter,gfs,vect,fieldname);
    vtkwriter.write(filename,Dune::VTK::ascii);
  }
}


int main(int argc, char **argv)
{
  using Dune::PDELab::Backend::native;

  // initialize MPI, finalize is done automatically on exit
  Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc,argv);

  driver(helper);

  return 0;
}
