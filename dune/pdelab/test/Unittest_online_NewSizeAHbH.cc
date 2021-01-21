#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/pdelab.hh>

#include <dune/pdelab/backend/istl/geneo/OfflineOnline/geneobasisOnline.hh>
#include <dune/pdelab/backend/istl/geneo/OfflineOnline/OnlineTools.hh>
#include <dune/pdelab/backend/istl/geneo/OfflineOnline/SubDomainGmshReader.hh>
#include <dune/pdelab/backend/istl/geneo/OfflineOnline/GenericEllipticProblem.hh>

void driver(std::string path_to_storage, std::vector<int> targeted, Dune::MPIHelper& helper) {

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
  std::vector<int> gmsh2dune = Dune::SubDomainGmshReader<GRID>::read_and_return(factory, path_to_storage + std::to_string(targeted[0]) + "_subdomain.msh", true, false);
  std::unique_ptr<GRID> grid (factory.createGrid());

  typedef typename GRID::LeafGridView GV;
  auto gv = grid->leafGridView();

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
  typedef GenericEllipticProblem<ES,K> Problem;
  Problem problem;
  typedef Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<Problem> BCType;
  BCType bctype(es,problem);

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
  typedef Dune::PDELab::Backend::Vector<GFS,K> V;
  V x(gfs,0.0);
  // Extract domain boundary constraints from problem definition, apply trace to solution vector
  typedef Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<Problem> G;
  G g(es,problem);
  Dune::PDELab::interpolate(g,gfs,x);
  // Set up constraints containers with boundary constraints, but without processor constraints
  typedef typename GFS::template ConstraintsContainer<K>::Type CC;
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
  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,K,K,K,CC,CC> GO;
  auto go = GO(gfs,cc,gfs,cc,lop,MBE(nonzeros));
  // Assemble fine grid matrix defined without processor constraints
  typedef typename GO::Jacobian M;
  M A(go);
  go.jacobian(x,A);
  // set up and assemble right hand side w.r.t. l(v)-a(u_g,v)
  V fine_b(gfs,0.0);
  go.residual(x,fine_b);

  // ~~~~~~~~~~~~~~~~~~
//  Solving process begin here: First some parameters
  // ~~~~~~~~~~~~~~~~~~
  double eigenvalue_threshold = -1;
  // const int algebraic_overlap = 0;
  int nev = 4;
  int nev_arpack = nev;
  // double shift = 0.001;

  using Dune::PDELab::Backend::native;

  std::size_t v_size = A.N();

  /* Define vtk writer utils */
  typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
  typedef Dune::PDELab::VTKGridFunctionAdapter<DGF> ADAPT;

  // ~~~~~~~~~~~~~~~~~~
//  Load a vector describing local basis sizes (number of EV) and creating the vector of offsets (to reach indices in the coarse space)
  // ~~~~~~~~~~~~~~~~~~
  vector1i local_basis_sizes = FromOffline<vector1i>(path_to_storage, "localBasisSizes");
  const int number_of_rank_used_offline = local_basis_sizes.size();

  vector1i local_offset(number_of_rank_used_offline+1);
  local_offset[0]=0;
  for (int i=0; i<local_basis_sizes.size();i++) {
    local_offset[i+1] = local_basis_sizes[i]+local_offset[i];
  }

  // ~~~~~~~~~~~~~~~~~~
  // Load the indices transformation
  // ~~~~~~~~~~~~~~~~~~
  vector1i offlineDoF2GI = FromOffline<vector1i>(path_to_storage, "GI", targeted[0]);
  vector1i DofOffline_to_DofOnline = offlineDoF2GI2gmsh2onlineDoF<vector1i>(targeted[0], gmsh2dune, offlineDoF2GI, path_to_storage);

  // ~~~~~~~~~~~~~~~~~~
  // Load Neighbour ranks
  // ~~~~~~~~~~~~~~~~~~
  vector1i NR = FromOffline<vector1i>(path_to_storage, "neighborRanks", targeted[0]);

  // ~~~~~~~~~~~~~~~~~~
  // Load PoU
  // ~~~~~~~~~~~~~~~~~~
  Vector PoU = FromOffline<Vector>(path_to_storage, "PoU", targeted[0]);
  Vector nPoU(v_size);
  for (int i=0; i<v_size; i++){
    nPoU[DofOffline_to_DofOnline[i]] = PoU[i];
  }

  // ~~~~~~~~~~~~~~~~~~
//  Subdomain basis computation and loading for neighbours
  // ~~~~~~~~~~~~~~~~~~

  // First : compute the new targeted subdomain basis
  std::shared_ptr<Dune::PDELab::SubdomainBasis<Vector>> online_subdomainbasis;
  online_subdomainbasis = std::make_shared<Dune::PDELab::GenEOBasisOnline<GO, Matrix, Vector>>(native(A), nPoU, eigenvalue_threshold, nev, nev_arpack);

  // ~~~~~~~~~~~~~~~~~~
//  Particular solution
  // ~~~~~~~~~~~~~~~~~~
  auto PartSol = ParticularSolution<Vector, Matrix>(native(A));
  PartSol.exactRHS(native(fine_b));
  PartSol.solveAndAppend(*online_subdomainbasis);


  int delta_basis_size = online_subdomainbasis->basis_size() - local_basis_sizes[targeted[0]];

  // Then : load other subdomain basis from offline and transfer them in the targeted subdomain space
  std::vector<std::shared_ptr<Dune::PDELab::SubdomainBasis<Vector>>> neighbour_subdomainbasis(NR.size());

  for (int iter_over_subdomains=0; iter_over_subdomains<NR.size(); iter_over_subdomains++) {

    int int_Nnumber = NR[iter_over_subdomains];
    vector1i offlineNeighbourDoF2GI = FromOffline<vector1i>(path_to_storage, "GI", int_Nnumber);

    neighbour_subdomainbasis[iter_over_subdomains] = std::make_shared<Dune::PDELab::NeighbourBasis<GO, Matrix, Vector, vector1i>>(path_to_storage, local_basis_sizes[NR[iter_over_subdomains]], NR[iter_over_subdomains], offlineDoF2GI, offlineNeighbourDoF2GI);

    for (int i=0; i<local_basis_sizes[iter_over_subdomains]; i++){
      auto tmp = *neighbour_subdomainbasis[iter_over_subdomains]->get_basis_vector(i);
      for (int j=0; j<v_size; j++){
        (*neighbour_subdomainbasis[iter_over_subdomains]->get_basis_vector(i))[DofOffline_to_DofOnline[j]] = tmp[j];
      }
    }
  }

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  Update AH & bH
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // Load the coarse matrix
  CoarseMatrix AH = FromOffline<CoarseMatrix>(path_to_storage, "OfflineAH");
  UpdateAHNewSize<Vector, Matrix, CoarseMatrix, vector1i>(AH, native(A), online_subdomainbasis, neighbour_subdomainbasis, local_basis_sizes, local_offset, targeted[0], delta_basis_size, path_to_storage);

  // Load the coarse vector
  CoarseVector bH = FromOffline<CoarseVector>(path_to_storage, "OfflineCoarseb");
  UpdatebHNewSize<Vector, CoarseVector>(bH, native(fine_b), online_subdomainbasis, local_basis_sizes, local_offset, targeted[0], delta_basis_size);

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  Coarse space solve
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  int global_basis_size = std::accumulate(local_basis_sizes.begin(), local_basis_sizes.end(), 0.0) + delta_basis_size;

  // Objective :  find x = RH^T * AH^-1 * RH * b
  // Use of UMFPack to solve the problem [AH * coarse_v = RH * b] instead of inversing AH :: objective is to have [coarse_v = AH^-1 *  RH * b]
  Dune::UMFPack<CoarseMatrix> coarse_solver(AH, false); // Is there something better than UMFPACK?
  CoarseVector coarse_sol(global_basis_size, global_basis_size);
  Dune::InverseOperatorResult result;
  coarse_solver.apply(coarse_sol,bH,result);

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  Saves
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // Save the coarse solution to further prolongate it over the offline domaine in post-processing
  std::string filename_coarse_sol = path_to_storage + "coarse-sol_w_subdomain-" + std::to_string(targeted[0]) + "_recomputed.mm";
  Dune::storeMatrixMarket(coarse_sol, filename_coarse_sol, 15);

  // Plot a part of the solution over the online domain
  V prolongated(gfs,0.0);
  // Prolongate result
  for (int basis_index = 0; basis_index < local_basis_sizes[targeted[0]] + delta_basis_size; basis_index++) {
    Vector local_result(*online_subdomainbasis->get_basis_vector(basis_index));
    // Vector local_result(*offline_subdomainbasis->get_basis_vector(basis_index));
    local_result *= coarse_sol[local_offset[targeted[0]] + basis_index];
    native(prolongated) += local_result;
  }

  Dune::VTKWriter<GV> vtkwriter(gv);
  /* Write a field in the vtu file */
  DGF xdgf(gfs,prolongated);
  auto adapt = std::make_shared<ADAPT>(xdgf,"Solution");
  vtkwriter.addVertexData(adapt);
  vtkwriter.write("onlineSol",Dune::VTK::ascii);

  // TODO find a good test to do
}


int main(int argc, char **argv)
{
  using Dune::PDELab::Backend::native;

  // Offline folder
  std::string path_to_storage = "Offline/";

  // Define what subdomain need to be solved
  std::vector<int> targeted = {1}; // Subdomains that need a second solve

  // initialize MPI, finalize is done automatically on exit
  Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc,argv);

  driver(path_to_storage, targeted, helper);

  return 0;
}
