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
//  Type definitions
  // ~~~~~~~~~~~~~~~~~~
  const int components = 1;
  using K = double;
  using Vector = Dune::BlockVector<Dune::FieldVector<K,components>>;
  using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<K,components,components>>;
  typedef Dune::BlockVector<Dune::FieldVector<K, 1>> CoarseVector;
  typedef Dune::BCRSMatrix<Dune::FieldMatrix<K, 1, 1>> CoarseMatrix;
  typedef Dune::BlockVector<Dune::FieldVector<int, 1>> vector1i;

  const unsigned int dim = 2;
  const unsigned int degree = 1;
  const std::size_t nonzeros = std::pow(2*degree+1,dim);
  typedef double NumberType;

  using Dune::PDELab::Backend::native;

  // Load the coarse matrix
  CoarseMatrix AH = FromOffline<CoarseMatrix>(path_to_storage, "OfflineAH");

  for (int row_id = 0; row_id < AH.N(); row_id++){
    double value = AH[row_id][0];
    if (std::abs(value) > 1e-6){std::cout << value;}
    else{std::cout << 0.0;}
    for (int i = 1; i < AH.M(); i++){
      value = AH[row_id][i];
      if (std::abs(value) > 1e-6){std::cout << ", " << value;}
      else {std::cout << ", " << 0.0;}
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  // Load the coarse vector
  CoarseVector bH = FromOffline<CoarseVector>(path_to_storage, "OfflinebH");

  std::cout << "[" << bH[0];
  for (int i = 1; i < bH.N(); i++){
    std::cout << "  " <<  bH[i];
  }
  std::cout << "]" << std::endl;
  std::cout << std::endl;


  // Offline localBasisSizes
  vector1i local_basis_sizes = FromOffline<vector1i>(path_to_storage, "localBasisSizes");
  const int number_of_rank_used_offline = local_basis_sizes.size();
  std::string DefectID = "Defect1";
  vector1i changes_in_local_basis_sizes = FromOffline<vector1i>(path_to_storage, DefectID+"_changeBasisSizes");


  for( const auto Target : targeted) {
    std::cout << "Updating... with subdomain " << Target << " contribution." << std::endl;

    //  Grid set up
    typedef Dune::UGGrid<dim> GRID;
    Dune::GridFactory<GRID> factory;
    std::vector<int> gmsh2dune = Dune::SubDomainGmshReader<GRID>::read_and_return(factory, path_to_storage + std::to_string(Target) + "_subdomain.msh", true, false);
    std::unique_ptr<GRID> grid (factory.createGrid());

    typedef typename GRID::LeafGridView GV;
    auto gv = grid->leafGridView();

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

    std::size_t v_size = A.N();

    /* Define vtk writer utils */
    typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
    typedef Dune::PDELab::VTKGridFunctionAdapter<DGF> ADAPT;

    // ~~~~~~~~~~~~~~~~~~
  //  Load a vector describing local basis sizes (number of EV) and creating the vector of offsets (to reach indices in the coarse space)
    // ~~~~~~~~~~~~~~~~~~
    local_basis_sizes[Target] += changes_in_local_basis_sizes[Target];

    vector1i local_offset(number_of_rank_used_offline+1);
    local_offset[0]=0;
    for (int i=0; i<local_basis_sizes.size();i++) {
      local_offset[i+1] = local_basis_sizes[i]+local_offset[i];
    }

    // ~~~~~~~~~~~~~~~~~~
    // Load the indices transformation
    // ~~~~~~~~~~~~~~~~~~
    vector1i offlineDoF2GI = FromOffline<vector1i>(path_to_storage, "GI", Target);
    vector1i DofOffline_to_DofOnline = offlineDoF2GI2gmsh2onlineDoF<vector1i>(Target, gmsh2dune, offlineDoF2GI, path_to_storage);

    // ~~~~~~~~~~~~~~~~~~
    // Load recomputed subdomain basis of the targeted subdomain
    // ~~~~~~~~~~~~~~~~~~
    std::shared_ptr<Dune::PDELab::SubdomainBasis<Vector>> online_subdomainbasis;
    online_subdomainbasis = std::make_shared<Dune::PDELab::GenEOBasisFromFiles<GO, Matrix, Vector, vector1i>>(path_to_storage, local_basis_sizes[Target], Target, DofOffline_to_DofOnline, DefectID);

    // ~~~~~~~~~~~~~~~~~~
    // Load PoU
    // ~~~~~~~~~~~~~~~~~~
    Vector PoU = FromOffline<Vector>(path_to_storage, "PoU", Target);
    Vector nPoU(v_size);
    for (int i=0; i<v_size; i++){
      nPoU[DofOffline_to_DofOnline[i]] = PoU[i];
    }

    // ~~~~~~~~~~~~~~~~~~
    // Load Neighbour ranks
    // ~~~~~~~~~~~~~~~~~~
    vector1i NR = FromOffline<vector1i>(path_to_storage, "neighborRanks", Target);

    std::vector<std::shared_ptr<Dune::PDELab::SubdomainBasis<Vector>>> neighbour_subdomainbasis(NR.size());
    for (int iter_over_subdomains=0; iter_over_subdomains<NR.size(); iter_over_subdomains++) {

      int subdomainNumber = NR[iter_over_subdomains];
      vector1i offlineNeighbourDoF2GI = FromOffline<vector1i>(path_to_storage, "GI", subdomainNumber);

      auto iterator = std::find(targeted.begin(), targeted.end(), subdomainNumber);
      if(iterator!=targeted.end()){
        neighbour_subdomainbasis[iter_over_subdomains] = std::make_shared<Dune::PDELab::NeighbourBasis<GO, Matrix, Vector, vector1i>>(path_to_storage, local_basis_sizes[NR[iter_over_subdomains]], NR[iter_over_subdomains], offlineDoF2GI, offlineNeighbourDoF2GI, DefectID);
      } else {
        neighbour_subdomainbasis[iter_over_subdomains] = std::make_shared<Dune::PDELab::NeighbourBasis<GO, Matrix, Vector, vector1i>>(path_to_storage, local_basis_sizes[NR[iter_over_subdomains]], NR[iter_over_subdomains], offlineDoF2GI, offlineNeighbourDoF2GI);
      }

      for (int i=0; i<neighbour_subdomainbasis[iter_over_subdomains]->basis_size(); i++){
        auto tmp = *neighbour_subdomainbasis[iter_over_subdomains]->get_basis_vector(i);
        for (int j=0; j<v_size; j++){
          (*neighbour_subdomainbasis[iter_over_subdomains]->get_basis_vector(i))[DofOffline_to_DofOnline[j]] = tmp[j];
        }
      }
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //  Update AH & bH
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    int delta_basis_size = changes_in_local_basis_sizes[Target];

    if (delta_basis_size==0){
      UpdateAH<Vector, Matrix, CoarseMatrix, vector1i>(AH, native(A), online_subdomainbasis, neighbour_subdomainbasis, local_basis_sizes, local_offset, NR, Target);
    } else {
      UpdateAHNewSize<Vector, Matrix, CoarseMatrix, vector1i>(AH, native(A), online_subdomainbasis, neighbour_subdomainbasis, local_basis_sizes, local_offset, Target, delta_basis_size, path_to_storage);
    }


    for (int row_id = 0; row_id < AH.N(); row_id++){
      double value = AH[row_id][0];
      if (std::abs(value) > 1e-6){std::cout << value;}
      else{std::cout << 0.0;}
      for (int i = 1; i < AH.M(); i++){
        value = AH[row_id][i];
        if (std::abs(value) > 1e-6){std::cout << ", " << value;}
        else {std::cout << ", " << 0.0;}
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;

    if (delta_basis_size==0){
      UpdatebH<Vector, CoarseVector, vector1i>(bH, native(fine_b), online_subdomainbasis, local_basis_sizes, local_offset, Target);
    } else {
      UpdatebHNewSize<Vector, CoarseVector, vector1i>(bH, native(fine_b), online_subdomainbasis, local_basis_sizes, local_offset, Target, delta_basis_size);
    }

    std::cout << "[" << bH[0];
    for (int i = 1; i < bH.N(); i++){
      std::cout << "  " <<  bH[i];
    }
    std::cout << "]" << std::endl;
    std::cout << std::endl;

    int global_basis_size = std::accumulate(local_basis_sizes.begin(), local_basis_sizes.end(), 0.0);

    // Objective :  find x = RH^T * AH^-1 * RH * b
    // Use of UMFPack to solve the problem [AH * coarse_v = RH * b] instead of inversing AH :: objective is to have [coarse_v = AH^-1 *  RH * b]
    Dune::UMFPack<CoarseMatrix> coarse_solver(AH, false); // Is there something better than UMFPACK?
    CoarseVector coarse_sol(global_basis_size, global_basis_size);
    Dune::InverseOperatorResult result;
    coarse_solver.apply(coarse_sol,bH,result);

    // Plot a part of the solution over the online domain
    V prolongated(gfs,0.0);
    // Prolongate result
    for (int basis_index = 0; basis_index < local_basis_sizes[Target]; basis_index++) {
      Vector local_result(*online_subdomainbasis->get_basis_vector(basis_index));
      // Vector local_result(*offline_subdomainbasis->get_basis_vector(basis_index));
      local_result *= coarse_sol[local_offset[Target] + basis_index];
      native(prolongated) += local_result;
    }

    Dune::VTKWriter<GV> vtkwriter(gv);
    /* Write a field in the vtu file */
    DGF xdgf(gfs,prolongated);
    auto adapt = std::make_shared<ADAPT>(xdgf,"Solution");
    vtkwriter.addVertexData(adapt);
    vtkwriter.write("onlineSolSub-"+std::to_string(Target),Dune::VTK::ascii);

  }

}


int main(int argc, char **argv)
{
  using Dune::PDELab::Backend::native;

  // Offline folder
  std::string path_to_storage = "Offline/";

  // Define what subdomain need to be solved
  std::vector<int> targeted = {0,1}; // Subdomains that need a second solve

  // initialize MPI, finalize is done automatically on exit
  Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc,argv);

  driver(path_to_storage, targeted, helper);

  return 0;
}
