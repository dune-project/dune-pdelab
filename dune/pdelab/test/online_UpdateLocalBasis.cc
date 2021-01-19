#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/pdelab.hh>

#include <dune/pdelab/backend/istl/geneo/OfflineOnline/geneobasisOnline.hh>
#include <dune/pdelab/backend/istl/geneo/OfflineOnline/OnlineTools.hh>
#include <dune/pdelab/backend/istl/geneo/OfflineOnline/SubDomainGmshReader.hh>
#include <dune/pdelab/backend/istl/geneo/OfflineOnline/GenericEllipticProblem.hh>

void driver(std::string path_to_storage, std::vector<int> targeted, Dune::MPIHelper& helper) {

  using Dune::PDELab::Backend::native;

  //  Type definitions
  const int components = 1;
  using K = double;
  using Vector = Dune::BlockVector<Dune::FieldVector<K,components>>;
  using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<K,components,components>>;
  typedef Dune::BlockVector<Dune::FieldVector<K, 1>> CoarseVector;
  typedef Dune::BCRSMatrix<Dune::FieldMatrix<K, 1, 1>> CoarseMatrix;
  typedef Dune::BlockVector<Dune::FieldVector<int, 1>> vector1i;

  // define parameters
  const unsigned int dim = 2;
  const unsigned int degree = 1;
  const std::size_t nonzeros = std::pow(2*degree+1,dim);
  typedef double NumberType;

  vector1i local_basis_sizes = FromOffline<vector1i>(path_to_storage, "localBasisSizes");
  const int number_of_rank_used_offline = local_basis_sizes.size();
  std::string DefectID = "Defect1";
  vector1i changes_in_local_basis_sizes(local_basis_sizes.size());

  for( const auto Target : targeted) {

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
  //  Solving process begin here: First some parameters
    // ~~~~~~~~~~~~~~~~~~
    double eigenvalue_threshold = -1;
    // const int algebraic_overlap = 0;
    int nev = 6;
    int nev_arpack = nev;
    // double shift = 0.001;

    // ~~~~~~~~~~~~~~~~~~
    // Load the indices transformation
    // ~~~~~~~~~~~~~~~~~~
    vector1i offlineDoF2GI = FromOffline<vector1i>(path_to_storage, "GI", Target);
    vector1i DofOffline_to_DofOnline = offlineDoF2GI2gmsh2onlineDoF<vector1i>(Target, gmsh2dune, offlineDoF2GI, path_to_storage);

    // ~~~~~~~~~~~~~~~~~~
    // Load PoU
    // ~~~~~~~~~~~~~~~~~~
    Vector PoU = FromOffline<Vector>(path_to_storage, "PoU", Target);
    Vector nPoU(v_size);
    for (int i=0; i<v_size; i++){
      nPoU[DofOffline_to_DofOnline[i]] = PoU[i];
    }

    // ~~~~~~~~~~~~~~~~~~
  //  Subdomain basis computation and loading for neighbours
    // ~~~~~~~~~~~~~~~~~~
    std::shared_ptr<Dune::PDELab::SubdomainBasis<Vector>> online_subdomainbasis;
    online_subdomainbasis = std::make_shared<Dune::PDELab::GenEOBasisOnline<GO, Matrix, Vector>>(native(A), nPoU, eigenvalue_threshold, nev, nev_arpack);
    // ~~~~~~~~~~~~~~~~~~
  //  Particular solution
    // ~~~~~~~~~~~~~~~~~~
    // auto PartSol = ParticularSolution<Vector, Matrix>(native(A));
    // PartSol.exactRHS(native(fine_b));
    // PartSol.solveAndAppend(*online_subdomainbasis);

    // Change the targeted basis size for the database
    changes_in_local_basis_sizes[Target] = online_subdomainbasis->basis_size() - local_basis_sizes[Target];

    // ~~~~~~~~~~~~~~~~~~
    // Save local basis
    // ~~~~~~~~~~~~~~~~~~
    for (int basis_index = 0; basis_index < online_subdomainbasis->basis_size(); basis_index++) {
      std::string filename_EV = path_to_storage + std::to_string(Target) + "_" + DefectID + "EV_" + std::to_string(basis_index) + ".mm";
      Vector backToOfflineNumbering(online_subdomainbasis->get_basis_vector(basis_index)->size());
      for (int i=0; i<backToOfflineNumbering.size(); i++){
        backToOfflineNumbering[i] = (*online_subdomainbasis->get_basis_vector(basis_index))[DofOffline_to_DofOnline[i]];
      }
      Dune::storeMatrixMarket(backToOfflineNumbering, filename_EV, 16);
    }
  }

  std::string filename_ls = path_to_storage + DefectID + "_changeBasisSizes.mm";
  Dune::storeMatrixMarket(changes_in_local_basis_sizes, filename_ls);

}


int main(int argc, char **argv)
{
  using Dune::PDELab::Backend::native;

  // Offline folder
  std::string path_to_storage = "Offline/";

  // Define what subdomain need to be solved
  std::vector<int> targeted = {2}; // Subdomains that need a second solve

  // initialize MPI, finalize is done automatically on exit
  Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc,argv);

  driver(path_to_storage, targeted, helper);

  return 0;
}
