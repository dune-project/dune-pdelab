#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/pdelab.hh>

#include <dune/pdelab/backend/istl/geneo/OfflineOnline/geneobasisOnline.hh>
#include <dune/pdelab/backend/istl/geneo/OfflineOnline/SubDomainGmshReader.hh>
#include <dune/pdelab/backend/istl/geneo/OfflineOnline/OnlineTools.hh>
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
  typedef GenericEllipticProblem<ES,NumberType> Problem;
  Problem problem;
  typedef Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<Problem> BCType;
  BCType bctype(es,problem);

  // ~~~~~~~~~~~~~~~~~~
//  FE fine Space definition
  // ~~~~~~~~~~~~~~~~~~
  typedef typename ES::Grid::ctype DF;
  // instantiate finite element maps
  typedef Dune::PDELab::QkLocalFiniteElementMap<ES,DF,NumberType,1> FEM;
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
  V fine_b(gfs,0.0);
  go.residual(x,fine_b);

  // std::cout << fine_b.N() << std::endl;

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
  // Load PoU
  // ~~~~~~~~~~~~~~~~~~
  Vector PoU = FromOffline<Vector>(path_to_storage, "PoU", targeted[0]);
  Vector nPoU(v_size);
  for (int i=0; i<v_size; i++){
    nPoU[DofOffline_to_DofOnline[i]] = PoU[i];
  }

  // ~~~~~~~~~~~~~~~~~~
  // Multiply fine b by PoU
  // ~~~~~~~~~~~~~~~~~~
  // for (int i=0; i<v_size; i++){
  //   native(fine_b)[i] *= nPoU[i];
  // }

  Vector OffFineb = FromOffline<Vector>(path_to_storage, "fineb", targeted[0]);
  assert(OffFineb.N() == v_size);
  Vector nfineb(v_size);
  for (int i=0; i<v_size; i++){
    nfineb[DofOffline_to_DofOnline[i]] = OffFineb[i];
  }

  Dune::VTKWriter<GV> vtkwriter(gv);
  V toplot(gfs,0.0);
  native(toplot) = nfineb;
  /* Write a field in the vtu file */
  DGF xdgf(gfs,toplot);
  auto adapt = std::make_shared<ADAPT>(xdgf,"offline");
  vtkwriter.addVertexData(adapt);
  DGF xdgf2(gfs,fine_b);
  auto adapt2 = std::make_shared<ADAPT>(xdgf2,"online");
  vtkwriter.addVertexData(adapt2);
  vtkwriter.write("RHS",Dune::VTK::ascii);

  // ~~~~~~~~~~~~~~~~~~
//  Subdomain basis computation and loading for neighbours
  // ~~~~~~~~~~~~~~~~~~

  std::shared_ptr<Dune::PDELab::SubdomainBasis<Vector>> online_subdomainbasis, offline_subdomainbasis;
  online_subdomainbasis = std::make_shared<Dune::PDELab::GenEOBasisOnline<GO, Matrix, Vector>>(native(A), nPoU, eigenvalue_threshold, nev, nev_arpack);

  // int basis_size = local_basis_sizes[targeted[0]];
  // offline_subdomainbasis = std::make_shared<Dune::PDELab::GenEOBasisFromFiles<GO, Matrix, Vector>>(path_to_storage, basis_size, targeted[0]);
  // for (int i=0; i<basis_size; i++){
  //   Vector tmp = *offline_subdomainbasis->get_basis_vector(i);
  //   for (int j=0; j<v_size; j++){
  //     (*offline_subdomainbasis->get_basis_vector(i))[DofOffline_to_DofOnline[j]] = tmp[j];
  //   }
  // }

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  Modify bH
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  CoarseVector bH = FromOffline<CoarseVector>(path_to_storage, "OfflineCoarseb");

  std::cout << "Offline bH: " << std::endl;
  std::cout << "[" << bH[0];
  for (int i = 1; i < bH.N(); i++){
    std::cout << "  " <<  bH[i];
  }
  std::cout << "]" << std::endl;

  UpdatebH<Vector, CoarseVector>(bH, native(fine_b), online_subdomainbasis, local_basis_sizes, local_offset, targeted[0]);
  // UpdatebH<Vector, CoarseVector>(bH, nfineb, offline_subdomainbasis, local_basis_sizes, local_offset, targeted[0]);

  std::cout << "Online bH: " << std::endl;
  std::cout << "[" << bH[0];
  for (int i = 1; i < bH.N(); i++){
    std::cout << "  " <<  bH[i];
  }
  std::cout << "]" << std::endl;

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
