#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/pdelab.hh>
#include <dune/pdelab/backend/istl/geneo/nonoverlapping/geneobasisOnline.hh>
#include <dune/pdelab/backend/istl/geneo/nonoverlapping/geneobasisfromfiles.hh>

/*
 * Defining a Darcy problem with alternating layers of permeability and a high contrast
 */
template<typename GV, typename RF>
class GenericEllipticProblem
{
  typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;

public:
  typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

  //! tensor diffusion coefficient
  typename Traits::PermTensorType
  A (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::DomainType xglobal = e.geometry().global(x);

    RF perm1 = 1e0;
    RF perm2 = 1e0; // FIXME we want high contrast
    RF layer_thickness = 1.0 / 40.0;

    RF coeff = (int)std::floor(xglobal[0] / layer_thickness) % 2 == 0 ? perm1 : perm2;

    typename Traits::PermTensorType I;
    I[0][0] = coeff;
    I[0][1] = 0;
    I[1][0] = 0;
    I[1][1] = coeff;
    return I;
  }

  //! tensor diffusion constant per cell? return false if you want more than one evaluation of A per cell.
  static constexpr bool permeabilityIsConstantPerCell()
  {
    return true;
  }

  //! velocity field
  typename Traits::RangeType
  b (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::RangeType v(0.0);
    return v;
  }

  //! sink term
  typename Traits::RangeFieldType
  c (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return 0.0;
  }

  //! source term
  typename Traits::RangeFieldType
  f (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return 1.0;
  }

  BCType
  bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    typename Traits::DomainType xglobal = is.geometry().global(x);
    if (!(xglobal[0]<1E-6 || xglobal[0]>1.0-1E-6))
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
    else
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
  }

  //! Dirichlet boundary condition value
  typename Traits::RangeFieldType
  g (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::DomainType xglobal = e.geometry().global(x);
    /*if (xglobal[0] > 1.0-1E-6)
      return 1.0;
    else*/
      return 0.0;
  }

  //! flux boundary condition
  typename Traits::RangeFieldType
  j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return 0.0;
  }

  //! outflow boundary condition
  typename Traits::RangeFieldType
  o (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return 0.0;
  }
};

void driver(std::string basis_type, std::string part_unity_type, Dune::MPIHelper& helper) {

  // ~~~~~~~~~~~~~~~~~~
  // ReCreate offline hierarchy
  // ~~~~~~~~~~~~~~~~~~
  std::string path_to_storage = "Offline/";

  // ~~~~~~~~~~~~~~~~~~
  // Define what subdomain need to be solved
  // ~~~~~~~~~~~~~~~~~~
  std::vector<int> targeted = {0}; // Subdomains that need a second solve

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
  Dune::GmshReader<GRID>::read(factory, path_to_storage + std::to_string(targeted[0]) + "_subdomain.msh", true, true);
  std::unique_ptr<GRID> grid (factory.createGrid());

  typedef typename GRID::LeafGridView GV;
  auto gv = grid->leafGridView();

  // grid->loadBalance();


  // ~~~~~~~~~~~~~~~~~~
//  Type definitions
  // ~~~~~~~~~~~~~~~~~~
  const int components = 1;
  using K = double;
  using Vector = Dune::BlockVector<Dune::FieldVector<K,components>>;
  using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<K,components,components>>;

  typedef Dune::BlockVector<Dune::FieldVector<K, 1>> CoarseVector;
  typedef Dune::BCRSMatrix<Dune::FieldMatrix<K, 1, 1>> CoarseMatrix;

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
  go.residual(x,d); // The rhs is loaded from offline directly restricted in the coarse spaces

  // ~~~~~~~~~~~~~~~~~~
//  Solving process begin here: First some parameters
  // ~~~~~~~~~~~~~~~~~~
  double eigenvalue_threshold = -1;
  const int algebraic_overlap = 0;
  int nev = 30;
  int nev_arpack = nev;
  // double shift = 0.001;

  using Dune::PDELab::Backend::native;

//   // ~~~~~~~~~~~~~~~~~~
// //  ADAPTER :: TODO remove adapter which is created for a parallel solve
//   // ~~~~~~~~~~~~~~~~~~
//   Dune::NonoverlappingOverlapAdapter<GV, Vector, Matrix> adapter(gv, native(A), nonzeros, algebraic_overlap);
//   using Attribute = Dune::EPISAttribute;
//   Dune::AllSet<Attribute> allAttribute;
//   auto allinterface = std::shared_ptr<Dune::Interface>(new Dune::Interface());
//   allinterface->build(*adapter.getRemoteIndices(), allAttribute, allAttribute); // all to all communication
//   // build up buffered communicator allowing communication over a dof vector
//   auto communicator = std::shared_ptr<Dune::BufferedCommunicator>(new Dune::BufferedCommunicator());
//   communicator->build<Vector>(*allinterface);

//   auto communicatorWithRank = std::shared_ptr<DuneWithRank::BufferedCommunicator>(new DuneWithRank::BufferedCommunicator());
//   communicatorWithRank->build<Vector>(*allinterface);

//   // ~~~~~~~ Version from subdomainprojectedcoarsespace.hh ~~~~~~~
//   // using Attribute = EPISAttribute;
//   // Dune::AllSet<Attribute> allAttribute;
//   // auto allinterface = std::shared_ptr<Dune::Interface>(new Dune::Interface());
//   // allinterface->build(*adapter_.getRemoteIndices(),allAttribute,allAttribute); // all to all communication
//   // auto communicator = std::shared_ptr<DuneWithRank::BufferedCommunicator>(new DuneWithRank::BufferedCommunicator());
//   // communicator->build<X>(*allinterface);

//   auto geneo_matrices = setupGenEOMatrices(go, adapter, A);
//   std::shared_ptr<Matrix> A_extended = std::get<0>(geneo_matrices);
//   std::shared_ptr<Matrix> A_overlap_extended = std::get<1>(geneo_matrices);
//   std::shared_ptr<Vector> part_unity = std::get<2>(geneo_matrices);

  // ~~~~~~~~~~~~~~~~~~
//  Load a vector describing local basis sizes (number of EV) and creating the vector of offsets (to reach indices in the coarse space)
  // ~~~~~~~~~~~~~~~~~~
  Vector lb;
  std::ifstream file_lb;
  std::string filename_lb = path_to_storage + "localBasisSizes.mm";
  file_lb.open(filename_lb.c_str(), std::ios::in);
  Dune::readMatrixMarket(lb,file_lb);
  file_lb.close();


  const int number_of_rank_used_offline = lb.size();


  std::vector<int> local_basis_sizes(number_of_rank_used_offline), local_offset(number_of_rank_used_offline+1);
  local_offset[0]=0;
  for (int i=0; i<lb.size();i++) {
    local_basis_sizes[i] = lb[i];
    local_offset[i+1] = lb[i]+local_offset[i];
  }


  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //  Load the coarse matrix
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  CoarseMatrix AH;
  std::ifstream file_AH;
  std::string filename_AH = path_to_storage + "OfflineAH.mm";
  file_AH.open(filename_AH.c_str(), std::ios::in);
  Dune::readMatrixMarket(AH,file_AH);
  file_AH.close();
  // if(adapter.gridView().comm().rank()==0) {
  //   std::cout << AH.N() << std::endl;
  //   std::cout << AH.M() << std::endl;}


  // ~~~~~~~~~~~~~~~~~~
  // Load PoU
  // ~~~~~~~~~~~~~~~~~~
  // TODO : Only load the one associated to targeted[0] :: need to load more if needed
  Vector PoU(A.N());
  std::string filename_PoU = path_to_storage + std::to_string(targeted[0]) + "_PoU.mm";
  std::ifstream file_PoU;
  file_PoU.open(filename_PoU.c_str(), std::ios::in);
  Dune::readMatrixMarket(PoU,file_PoU);
  file_PoU.close();
  // std::shared_ptr<Vector> part_unity = std::make_shared<Vector>(PoU.N(), PoU);
  // std::shared_ptr<Matrix> A_ptr = std::make_shared<Matrix>(A);

  // ~~~~~~~~~~~~~~~~~~
  // Load Neighbour ranks
  // ~~~~~~~~~~~~~~~~~~
  // TODO : Only load the one associated to targeted[0] :: need to load more if needed
  Vector NR(A.N());
  std::string filename_NR = path_to_storage + std::to_string(targeted[0]) + "_neighborRanks.mm";
  std::ifstream file_NR;
  file_NR.open(filename_NR.c_str(), std::ios::in);
  Dune::readMatrixMarket(NR,file_NR);
  file_NR.close();


  // ~~~~~~~~~~~~~~~~~~
//  Subdomain basis computation or loading
  // ~~~~~~~~~~~~~~~~~~
  std::vector<std::shared_ptr<Dune::PDELab::SubdomainBasis<Vector>>> subdomainbasis(number_of_rank_used_offline);
  // Load from offline the EV basis & recompute if needed
  for (int iter_over_subdomains=0; iter_over_subdomains<number_of_rank_used_offline; iter_over_subdomains++) {
    std::vector<int>::iterator it = std::find(std::begin(targeted), std::end(targeted), iter_over_subdomains);

    if (it != targeted.end()) { // Recompute subdomain basis for the targeted subdomain
      subdomainbasis[iter_over_subdomains] = std::make_shared<Dune::PDELab::GenEOBasisOnline<GO, Matrix, Vector>>(native(A), PoU, eigenvalue_threshold, nev, nev_arpack);

      // TODO:: Problem with the PARTICULAR SOLUTION
      // This new subdomain has no particular solution so the size of the new EV is reduced by 1
      // Several possibilities:
      // get only the particular solution from offline and add it to the subdomain basis:
      // recompute the particular solution for this subdomain:
      // not use any particular solution for this subdomain but the vector local_basis_sizes and local offset have to change:

      // local_basis_sizes[targeted[0]] -= 1;
      // for (int j=targeted[0]+1; j<local_basis_sizes.size(); j++){
      //   local_offset[j]-=1;
      // }

    } else { // Load other subdomain basis from Offline
      // TODO :: only load neighbour subdomain basis, see exactly what we need after
      subdomainbasis[iter_over_subdomains] = std::make_shared<Dune::PDELab::GenEOBasisFromFiles<GO, Matrix, Vector>>(path_to_storage, local_basis_sizes[iter_over_subdomains], iter_over_subdomains, 2);
    }
  }


  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  Modify AH
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  int max_local_basis_size = *std::max_element(local_basis_sizes.begin(),local_basis_sizes.end());
  int my_offset = local_offset[targeted[0]];

  // ~~~~~~~~~~~~~~~~~~
  //  Create a vector of AH entries modification
  // ~~~~~~~~~~~~~~~~~~

  // Set up container for storing rows of coarse matrix associated with current rank
  std::vector<std::vector<std::vector<Matrix::field_type> > > local_rows;
  local_rows.resize(local_basis_sizes[targeted[0]]);
  for (int basis_index = 0; basis_index < local_basis_sizes[targeted[0]]; basis_index++) {
    local_rows[basis_index].resize(NR.size()+1);
  }

  for (int basis_index_remote = 0; basis_index_remote < max_local_basis_size; basis_index_remote++) {

    // Compute local products of basis functions with discretization matrix
    if (basis_index_remote < local_basis_sizes[targeted[0]]) {
      auto basis_vector = *subdomainbasis[targeted[0]]->get_basis_vector(basis_index_remote);
      Vector Atimesv(A.N());
      native(A).mv(basis_vector, Atimesv);
      for (int basis_index = 0; basis_index < local_basis_sizes[targeted[0]]; basis_index++) {
        Matrix::field_type entry = *subdomainbasis[targeted[0]]->get_basis_vector(basis_index)*Atimesv;
        local_rows[basis_index][NR.size()].push_back(entry);
      }
    }

    // Compute products of discretization matrix with local and remote vectors
    // for (std::size_t neighbor_id = 0; neighbor_id < NR.size(); neighbor_id++) {
    //   if (basis_index_remote >= local_basis_sizes[NR[neighbor_id]])
    //     continue;
    //   auto basis_vector = *subdomainbasis[NR[neighbor_id]]->get_basis_vector(basis_index_remote);
    //   Vector Atimesv(A.N());
    //   native(A).mv(basis_vector, Atimesv);
    //   for (int basis_index = 0; basis_index < local_basis_sizes[targeted[0]]; basis_index++) {
    //     Matrix::field_type entry = *subdomainbasis[targeted[0]]->get_basis_vector(basis_index)*Atimesv;
    //     local_rows[basis_index][neighbor_id].push_back(entry);
    //   }
    // }
  }

  // ~~~~~~~~~~~~~~~~~~
  //  Modify AH entries
  // ~~~~~~~~~~~~~~~~~~

  int row_id = local_offset[targeted[0]];
  // Modify AH entries with just computed local_rows
  for (int basis_index = 0; basis_index < local_basis_sizes[targeted[0]]; basis_index++) {
    // Communicate number of entries in this row
    int couplings = local_basis_sizes[targeted[0]];
    for (int neighbor_id : NR) {
      couplings += local_basis_sizes[neighbor_id];
    }

    // Communicate row's pattern
    int entries_pos[couplings];
    int cnt = 0;
    for (int basis_index2 = 0; basis_index2 < local_basis_sizes[targeted[0]]; basis_index2++) {
      entries_pos[cnt] = my_offset + basis_index2;
      cnt++;
    }
    // for (std::size_t neighbor_id = 0; neighbor_id < NR.size(); neighbor_id++) {
    //   int neighbor_offset = local_offset[NR[neighbor_id]];
    //   for (int basis_index2 = 0; basis_index2 < local_basis_sizes[NR[neighbor_id]]; basis_index2++) {
    //     entries_pos[cnt] = neighbor_offset + basis_index2;
    //     cnt++;
    //   }
    // }

    // Communicate actual entries
    Matrix::field_type entries[couplings];
    cnt = 0;
    for (int basis_index2 = 0; basis_index2 < local_basis_sizes[targeted[0]]; basis_index2++) {
      entries[cnt] = local_rows[basis_index][NR.size()][basis_index2];
      cnt++;
    }
    // for (std::size_t neighbor_id = 0; neighbor_id < NR.size(); neighbor_id++) {
    //   for (int basis_index2 = 0; basis_index2 < local_basis_sizes[NR[neighbor_id]]; basis_index2++) {
    //     entries[cnt] = local_rows[basis_index][neighbor_id][basis_index2];
    //     cnt++;
    //   }
    // }

    // Set matrix entries
    for (int i = 0; i < couplings; i++){
      // std::cout << "ici entries[i]:" << entries[i] << std::endl;
      // std::cout << "ici entries_pos[i]:" << entries_pos[i] << std::endl;
      AH[row_id][entries_pos[i]] = entries[i];
    }

    row_id++;
  }

  // ~~~~~~~~~~~~~~~~~~
//  Load the coarse right hand side
  // ~~~~~~~~~~~~~~~~~~
  // // Use the correct rhs to solve the final system
  // adapter.extendVector(native(d), b);
  // // Initializate a load vector (right hand side) in the coarse space : coarse_d = RH * b
  CoarseVector coarse_d;
  // coarse_space->restrict(b,coarse_d);
  // From Offline
  std::ifstream file_cb;
  std::string filename_cb = path_to_storage + "OfflineCoarseb.mm";
  file_cb.open(filename_cb.c_str(), std::ios::in);
  Dune::readMatrixMarket(coarse_d,file_cb);
  file_cb.close();

  // ~~~~~~~~~~~~~~~~~~
//  Solve the coarse space system
  // ~~~~~~~~~~~~~~~~~~
  // Objective :  find x = RH^T * AH^-1 * RH * b
  // Use of UMFPack to solve the problem [AH * coarse_v = RH * b] instead of inversing AH :: objective is to have [coarse_v = AH^-1 *  RH * b]
  Dune::UMFPack<CoarseMatrix> coarse_solver(AH, false);
  CoarseVector coarse_v(AH.N(), AH.M());
  Dune::InverseOperatorResult result;
  coarse_solver.apply(coarse_v,coarse_d,result);

  // ~~~~~~~~~~~~~~~~~~
//  Prolongate the solution and restrict it in order to have vsol on the nonoverlapping fine space : vsol = RH^T * coarse_v = RH^T * AH^-1 * RH * b
  // ~~~~~~~~~~~~~~~~~~
  Vector v_fine_ovlp(A.N());

  V v(gfs, 0.0); // TODO find the correct size of v because gfs changed an fill it properly

  // Prolongate result
  v_fine_ovlp = 0.0;
  for (int iter_over_subdomains=0; iter_over_subdomains<number_of_rank_used_offline; iter_over_subdomains++) {
    for (int basis_index = 0; basis_index < local_basis_sizes[iter_over_subdomains]; basis_index++) {
      Vector local_result(*subdomainbasis[iter_over_subdomains]->get_basis_vector(basis_index));
      local_result *= coarse_v[local_offset[iter_over_subdomains] + basis_index];
      v_fine_ovlp += local_result;
    }

    Vector Restrict(A.N()); // Overlapped subdomain to nonoverlapped subdomain
    std::string filename_restrict = path_to_storage + std::to_string(targeted[0]) + "_restrict.mm";
    std::ifstream file_restrict;
    file_restrict.open(filename_restrict.c_str(), std::ios::in);
    Dune::readMatrixMarket(Restrict,file_restrict);
    file_restrict.close();

    for (int i=0; i<Restrict.size(); i++)
      std::cout << Restrict[i] << std::endl;

  }




//   // ~~~~~~~~~~~~~~~~~~
// //  Write solution to VTK
//   // ~~~~~~~~~~~~~~~~~~
//   Dune::VTKWriter<GV> vtkwriter(gv);
//   typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
//   DGF xdgf(gfs,x);
//   typedef Dune::PDELab::VTKGridFunctionAdapter<DGF> ADAPT;
//   auto adapt = std::make_shared<ADAPT>(xdgf,"solution");
//   vtkwriter.addVertexData(adapt);
//   vtkwriter.write("nonovlptestgeneo_basis_" + basis_type + "_part_unity_" + part_unity_type);

  // ~~~~~~~~~~~~~~~~~~
  // Visualise all the basis in vtk format
  // ~~~~~~~~~~~~~~~~~~
  // for (int basis_index = 0; basis_index < subdomainbasis->basis_size(); basis_index++) {
  //     Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,Dune::refinementLevels(0));
  //     V vect(gfs, 0.0);
  //     adapter.restrictVector(native(*subdomainbasis->get_basis_vector(basis_index)), native(vect));

  //     int rank = adapter.gridView().comm().rank();
  //     std::string filename = "BasisVector_"+std::to_string(basis_index);

  //     Dune::PDELab::vtk::DefaultFunctionNameGenerator fieldname;
  //     fieldname.prefix("EV");

  //     Dune::PDELab::addSolutionToVTKWriter(vtkwriter,gfs,vect,fieldname);
  //     vtkwriter.write(filename,Dune::VTK::ascii);}
}


int main(int argc, char **argv)
{
  using Dune::PDELab::Backend::native;

  // initialize MPI, finalize is done automatically on exit
  Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc,argv);

  driver("geneo", "standard", helper);

  return 0;
}
