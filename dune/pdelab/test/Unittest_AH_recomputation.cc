#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/pdelab.hh>

#include <dune/pdelab/backend/istl/geneo/OfflineOnline/geneobasisOnline.hh>
#include <dune/pdelab/backend/istl/geneo/OfflineOnline/OnlineTools.hh>
#include <dune/pdelab/backend/istl/geneo/OfflineOnline/SubDomainGmshReader.hh>

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
    RF perm2 = 1e5; // FIXME we want high contrast
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

  // ~~~~~~~~~~~~~~~~~~
//  Solving process begin here: First some parameters
  // ~~~~~~~~~~~~~~~~~~
  double eigenvalue_threshold = -1;
  // const int algebraic_overlap = 0;
  int nev = 3;
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

  // ~~~~~~~~~~~~~~~~~~
  // Load the indices transformation
  // ~~~~~~~~~~~~~~~~~~
  Vector offlineDoF2GI;
  std::string filename_off2GI = path_to_storage + std::to_string(targeted[0]) + "_GI.mm";
  std::ifstream file_off2GI;
  file_off2GI.open(filename_off2GI.c_str(), std::ios::in);
  Dune::readMatrixMarket(offlineDoF2GI,file_off2GI);
  file_off2GI.close();

  std::vector<int> DofOffline_to_DofOnline = offlineDoF2GI2gmsh2onlineDoF<Vector>(targeted[0], gmsh2dune, offlineDoF2GI, path_to_storage);

  // ~~~~~~~~~~~~~~~~~~
  // Load Neighbour ranks
  // ~~~~~~~~~~~~~~~~~~
  Vector NR;
  std::string filename_NR = path_to_storage + std::to_string(targeted[0]) + "_neighborRanks.mm";
  std::ifstream file_NR;
  file_NR.open(filename_NR.c_str(), std::ios::in);
  Dune::readMatrixMarket(NR,file_NR);
  file_NR.close();

  // ~~~~~~~~~~~~~~~~~~
  // Load PoU
  // ~~~~~~~~~~~~~~~~~~
  Vector PoU;
  std::string filename_PoU = path_to_storage + std::to_string(targeted[0]) + "_PoU.mm";
  std::ifstream file_PoU;
  file_PoU.open(filename_PoU.c_str(), std::ios::in);
  Dune::readMatrixMarket(PoU,file_PoU);
  file_PoU.close();

  Vector nPoU(v_size);
  for (int i=0; i<v_size; i++){
    nPoU[DofOffline_to_DofOnline[i]] = PoU[i];
  }

  // ~~~~~~~~~~~~~~~~~~
//  Subdomain basis computation or loading
  // ~~~~~~~~~~~~~~~~~~

  int basis_size = local_basis_sizes[targeted[0]];

  // First : compute the new targeted subdomain basis
  std::shared_ptr<Dune::PDELab::SubdomainBasis<Vector>> online_subdomainbasis;
  online_subdomainbasis = std::make_shared<Dune::PDELab::GenEOBasisOnline<GO, Matrix, Vector>>(native(A), nPoU, eigenvalue_threshold, nev, nev_arpack);

  // Then : load other subdomain basis from offline and transfer them in the targeted subdomain space
  std::vector<std::shared_ptr<Dune::PDELab::SubdomainBasis<Vector>>> neighbour_subdomainbasis(NR.size());

  for (int iter_over_subdomains=0; iter_over_subdomains<NR.size(); iter_over_subdomains++) {

    Vector offlineNeighbourDoF2GI;
    int int_Nnumber = NR[iter_over_subdomains];
    std::string filename_ = path_to_storage + std::to_string(int_Nnumber) + "_GI.mm";
    std::ifstream file_;
    file_.open(filename_.c_str(), std::ios::in);
    Dune::readMatrixMarket(offlineNeighbourDoF2GI,file_);
    file_.close();

    neighbour_subdomainbasis[iter_over_subdomains] = std::make_shared<Dune::PDELab::NeighbourBasis<GO, Matrix, Vector>>(path_to_storage, local_basis_sizes[NR[iter_over_subdomains]], NR[iter_over_subdomains], offlineDoF2GI, offlineNeighbourDoF2GI, 2);

    for (int i=0; i<local_basis_sizes[iter_over_subdomains]; i++){
      auto tmp = *neighbour_subdomainbasis[iter_over_subdomains]->get_basis_vector(i);
      for (int j=0; j<v_size; j++){
        (*neighbour_subdomainbasis[iter_over_subdomains]->get_basis_vector(i))[DofOffline_to_DofOnline[j]] = tmp[j];
      }
    }
  }

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  Modify AH
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //  Load the coarse matrix
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  CoarseMatrix AH;
  std::ifstream file_AH;
  std::string filename_AH = path_to_storage + "OfflineAH.mm";
  file_AH.open(filename_AH.c_str(), std::ios::in);
  Dune::readMatrixMarket(AH,file_AH);
  file_AH.close();

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
      auto basis_vector = *online_subdomainbasis->get_basis_vector(basis_index_remote);
      Vector Atimesv(A.N());
      native(A).mv(basis_vector, Atimesv);
      for (int basis_index = 0; basis_index < local_basis_sizes[targeted[0]]; basis_index++) {
        Matrix::field_type entry = *online_subdomainbasis->get_basis_vector(basis_index)*Atimesv;
        local_rows[basis_index][NR.size()].push_back(entry);
      }
    }

    // Compute products of discretization matrix with local and remote vectors
    for (std::size_t neighbor_id = 0; neighbor_id < NR.size(); neighbor_id++) {
      if (basis_index_remote >= local_basis_sizes[NR[neighbor_id]])
        continue;
      auto basis_vector = *neighbour_subdomainbasis[neighbor_id]->get_basis_vector(basis_index_remote);
      Vector Atimesv(A.N());
      native(A).mv(basis_vector, Atimesv);
      for (int basis_index = 0; basis_index < local_basis_sizes[targeted[0]]; basis_index++) {
        Matrix::field_type entry = *online_subdomainbasis->get_basis_vector(basis_index)*Atimesv;
        local_rows[basis_index][neighbor_id].push_back(entry);
      }
    }
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
    for (std::size_t neighbor_id = 0; neighbor_id < NR.size(); neighbor_id++) {
      int neighbor_offset = local_offset[NR[neighbor_id]];
      for (int basis_index2 = 0; basis_index2 < local_basis_sizes[NR[neighbor_id]]; basis_index2++) {
        entries_pos[cnt] = neighbor_offset + basis_index2;
        cnt++;
      }
    }

    // Communicate actual entries
    Matrix::field_type entries[couplings];
    cnt = 0;
    for (int basis_index2 = 0; basis_index2 < local_basis_sizes[targeted[0]]; basis_index2++) {
      entries[cnt] = local_rows[basis_index][NR.size()][basis_index2];
      cnt++;
    }
    for (std::size_t neighbor_id = 0; neighbor_id < NR.size(); neighbor_id++) {
      for (int basis_index2 = 0; basis_index2 < local_basis_sizes[NR[neighbor_id]]; basis_index2++) {
        entries[cnt] = local_rows[basis_index][neighbor_id][basis_index2];
        cnt++;
      }
    }

    // Plot what was AH
    std::cout << "Offline AH["<< row_id <<"]\\/ " << std::endl;
    std::cout << "[" << AH[row_id][entries_pos[0]];
    for (int i = 1; i < couplings; i++){
      std::cout << "  " << AH[row_id][entries_pos[i]];
    }
    std::cout << "]" << std::endl;


    // Set matrix entries
    for (int i = 0; i < couplings; i++){
      // std::cout << "ici entries[i]:" << entries[i] << std::endl;
      // std::cout << "Before AH[" << row_id << "]["<< entries_pos[i] <<"]= " << AH[row_id][entries_pos[i]] << std::endl;
      AH[row_id][entries_pos[i]] = entries[i];
      // std::cout << "Then AH[" << row_id << "]["<< entries_pos[i] <<"]= " << AH[row_id][entries_pos[i]] << std::endl;
    }

    // Plot what was AH
    std::cout << "[" << AH[row_id][entries_pos[0]];
    for (int i = 1; i < couplings; i++){
      std::cout << "  " << AH[row_id][entries_pos[i]];
    }
    std::cout << "]" << std::endl;
    std::cout << "Online AH["<< row_id <<"]/\\ " << std::endl;
    std::cout << std::endl;

    row_id++;
  }

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
