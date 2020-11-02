#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/pdelab.hh>
#include <dune/pdelab/backend/istl/geneo/nonoverlapping/geneobasisOnline.hh>
#include <dune/pdelab/backend/istl/geneo/nonoverlapping/geneobasisfromfiles.hh>

#include <dune/pdelab/test/SubDomainGmshReader.hh>

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
  // set up and assemble right hand side w.r.t. l(v)-a(u_g,v)
  V d(gfs,0.0);
  go.residual(x,d); // The rhs is loaded from offline directly restricted in the coarse spaces

  // ~~~~~~~~~~~~~~~~~~
//  Solving process begin here: First some parameters
  // ~~~~~~~~~~~~~~~~~~
  double eigenvalue_threshold = -1;
  const int algebraic_overlap = 0;
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

  Vector online2GI(v_size);
  std::string filename_on2GI = path_to_storage + std::to_string(targeted[0]) + "_LocalToGlobalNode.mm";
  std::ifstream file_on2GI;
  file_on2GI.open(filename_on2GI.c_str(), std::ios::in);
  Dune::readMatrixMarket(online2GI,file_on2GI);
  file_on2GI.close();

  Vector offline2GI(v_size);
  std::string filename_off2GI = path_to_storage + std::to_string(targeted[0]) + "_GI.mm";
  std::ifstream file_off2GI;
  file_off2GI.open(filename_off2GI.c_str(), std::ios::in);
  Dune::readMatrixMarket(offline2GI,file_off2GI);
  file_off2GI.close();

  // /* First method to deal with online2GI and offline2GI */
  // std::vector<int> GI2gmsh;
  // auto result = std::max_element(online2GI.begin(), online2GI.end());
  // GI2gmsh.resize(online2GI[std::distance(online2GI.begin(), result)]);
  // for(int i=0; i<GI2gmsh.size(); i++)
  //   GI2gmsh[i]=-1;
  // for(int i=0; i<v_size; i++){
  //   GI2gmsh[online2GI[i]] = i;
  // }
  // /* End first method */

  /* Second method to deal with online2GI and offline2GI */
  std::vector<std::pair<int, int>> online(v_size), offline(v_size);
  for(int i=0; i<v_size; i++){
    online[i] = std::make_pair(online2GI[i],i);
    offline[i] = std::make_pair(offline2GI[i],i);
  }
  std::sort(online.begin(), online.end());
  std::sort(offline.begin(), offline.end());
  /* End second method */

  std::vector<int> DofOffline_to_DofOnline(v_size);
  // Change order from Dof offline to offline global
  for(int i=0; i<v_size; i++){ // go through dof offline ordering
    // DofOffline_to_DofOnline[i] = GI2gmsh[offline2GI[i]];
    DofOffline_to_DofOnline[offline[i].second] = online[i].second;
  }

  auto tmp = DofOffline_to_DofOnline;
  // We know that a re-ordering of nodes is done during the GmshReader operation
  for(int i=0; i<v_size; i++){ // go through dof offline ordering
    DofOffline_to_DofOnline[i] = gmsh2dune[tmp[i]];
  }

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
  // Load PoU
  // ~~~~~~~~~~~~~~~~~~
  // TODO : For now, only load the one associated to targeted[0] :: need to load more for neighbour
  Vector PoU(A.N());
  std::string filename_PoU = path_to_storage + std::to_string(targeted[0]) + "_PoU.mm";
  std::ifstream file_PoU;
  file_PoU.open(filename_PoU.c_str(), std::ios::in);
  Dune::readMatrixMarket(PoU,file_PoU);
  file_PoU.close();
  // std::shared_ptr<Vector> part_unity = std::make_shared<Vector>(PoU.N(), PoU);
  // std::shared_ptr<Matrix> A_ptr = std::make_shared<Matrix>(A);

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

  for (int i=0; i<NR.size(); i++) {
    subdomainbasis[iter_over_subdomains] = std::make_shared<Dune::PDELab::GenEOBasisFromFiles<GO, Matrix, Vector>>(path_to_storage, local_basis_sizes[iter_over_subdomains], iter_over_subdomains, 2);

    for (int i=0; i<local_basis_sizes[iter_over_subdomains]; i++){
      auto tmp = *subdomainbasis[iter_over_subdomains]->get_basis_vector(i);
      for (int j=0; j<v_size; j++){
        (*subdomainbasis[iter_over_subdomains]->get_basis_vector(i))[DofOffline_to_DofOnline[j]] = tmp[j];
      }
    }
  }




  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  Modify AH
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //  Load the coarse matrix
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // CoarseMatrix AH;
  // std::ifstream file_AH;
  // std::string filename_AH = path_to_storage + "OfflineAH.mm";
  // file_AH.open(filename_AH.c_str(), std::ios::in);
  // Dune::readMatrixMarket(AH,file_AH);
  // file_AH.close();
  // // std::cout << AH.N() << std::endl;
  // // std::cout << AH.M() << std::endl;}

  // int max_local_basis_size = *std::max_element(local_basis_sizes.begin(),local_basis_sizes.end());
  // int my_offset = local_offset[targeted[0]];

  // // ~~~~~~~~~~~~~~~~~~
  // //  Create a vector of AH entries modification
  // // ~~~~~~~~~~~~~~~~~~

  // // Set up container for storing rows of coarse matrix associated with current rank
  // std::vector<std::vector<std::vector<Matrix::field_type> > > local_rows;
  // local_rows.resize(local_basis_sizes[targeted[0]]);
  // for (int basis_index = 0; basis_index < local_basis_sizes[targeted[0]]; basis_index++) {
  //   local_rows[basis_index].resize(NR.size()+1);
  // }

  // for (int basis_index_remote = 0; basis_index_remote < max_local_basis_size; basis_index_remote++) {

  //   // Compute local products of basis functions with discretization matrix
  //   if (basis_index_remote < local_basis_sizes[targeted[0]]) {
  //     auto basis_vector = *subdomainbasis[targeted[0]]->get_basis_vector(basis_index_remote);
  //     Vector Atimesv(A.N());
  //     native(A).mv(basis_vector, Atimesv);
  //     for (int basis_index = 0; basis_index < local_basis_sizes[targeted[0]]; basis_index++) {
  //       Matrix::field_type entry = *subdomainbasis[targeted[0]]->get_basis_vector(basis_index)*Atimesv;
  //       local_rows[basis_index][NR.size()].push_back(entry);
  //     }
  //   }

  //   // Compute products of discretization matrix with local and remote vectors
  //   // for (std::size_t neighbor_id = 0; neighbor_id < NR.size(); neighbor_id++) {
  //   //   if (basis_index_remote >= local_basis_sizes[NR[neighbor_id]])
  //   //     continue;
  //   //   auto basis_vector = *subdomainbasis[NR[neighbor_id]]->get_basis_vector(basis_index_remote);
  //   //   Vector Atimesv(A.N());
  //   //   native(A).mv(basis_vector, Atimesv);
  //   //   for (int basis_index = 0; basis_index < local_basis_sizes[targeted[0]]; basis_index++) {
  //   //     Matrix::field_type entry = *subdomainbasis[targeted[0]]->get_basis_vector(basis_index)*Atimesv;
  //   //     local_rows[basis_index][neighbor_id].push_back(entry);
  //   //   }
  //   // }
  // }


  // // ~~~~~~~~~~~~~~~~~~
  // //  Modify AH entries
  // // ~~~~~~~~~~~~~~~~~~

  // int row_id = local_offset[targeted[0]];
  // // Modify AH entries with just computed local_rows
  // for (int basis_index = 0; basis_index < local_basis_sizes[targeted[0]]; basis_index++) {
  //   // Communicate number of entries in this row
  //   int couplings = local_basis_sizes[targeted[0]];
  //   for (int neighbor_id : NR) {
  //     couplings += local_basis_sizes[neighbor_id];
  //   }

  //   // Communicate row's pattern
  //   int entries_pos[couplings];
  //   int cnt = 0;
  //   for (int basis_index2 = 0; basis_index2 < local_basis_sizes[targeted[0]]; basis_index2++) {
  //     entries_pos[cnt] = my_offset + basis_index2;
  //     cnt++;
  //   }
  //   // for (std::size_t neighbor_id = 0; neighbor_id < NR.size(); neighbor_id++) {
  //   //   int neighbor_offset = local_offset[NR[neighbor_id]];
  //   //   for (int basis_index2 = 0; basis_index2 < local_basis_sizes[NR[neighbor_id]]; basis_index2++) {
  //   //     entries_pos[cnt] = neighbor_offset + basis_index2;
  //   //     cnt++;
  //   //   }
  //   // }

  //   // Communicate actual entries
  //   Matrix::field_type entries[couplings];
  //   cnt = 0;
  //   for (int basis_index2 = 0; basis_index2 < local_basis_sizes[targeted[0]]; basis_index2++) {
  //     entries[cnt] = local_rows[basis_index][NR.size()][basis_index2];
  //     cnt++;
  //   }
  //   // for (std::size_t neighbor_id = 0; neighbor_id < NR.size(); neighbor_id++) {
  //   //   for (int basis_index2 = 0; basis_index2 < local_basis_sizes[NR[neighbor_id]]; basis_index2++) {
  //   //     entries[cnt] = local_rows[basis_index][neighbor_id][basis_index2];
  //   //     cnt++;
  //   //   }
  //   // }


  //   // Set matrix entries
  //   for (int i = 0; i < couplings; i++){
  //     // std::cout << "ici entries[i]:" << entries[i] << std::endl;
  //     std::cout << "Before AH[row_id][entries_pos[i]]= " << AH[row_id][entries_pos[i]] << std::endl;
  //     AH[row_id][entries_pos[i]] = entries[i];
  //     std::cout << "Then AH[row_id][entries_pos[i]]= " << AH[row_id][entries_pos[i]] << std::endl;
  //   }

  //   row_id++;
  // }

}


int main(int argc, char **argv)
{
  using Dune::PDELab::Backend::native;

  // initialize MPI, finalize is done automatically on exit
  Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc,argv);

  driver("geneo", "standard", helper);

  return 0;
}
