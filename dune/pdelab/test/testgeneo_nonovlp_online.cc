#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/pdelab.hh>

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

template <typename V_>
struct AddGatherScatter
{
    static typename V_::value_type gather(const V_ &a, int i)
    {
        return a[i]; // I am sending my value
    }
    static void scatter(V_ &a, typename V_::value_type v, int i)
    {
        a[i] += v; // add what I receive to my value
    }
};

template<typename V_>
struct MultiGatherScatter
{
  static typename V_::value_type gather(const V_& a, int i)
  {
    return (*a.localVector_)[i];
  }
  static void scatterWithRank (V_& a, typename V_::value_type v, int i, int proc)
  {
    (*a.getVectorForRank(proc))[i]=v;
  }
};

template<typename GV, typename Vector, typename Matrix, typename rank_type>
class MultiVectorBundle { // TODO : recreate this function without using adapter which will disappear when we will solve online sequentially
public:
  typedef typename Vector::value_type value_type;

  MultiVectorBundle(Dune::NonoverlappingOverlapAdapter<GV, Vector, Matrix>& adapter, std::vector<rank_type> neighbor_ranks)
  : neighboringRanks_(neighbor_ranks),
    neighbor_basis(0)
  {
    for (auto& rank : neighboringRanks_) {
      neighbor_basis.push_back(std::make_shared<Vector>(adapter.getExtendedSize()));
    }
  }

  std::shared_ptr<Vector> getVectorForRank(int rank) {
    for (int i = 0; i < neighboringRanks_.size(); i++) {
      if (neighboringRanks_[i] == rank) {
        return neighbor_basis[i];
      }
    }
    DUNE_THROW(Dune::Exception, "Trying to access vector for unknown neighbor rank!");
    return nullptr;
  }

  void print() {
    for (size_t i = 0; i < neighboringRanks_.size(); i++) {
      Dune::printvector(std::cout, *(neighbor_basis[i]), "remote vec from neighbor " + std::to_string(neighboringRanks_[i]), "");
    }
  }

  std::shared_ptr<Vector> localVector_;
  std::vector<std::shared_ptr<Vector> > neighbor_basis;
  std::vector<int> neighboringRanks_;
};

template <class X>
class SubdomainBasisFromOffline {

  public:

  NonoverlappingGenEOBasisFromFiles(std::string& basename, int verbose = 0) {

    // std::ostringstream osrank;
    // osrank << adapter.gridView().comm().rank();
    // int basis_size;

    // // Get the basis size from file
    // std::string filename_basis_size = basename+ "_" + osrank.str() + "_size.txt";
    // std::ifstream input_basis_size;
    // input_basis_size.open(filename_basis_size, std::ios::in);
    // if (!input_basis_size.is_open())
    //   DUNE_THROW(IOError, "Could not open file: " << filename_basis_size);
    // input_basis_size >> basis_size;
    // input_basis_size.close();
    // this->local_basis.resize(basis_size);

    // for (int basis_index = 0; basis_index < this->local_basis.size(); basis_index++) {

    //   std::shared_ptr<X> ev = std::make_shared<X>();
    //   std::ostringstream rfilename;
    //   rfilename<< basename <<  "_" << basis_index  << "_" << adapter.gridView().comm().rank() << ".mm";
    //   std::ifstream file;
    //   file.open(rfilename.str().c_str(), std::ios::in);
    //   if(!file)
    //     DUNE_THROW(IOError, "Could not open file: " << rfilename.str().c_str());
    //   Dune::readMatrixMarket(*ev,file);
    //   file.close();

    //   this->local_basis[basis_index] = ev;
    // }
  }

};

void driver(std::string basis_type, std::string part_unity_type, Dune::MPIHelper& helper) {

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

  grid->loadBalance();


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
  // V d(gfs,0.0);
  // go.residual(x,d); // The rhs is loaded from offline directly restricted in the coarse spaces

  // ~~~~~~~~~~~~~~~~~~
//  Solving process begin here: First some parameters
  // ~~~~~~~~~~~~~~~~~~
  double eigenvalue_threshold = -1;
  const int algebraic_overlap = 1;
  int nev = 30;
  int nev_arpack = nev;
  // double shift = 0.001;

  std::vector<int> rank_to_be_solved = {0}; // Rank that need a second solve

  std::ifstream file;
  std::string path_to_storage = "Offline/";

  using Dune::PDELab::Backend::native;

  // ~~~~~~~~~~~~~~~~~~
//  ADAPTER :: TODO remove adapter which is created for a parallel solve
  // ~~~~~~~~~~~~~~~~~~
  Dune::NonoverlappingOverlapAdapter<GV, Vector, Matrix> adapter(gv, native(A), nonzeros, algebraic_overlap);
  using Attribute = Dune::EPISAttribute;
  Dune::AllSet<Attribute> allAttribute;
  auto allinterface = std::shared_ptr<Dune::Interface>(new Dune::Interface());
  allinterface->build(*adapter.getRemoteIndices(), allAttribute, allAttribute); // all to all communication
  // build up buffered communicator allowing communication over a dof vector
  auto communicator = std::shared_ptr<Dune::BufferedCommunicator>(new Dune::BufferedCommunicator());
  communicator->build<Vector>(*allinterface);

  auto communicatorWithRank = std::shared_ptr<DuneWithRank::BufferedCommunicator>(new DuneWithRank::BufferedCommunicator());
  communicatorWithRank->build<Vector>(*allinterface);

  // ~~~~~~~ Version from subdomainprojectedcoarsespace.hh ~~~~~~~
  // using Attribute = EPISAttribute;
  // Dune::AllSet<Attribute> allAttribute;
  // auto allinterface = std::shared_ptr<Dune::Interface>(new Dune::Interface());
  // allinterface->build(*adapter_.getRemoteIndices(),allAttribute,allAttribute); // all to all communication
  // auto communicator = std::shared_ptr<DuneWithRank::BufferedCommunicator>(new DuneWithRank::BufferedCommunicator());
  // communicator->build<X>(*allinterface);

  auto geneo_matrices = setupGenEOMatrices(go, adapter, A);
  std::shared_ptr<Matrix> A_extended = std::get<0>(geneo_matrices);
  std::shared_ptr<Matrix> A_overlap_extended = std::get<1>(geneo_matrices);
  std::shared_ptr<Vector> part_unity = std::get<2>(geneo_matrices);

  // ~~~~~~~~~~~~~~~~~~
//  Load a vector describing local basis sizes (number of EV) and creating the vector of offsets (to reach indices in the coarse space)
  // ~~~~~~~~~~~~~~~~~~
  std::string filename_lb = path_to_storage + "local_basis_sizes.txt";
  std::vector<int> local_basis_sizes(adapter.gridView().comm().size()), local_offset(adapter.gridView().comm().size()+1);
  file.open(filename_lb.c_str(), std::ios::in);
  int count=0;
  local_offset[0]=0;
  for (std::string line; std::getline(file, line); ) {
    int value = std::stoi(line);
    local_basis_sizes[count] = value;
    local_offset[count+1] = value+local_offset[count];
    if(adapter.gridView().comm().rank()==0)
      std::cout << "local_basis_sizes : " << local_basis_sizes[count] << ", local_offset : " << local_offset[count] << std::endl;
    count++;
  }
  file.close();

  // ~~~~~~~~~~~~~~~~~~
//  Subdomain basis computation or loading
  // ~~~~~~~~~~~~~~~~~~
  std::string basename = path_to_storage + "EV";
  std::shared_ptr<Dune::PDELab::SubdomainBasis<Vector>> subdomainbasis;
  // Load from offline the EV basis & recompute if needed
  std::vector<int>::iterator it = std::find(std::begin(rank_to_be_solved), std::end(rank_to_be_solved), adapter.gridView().comm().rank());
  if (it != rank_to_be_solved.end()) {
    subdomainbasis = std::make_shared<Dune::PDELab::NonoverlappingGenEOBasis<GO, Matrix, Vector>>(adapter, A_extended, A_overlap_extended, part_unity, eigenvalue_threshold, nev, nev_arpack);

    // Problem with the PARTICULAR SOLUTION
    // This new subdomain has no particular solution so the size of the new EV is reduced by 1
    // Several possibilities:
    // get only the particular solution from offline and add it to the subdomain basis:

    // recompute the particular solution for this subdomain:

    // not use any particular solution for this subdomain but the vector local_basis_sizes and local offset have to change:
    local_basis_sizes[adapter.gridView().comm().rank()] -= 1;
    for (int j=adapter.gridView().comm().rank()+1; j<local_basis_sizes.size(); j++) {
      local_offset[j]-=1;
    }

  } else {
    subdomainbasis = std::make_shared<Dune::PDELab::NonoverlappingGenEOBasisFromFiles<GV, Matrix, Vector>>(adapter, basename);
  }

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  Load the coarse matrix
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  CoarseMatrix AH;
  std::string filename = path_to_storage + "OfflineAH.mm";
  file.open(filename.c_str(), std::ios::in);
  Dune::readMatrixMarket(AH,file);
  file.close();
  // if(adapter.gridView().comm().rank()==0) {
  //   std::cout << AH.N() << std::endl;
  //   std::cout << AH.M() << std::endl;}

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  Modify AH
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  typedef int rank_type;
  rank_type affected_rank, ranks, my_rank, max_local_basis_size, my_offset, cnt;

  affected_rank = rank_to_be_solved[0]; // TODO for loop on rank_to_be_solved
  my_rank = adapter.gridView().comm().rank();

  std::vector<rank_type> neighbor_ranks = adapter.findNeighboringRanks(); // TODO replace when adapter will be removed :: from Offline

  ranks = local_basis_sizes.size();
  max_local_basis_size = *std::max_element(local_basis_sizes.begin(),local_basis_sizes.end());
  my_offset = local_offset[my_rank];

  // ~~~~~~~~~~~~~~~~~~
  //  Pre-step : create a vector with all subdomain basis
  // ~~~~~~~~~~~~~~~~~~
  // Objective access all subdomain basis from a only proc without using a communicator

  std::vector< std::vector< std::shared_ptr<Vector> > > all_subdomainbasis;
  all_subdomainbasis.resize(ranks);
  for (int rank_index=0; rank_index < all_subdomainbasis.size(); rank_index++) {
    // std::cout << "my_rank: " << my_rank << ", rank_index: " << rank_index << std::endl;
    if(my_rank == rank_index) {
      std::cout << "my_rank: " << my_rank << ", rank_index: " << rank_index << std::endl;
      std::cout << "ici all_subdomain filling" << std::endl;
      all_subdomainbasis[rank_index].resize(subdomainbasis->basis_size());
      std::cout << "all_subdomainbasis[rank_index].size(): " << all_subdomainbasis[rank_index].size() << std::endl;
      for (int basis_index=0; basis_index<all_subdomainbasis[rank_index].size(); basis_index++){
        all_subdomainbasis[rank_index][basis_index] = subdomainbasis->get_basis_vector(basis_index);
      }
    }
  }

  if (my_rank==0) {
    std::cout << "all_subdomainbasis[1].size() for rank 0: " << all_subdomainbasis[1].size() << std::endl;
  }

  // std::cout << "ici all_subdomain filling" << std::endl;
  // all_subdomainbasis[my_rank].resize(subdomainbasis->basis_size());
  // for (int basis_index=0; basis_index<all_subdomainbasis[my_rank].size(); basis_index++){
  //   all_subdomainbasis[my_rank][basis_index] = subdomainbasis->get_basis_vector(basis_index);
  // }

  // ~~~~~~~~~~~~~~~~~~
  //  Create a vector of AH entries modification
  // ~~~~~~~~~~~~~~~~~~

  if (my_rank == affected_rank) {
    std::cout << "Only the proc " << adapter.gridView().comm().rank() << " is used." << std::endl;

    std::cout << "local_basis_sizes.size(): " << local_basis_sizes.size() << std::endl;
    std::cout << "my_offset: " << my_offset << std::endl;

    std::cout << "Initializing local_rows: " << std::endl;

    // Set up container for storing rows of coarse matrix associated with current rank
    std::vector<std::vector<std::vector<Matrix::field_type> > > local_rows(local_basis_sizes[my_rank]);
    for (rank_type basis_index = 0; basis_index < local_basis_sizes[my_rank]; basis_index++) {
      local_rows[basis_index].resize(neighbor_ranks.size()+1);
    }

    std::cout << "... done " << std::endl;
    std::cout << "local_rows size: " << local_rows.size() << std::endl;
    std::cout << "local_rows[i] size: " << local_rows[my_rank].size() << std::endl;

    MultiVectorBundle<GV, Vector, Matrix, rank_type> bundle(adapter, neighbor_ranks);

    std::cout << "bundle is created" << std::endl;

    for (rank_type basis_index_remote = 0; basis_index_remote < max_local_basis_size; basis_index_remote++) {

      std::cout << basis_index_remote << std::endl;
      // Communicate one basis vectors of every subdomain to all of its neighbors in one go
      // If the current rank has already communicated all its basis vectors, just pass zeros
      // if (basis_index_remote < local_basis_sizes[my_rank]) { bundle.localVector_ = subdomainbasis->get_basis_vector(basis_index_remote); }
      // else { bundle.localVector_ = std::make_shared<Vector>(adapter.getExtendedSize()); }

      std::cout << "ici" << std::endl;

      // communicatorWithRank->forward<MultiGatherScatter<MultiVectorBundle<GV, Vector, Matrix, rank_type>>>(bundle,bundle); // make function known in other subdomains
      std::cout << "la" << std::endl;

      // Compute local products of basis functions with discretization matrix
      if (basis_index_remote < local_basis_sizes[my_rank]) {
        auto basis_vector = *subdomainbasis->get_basis_vector(basis_index_remote);
        Vector Atimesv(adapter.getExtendedSize());
        (*A_extended).mv(basis_vector, Atimesv);
        for (rank_type basis_index = 0; basis_index < local_basis_sizes[my_rank]; basis_index++) {
          Matrix::field_type entry = *subdomainbasis->get_basis_vector(basis_index)*Atimesv;
          local_rows[basis_index][neighbor_ranks.size()].push_back(entry);
        }
      }

      // Compute products of discretization matrix with local and remote vectors
      for (std::size_t neighbor_id = 0; neighbor_id < neighbor_ranks.size(); neighbor_id++) {
        if (basis_index_remote >= local_basis_sizes[neighbor_ranks[neighbor_id]])
          continue;
        std::cout << "neighbor_ranks[neighbor_id]: " << neighbor_ranks[neighbor_id] << std::endl;
        std::cout << "all_subdomainbasis.size(): " << all_subdomainbasis.size() << std::endl;
        std::cout << "all_subdomainbasis[neighbor_ranks[neighbor_id]].size(): " << all_subdomainbasis[neighbor_ranks[neighbor_id]].size() << std::endl;
        std::shared_ptr<Vector> basis_vector = all_subdomainbasis[neighbor_ranks[neighbor_id]][basis_index_remote]; // bundle.getVectorForRank(neighbor_ranks[neighbor_id]); //*bundle.neighbor_basis[neighbor_id];
        std::cout << "ici" << std::endl;
        Vector Atimesv(adapter.getExtendedSize());
        (*A_extended).mv(*basis_vector, Atimesv);
        for (rank_type basis_index = 0; basis_index < local_basis_sizes[my_rank]; basis_index++) {
          Matrix::field_type entry = *subdomainbasis->get_basis_vector(basis_index)*Atimesv;
          local_rows[basis_index][neighbor_id].push_back(entry);
        }
      }
    }

    std::cout << "Filling the vector local_rows of size " << local_rows.size() << " is done." << std::endl;


    // ~~~~~~~~~~~~~~~~~~
    //  Modify AH entries
    // ~~~~~~~~~~~~~~~~~~

    rank_type row_id = 0;
    // Modify AH entries with just computed local_rows
    for (rank_type basis_index = 0; basis_index < local_basis_sizes[my_rank]; basis_index++) {
      // Communicate number of entries in this row
      rank_type couplings = 0;
      couplings = local_basis_sizes[my_rank];
      for (rank_type neighbor_id : neighbor_ranks) {
        couplings += local_basis_sizes[neighbor_id];
      }
      adapter.gridView().comm().broadcast(&couplings, 1, my_rank);

      // Communicate row's pattern
      rank_type entries_pos[couplings];
      cnt = 0;
      for (rank_type basis_index2 = 0; basis_index2 < local_basis_sizes[my_rank]; basis_index2++) {
        entries_pos[cnt] = my_offset + basis_index2;
        cnt++;
      }
      for (std::size_t neighbor_id = 0; neighbor_id < neighbor_ranks.size(); neighbor_id++) {
        rank_type neighbor_offset = local_offset[neighbor_ranks[neighbor_id]];
        for (rank_type basis_index2 = 0; basis_index2 < local_basis_sizes[neighbor_ranks[neighbor_id]]; basis_index2++) {
          entries_pos[cnt] = neighbor_offset + basis_index2;
          cnt++;
        }
      }

      std::cout << "Filling the vector entries_pos is done." << std::endl;


      adapter.gridView().comm().broadcast(entries_pos, couplings, my_rank);

      // Communicate actual entries
      Matrix::field_type entries[couplings];
      cnt = 0;
      for (rank_type basis_index2 = 0; basis_index2 < local_basis_sizes[my_rank]; basis_index2++) {
        entries[cnt] = local_rows[basis_index][neighbor_ranks.size()][basis_index2];
        cnt++;
      }
      for (std::size_t neighbor_id = 0; neighbor_id < neighbor_ranks.size(); neighbor_id++) {
        for (rank_type basis_index2 = 0; basis_index2 < local_basis_sizes[neighbor_ranks[neighbor_id]]; basis_index2++) {
          entries[cnt] = local_rows[basis_index][neighbor_id][basis_index2];
          cnt++;
        }
      }
      adapter.gridView().comm().broadcast(entries, couplings, my_rank);

      std::cout << "Filling the vector entries is done." << std::endl;

      // Set matrix entries
      for (rank_type i = 0; i < couplings; i++) {
        AH[row_id][entries_pos[i]] = entries[i];
      }

      std::cout << "Filling AH of size " << AH.N() << "x" << AH.M() << " is done." << std::endl;
      std::cout << couplings << " entries have been changed." << std::endl;

      row_id++;
    }
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
  std::string filename_cb = path_to_storage + "OfflineCoarseb.mm";
  file.open(filename_cb.c_str(), std::ios::in);
  Dune::readMatrixMarket(coarse_d,file);
  file.close();

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
//  Prolongate the solution in order to have vsol on the fine space : vsol = RH^T * coarse_v = RH^T * AH^-1 * RH * b
  // ~~~~~~~~~~~~~~~~~~
  Vector v_fine_vovlp(adapter.getExtendedSize()); // TODO get rid of adapter in the sequential solve
  // Prolongate result
  v_fine_vovlp = 0.0;
  for (int basis_index = 0; basis_index < local_basis_sizes[gv.comm().rank()]; basis_index++) {
    Vector local_result(*subdomainbasis->get_basis_vector(basis_index));
    local_result *= coarse_v[local_offset[gv.comm().rank()] + basis_index];
    v_fine_vovlp += local_result;
  }

  communicator->forward<AddGatherScatter<Vector>>(v_fine_vovlp, v_fine_vovlp);

  V v(gfs, 0.0);
  adapter.restrictVector(v_fine_vovlp, v);
  native(x) -= v;

  // ~~~~~~~~~~~~~~~~~~
//  Write solution to VTK
  // ~~~~~~~~~~~~~~~~~~
  Dune::VTKWriter<GV> vtkwriter(gv);
  typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
  DGF xdgf(gfs,x);
  typedef Dune::PDELab::VTKGridFunctionAdapter<DGF> ADAPT;
  auto adapt = std::make_shared<ADAPT>(xdgf,"solution");
  vtkwriter.addVertexData(adapt);
  vtkwriter.write("nonovlptestgeneo_basis_" + basis_type + "_part_unity_" + part_unity_type);

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
