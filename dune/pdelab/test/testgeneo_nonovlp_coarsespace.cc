#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/pdelab.hh>

#include <dune/grid/utility/parmetisgridpartitioner.hh>

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


void driver(std::string basis_type, std::string part_unity_type, Dune::MPIHelper& helper) {


  unsigned int cells = 24;
  int overlap = 0;

  // define parameters
  const unsigned int dim = 2;
  const unsigned int degree = 1;
  const std::size_t nonzeros = std::pow(2*degree+1,dim);
  typedef double NumberType;

  // build a grid
  typedef Dune::UGGrid<dim> GM;
  Dune::FieldVector<double,dim> l(0.0);
  Dune::FieldVector<double,dim> u(1.0);
  std::array<unsigned int,dim> N;
  std::shared_ptr<GM> grid;
  N.fill(cells);
  grid = Dune::StructuredGridFactory<GM>::createCubeGrid(l,u,N);


  typedef typename GM::LeafGridView GV;
  auto gv = grid->leafGridView();

  // Transfer partitioning from ParMETIS to our grid
#if PARMETIS_MAJOR_VERSION
  std::vector<unsigned> part(Dune::ParMetisGridPartitioner<GV>::partition(gv, helper));
  grid->loadBalance(part, 0);
#else
  grid->loadBalance();
#endif

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


  // make a finite element space

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

  double eigenvalue_threshold = -1;
  const int algebraic_overlap = 1;
  int nev = 10;
  int nev_arpack = nev;
  // double shift = 0.001;

  int multiscale = 0;
  std::vector<int> proc_to_be_solved = {1, 2};

  using Dune::PDELab::Backend::native;

  // const GV &gv = go.trialGridFunctionSpace().gridView();

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


  std::string path_to_storage = "Offline/";
  std::string basename = path_to_storage + "EV";

  std::shared_ptr<Dune::PDELab::SubdomainBasis<Vector>> subdomainbasis;

  if (multiscale==1) { // first step of the multiscale FRAMEWORK: solving the pristine model & saving it to file
    subdomainbasis = std::make_shared<Dune::PDELab::NonoverlappingGenEOBasis<GO, Matrix, Vector>>(adapter, A_extended, A_overlap_extended, part_unity, eigenvalue_threshold, nev, nev_arpack);
    // Save the EV basis
    int rank = adapter.gridView().comm().rank();
    subdomainbasis->to_file(basename, rank);
  } else if (multiscale==2) { // other step of the multiscale FRAMEWORK: loading the subdomain basis from files
    std::vector<int>::iterator it = std::find(std::begin(proc_to_be_solved), std::end(proc_to_be_solved), adapter.gridView().comm().rank());
    if (it != proc_to_be_solved.end()) {
      subdomainbasis = std::make_shared<Dune::PDELab::NonoverlappingGenEOBasis<GO, Matrix, Vector>>(adapter, A_extended, A_overlap_extended, part_unity, eigenvalue_threshold, nev, nev_arpack);
    } else {
      subdomainbasis = std::make_shared<Dune::PDELab::NonoverlappingGenEOBasisFromFiles<GV, Matrix, Vector>>(adapter, basename);
    }
  } else if (multiscale==3) {
    // Test case for write/read database
    // Check numbers of digit initially and after the saving/reading procedure
    subdomainbasis = std::make_shared<Dune::PDELab::NonoverlappingGenEOBasis<GO, Matrix, Vector>>(adapter, A_extended, A_overlap_extended, part_unity, eigenvalue_threshold, nev, nev_arpack);

    int rank = adapter.gridView().comm().rank();
    subdomainbasis->to_file(basename, rank);

    auto fromfile_subdomainbasis = std::make_shared<Dune::PDELab::NonoverlappingGenEOBasisFromFiles<GV, Matrix, Vector>>(adapter, basename);
    fromfile_subdomainbasis->to_file(basename+"_rewritten", rank);
    // auto tmp = (*subdomainbasis->get_basis_vector(0));
    // tmp -= (*fromfile_subdomainbasis->get_basis_vector(0));
    // std::cout << tmp.two_norm() << std::endl;
  } else { // Classic case: no need to use multiscale FRAMEWORK
    subdomainbasis = std::make_shared<Dune::PDELab::NonoverlappingGenEOBasis<GO, Matrix, Vector>>(adapter, A_extended, A_overlap_extended, part_unity, eigenvalue_threshold, nev, nev_arpack);
  }

  // Extend the load vector d (right hand side) of the fine to space to its virtual overlap
  Vector b(adapter.getExtendedSize());

  // if(adapter.gridView().comm().rank()==0) std::cout<< b.size() << std::endl;
  //~ Add a solution to another rhs here fill with ones
  for(auto it = b.begin(); it!=b.end(); ++it){
    b[it.index()] += 1.0;
    // std::srand(std::time(0));
    // b[it.index()] = -1.0 + 2.0* (std::rand()+0.0) / (RAND_MAX + 1.0);
    // if(adapter.gridView().comm().rank()==0) std::cout<< b[it.index()] << std::endl;
  }

  //~ Add the true solution
  // adapter.extendVector(native(d), b);

  communicator->forward<AddGatherScatter<Vector>>(b,b); // make function known in other subdomains

  const int block_size = Vector::block_type::dimension;
  Matrix A_dirichlet = *A_extended;
  // Apply Dirichlet conditions to matrix on processor boundaries, inferred from partition of unity
  for (auto rIt=A_dirichlet.begin(); rIt!=A_dirichlet.end(); ++rIt){
      for(int block_i = 0; block_i < block_size; block_i++){
          if ((*part_unity)[rIt.index()][block_i] == .0){
              for (auto cIt=rIt->begin(); cIt!=rIt->end(); ++cIt){
                  for(int block_j = 0; block_j < block_size; block_j++){
                      (*cIt)[block_i][block_j] = (rIt.index() == cIt.index() && block_i == block_j) ? 1.0 : 0.0;
                  }
              }
          }
      }
  }

  Vector b_cpy(b);
  for (auto rIt=A_dirichlet.begin(); rIt!=A_dirichlet.end(); ++rIt) {
      for(int block_i = 0; block_i < block_size; block_i++){
          bool isDirichlet = true;
          for (auto cIt=rIt->begin(); cIt!=rIt->end(); ++cIt){
              for(int block_j = 0; block_j < block_size; block_j++){
                  if ((rIt.index() != cIt.index() || block_i!=block_j) && (*cIt)[block_i][block_j] != 0.0){
                      isDirichlet = false;
                      break;
                  }
              }
              if(!isDirichlet) break;
          }
          if (isDirichlet){
              b_cpy[rIt.index()] = .0;
              b[rIt.index()] = .0;
          }
      }
  }

  // Compute the fine solution for all subdomains
  // std::shared_ptr<Vector> ui(adapter.getExtendedSize());
  Vector ui(adapter.getExtendedSize());
  Dune::UMFPack<Matrix> subdomain_solver(A_dirichlet, false);
  Dune::InverseOperatorResult result1;
  subdomain_solver.apply(ui,b_cpy,result1);

  subdomainbasis->append(ui);

  // Visualise all the basis in vtk format
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

  auto coarse_space = std::make_shared<Dune::PDELab::NonoverlappingSubdomainProjectedCoarseSpace<GV, Matrix, Vector>>(adapter, gv, *A_extended, subdomainbasis, verbose);

  // ~~~~~~~~~~ Solve the coarse space system ~~~~~~~~~~
  // Objective :  find x = RH^T * AH^-1 * RH * b

  // Use the correct rhs to solve the final system
  adapter.extendVector(native(d), b);

  // Initializate a load vector (right hand side) in the coarse space : coarse_d = RH * b
  CoarseVector coarse_d(coarse_space->basis_size(), coarse_space->basis_size());
  coarse_space->restrict(b,coarse_d);

  // Get the coarse matrix from the coarse space
  std::shared_ptr<CoarseMatrix> AH = coarse_space->get_coarse_system();


  // Use of UMFPack to solve the problem [AH * coarse_v = RH * b] instead of inversing AH :: objective is to have [coarse_v = AH^-1 *  RH * b]
  Dune::UMFPack<CoarseMatrix> coarse_solver(*AH, false);
  CoarseVector coarse_v(coarse_space->basis_size(), coarse_space->basis_size());
  Dune::InverseOperatorResult result;
  coarse_solver.apply(coarse_v,coarse_d,result);

  // Prolongate the solution in order to have vsol on the fine space : vsol = RH^T * coarse_v = RH^T * AH^-1 * RH * b
  Vector vsol(adapter.getExtendedSize());
  coarse_space->prolongate(coarse_v, vsol);

  communicator->forward<AddGatherScatter<Vector>>(vsol, vsol);

  V v(gfs, 0.0);
  adapter.restrictVector(vsol, v);

  native(x) -= v;

  // Write solution to VTK
  Dune::VTKWriter<GV> vtkwriter(gv);
  typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
  DGF xdgf(gfs,x);
  typedef Dune::PDELab::VTKGridFunctionAdapter<DGF> ADAPT;
  auto adapt = std::make_shared<ADAPT>(xdgf,"solution");
  vtkwriter.addVertexData(adapt);
  vtkwriter.write("nonovlptestgeneo_basis_" + basis_type + "_part_unity_" + part_unity_type);
}


int main(int argc, char **argv)
{
  using Dune::PDELab::Backend::native;

  // initialize MPI, finalize is done automatically on exit
  Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc,argv);

  driver("geneo", "standard", helper);

  return 0;
}
