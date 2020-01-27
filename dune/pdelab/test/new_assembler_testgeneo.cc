#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/pdelab.hh>

#include <dune/grid/utility/parmetisgridpartitioner.hh>

#include <dune/istl/matrixmarket.hh>

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
      return 0.0; // TODO: support Dirichlet again
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

class GhostExcluder { // Who ya gonna call? Ghost excluders!!
public:

  template<typename LFSCache, typename Entity>
  bool assembleCell(const Entity& entity, const LFSCache & cache) const {
    return entity.partitionType() != Dune::PartitionType::GhostEntity;
  }
private:
};


template <typename Vector>
class InteriorExcluder {
public:

  InteriorExcluder(const Vector& partUnity)
   : partUnity_(partUnity) {
  }

  template<typename LFSCache, typename Entity>
  bool assembleCell(const Entity& entity, const LFSCache & cache) const {
    if (entity.partitionType() == Dune::PartitionType::GhostEntity)
      return false;

    for (std::size_t i = 0; i < cache.size(); i++)
    {
      if (partUnity_[cache.containerIndex(i)[0]] > 0.0 &&
        partUnity_[cache.containerIndex(i)[0]] < 1.0)
        return true;
    }
    return false;
  }
private:
  const Vector& partUnity_;
};


void driver(std::string basis_type, std::string part_unity_type, Dune::MPIHelper& helper) {

  Dune::Timer timer_detailed(true);

  unsigned int cells = 24;
  //int cells = 80;
  int overlap = 0;

  // define parameters
  //const unsigned int dim = 1;
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

  std::cout << "grid setup: " << timer_detailed.elapsed() << std::endl; timer_detailed.reset();

  typedef typename GM::LeafGridView GV;
  auto gv = grid->leafGridView();

  // Transfer partitioning from ParMETIS to our grid
  std::vector<unsigned> part(Dune::ParMetisGridPartitioner<GV>::partition(gv, helper));
  grid->loadBalance(part, 0);


  /*typedef Dune::YaspGrid<dim> GM;


  Dune::FieldVector<double,dim> L(1.0);
  std::array<int,dim> N;
  N.fill(cells);
  std::bitset<dim> B(false);

  auto grid = std::make_shared<GM>(L,N,B,overlap,Dune::MPIHelper::getCollectiveCommunication());
  typedef typename GM::LeafGridView GV_;
  auto gv_ = grid->leafGridView();

  std::cout << "grid setup: " << timer_detailed.elapsed() << std::endl; timer_detailed.reset();*/




  const int components = 1;
  using K = double;
  using Vector = Dune::BlockVector<Dune::FieldVector<K,components>>;
  using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<K,components,components>>;

  // make problem parameters
  typedef GenericEllipticProblem<GV,NumberType> Problem;
  Problem problem;
  typedef Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<Problem> BCType;
  BCType bctype(gv,problem);


  // make a finite element space

  typedef typename GV::Grid::ctype DF;
  // instantiate finite element maps
  typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,double,1> FEM;
  FEM fem(gv);

  // function space with no constraints on processor boundaries, needed for the GenEO eigenproblem
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,
                                          Dune::PDELab::ConformingDirichletConstraints,
                                          Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed,1>
                                          > GFS;
  GFS gfs(gv,fem);

  int verb=0;
  if (gfs.gridView().comm().rank()==0) verb=2;

  // make a degree of freedom vector on fine grid and initialize it with interpolation of Dirichlet condition
  typedef Dune::PDELab::Backend::Vector<GFS,NumberType> V;
  V x(gfs,0.0);

  // Extract domain boundary constraints from problem definition, apply trace to solution vector
  typedef Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<Problem> G;
  G g(gv,problem);
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

  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,NumberType,NumberType,NumberType,CC,CC,GhostExcluder> GO;
  GhostExcluder asmb_ghost_excluder;
  auto go = GO(gfs,cc,gfs,cc,lop,MBE(nonzeros),asmb_ghost_excluder);


  std::cout << "problem definition: " << timer_detailed.elapsed() << std::endl; timer_detailed.reset();


  // set up and assemble right hand side w.r.t. l(v)-a(u_g,v)
  V d(gfs,0.0);
  go.residual(x,d); // NOTE: Need operator without Dirichlet on proc boundaries for rhs setup!

  // types
  typedef GO::Jacobian M;


  // Assemble fine grid matrix defined without processor constraints
  M A(go);
  go.jacobian(x,A);

  std::cout << "fine assembly: " << timer_detailed.elapsed() << std::endl; timer_detailed.reset();

  using Dune::PDELab::Backend::native;
  int avg = nonzeros;
  int algebraic_overlap = 1;


  Dune::NonoverlappingOverlapAdapter<GV, Vector, Matrix> adapter(gv, native(A), avg, algebraic_overlap);

  std::shared_ptr<Matrix> A_extended = adapter.extendMatrix(native(A));
  std::shared_ptr<Vector> part_unity = Dune::makePartitionOfUnity(adapter, *A_extended);


  using Attribute = Dune::EPISAttribute;
  Dune::AllSet<Attribute> allAttribute;
  auto allinterface = std::shared_ptr<Dune::Interface>(new Dune::Interface());
  allinterface->build(*adapter.getRemoteIndices(),allAttribute,allAttribute); // all to all communication
  auto communicator = std::shared_ptr<Dune::BufferedCommunicator>(new Dune::BufferedCommunicator());
  communicator->build<Vector>(*allinterface);

  Dune::PDELab::MultiVectorBundle<GV, Vector, Matrix> remotePartUnities(adapter);
  remotePartUnities.localVector_ = part_unity;
  communicator->forward<Dune::PDELab::MultiGatherScatter<Dune::PDELab::MultiVectorBundle<GV, Vector, Matrix>>>(remotePartUnities,remotePartUnities); // make function known in other subdomains

  std::cout << "comm part_unity: " << timer_detailed.elapsed() << std::endl; timer_detailed.reset();


  // Assemble fine grid matrix defined only on overlap region
  auto part_unity_restricted = std::make_shared<Vector>(native(A).N());
  adapter.restrictVector(*part_unity, *part_unity_restricted);


  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,NumberType,NumberType,NumberType,CC,CC,InteriorExcluder<Vector>> GO_OVLP;
  InteriorExcluder<Vector> excluder(*part_unity_restricted);
  auto go_overlap = GO_OVLP(gfs,cc,gfs,cc,lop,MBE(nonzeros),excluder);


  M A_ovlp(go_overlap);
  go_overlap.jacobian(x,A_ovlp);

  std::cout << "local A_ovlp: " << timer_detailed.elapsed() << std::endl; timer_detailed.reset();


  // Choose an eigenvalue threshold according to Spillane et al., 2014.
  // This particular choice is a heuristic working very well for Darcy problems.
  // Theoretically, any value delivers robustness; practically, you may need to
  // choose another value to achieve a good balance between condition bound and
  // global basis size
  double eigenvalue_threshold = (double)overlap / (cells + overlap);


  int verbose = verb;
  using ScalarVector = Dune::BlockVector<Dune::FieldVector<K,1>>;



  M newmat(go);
  auto extended_matrices = adapter.lambdaMultiExtendMatrix(native(A_ovlp), native(A), [&](int i){
std::cout << "Lambda called for " << i << std::endl;
    std::shared_ptr<Vector> neighbor_part_unity = remotePartUnities.getVectorForRank(i); //.neighbor_basis[i];

    adapter.restrictVector(*neighbor_part_unity, *part_unity_restricted);

    //std::shared_ptr<EntitySetExcluder<Vector, GV_>> es_local_pou_excluder = std::make_shared<EntitySetPartUnityExcluder<Vector, GV_, GFS>> (gfs, part_unity_restricted);
    //gv.setExcluder(es_local_pou_excluder);
    //gv.setExcluder(ghost_excluder); // FIXME

    InteriorExcluder<Vector> excluder(*part_unity_restricted);
    typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,NumberType,NumberType,NumberType,CC,CC,InteriorExcluder<Vector>> GO_OVLP;
    auto go_overlap = GO_OVLP(gfs,cc,gfs,cc,lop,MBE(nonzeros),excluder);
    //typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,NumberType,NumberType,NumberType,CC,CC> GO_OVLP;
    //auto go_overlap = GO_OVLP(gfs,cc,gfs,cc,lop,MBE(nonzeros));



    for (auto rIt=native(newmat).begin(); rIt!=native(newmat).end(); ++rIt)
      for (auto cIt=rIt->begin(); cIt!=rIt->end(); ++cIt) {
        *cIt = 0.0;
      }

std::cout << "Jacobian for " << i << std::endl;
    go_overlap.jacobian(x,newmat);
    //gv.setExcluder(ghost_excluder);

    if (gfs.gridView().comm().rank()==0) {
      using Dune::PDELab::Backend::native;
      Dune::storeMatrixMarket(native(newmat), "comm_mat" + std::to_string(i) + ".mm");
    }

std::cout << "Return for " << i << std::endl;
    return stackobject_to_shared_ptr(native(newmat));
  });

  std::shared_ptr<Matrix> A_ovlp_extended = extended_matrices.first;
  A_extended = extended_matrices.second;


  std::cout << "Comm ovlp matrices: " << timer_detailed.elapsed() << std::endl; timer_detailed.reset();
  std::cout << "Basis setup" << std::endl;
  // Enforce problem's Dirichlet condition on PoU
  for (auto rIt=A_extended->begin(); rIt!=A_extended->end(); ++rIt) {
    bool isDirichlet = true;
    for (auto cIt=rIt->begin(); cIt!=rIt->end(); ++cIt)
    {
      if (rIt.index() != cIt.index() && *cIt != 0.0) {
        isDirichlet = false;
        break;
      }
    }
    if (isDirichlet) {
      (*part_unity)[rIt.index()] = .0;
    }
  }

  std::cout << "part_unity with Dirichlet: " << timer_detailed.elapsed() << std::endl; timer_detailed.reset();

  int nev = 2;
  auto subdomainbasis = std::make_shared<Dune::PDELab::NewGenEOBasis<GV, Matrix, Vector>>(adapter, *A_extended, *A_ovlp_extended, *part_unity, -1.0, nev);

  std::cout << "eigenproblems: " << timer_detailed.elapsed() << std::endl; timer_detailed.reset();


  auto coarse_space = std::make_shared<Dune::PDELab::NewSubdomainProjectedCoarseSpace<GV, Matrix, Vector>>(adapter, gv, *A_extended, subdomainbasis, verbose);

  // Apply Dirichlet conditions on processor boundaries, needed for Schwarz method
  for (auto rIt=A_extended->begin(); rIt!=A_extended->end(); ++rIt)
    for (auto cIt=rIt->begin(); cIt!=rIt->end(); ++cIt)
    {
      if ((*part_unity)[rIt.index()] == .0) {
        *cIt = rIt.index() == cIt.index() ? 1.0 : 0.0;
      }
    }

  Dune::PDELab::ISTL::NewTwoLevelOverlappingAdditiveSchwarz<GV, Matrix, Matrix, Vector,Vector> prec(adapter, *A_extended, coarse_space, true, verb);


  NonoverlappingOperator<GV, Matrix,Vector> linearOperator(gv,native(A));
  NonoverlappingScalarProduct<GV,Vector> scalarproduct(gv,native(x));
  Dune::CGSolver<Vector> solver(linearOperator,scalarproduct,prec,1e-6,500,verbose);

  Vector v(native(x));
  Dune::InverseOperatorResult stat;
  solver.apply(native(v),native(d),stat);
  native(x) -= v;


  // Write solution to VTK
  Dune::VTKWriter<GV> vtkwriter(gv);//gfs.gridView());
  typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
  DGF xdgf(gfs,x);
  typedef Dune::PDELab::VTKGridFunctionAdapter<DGF> ADAPT;
  auto adapt = std::make_shared<ADAPT>(xdgf,"solution");
  vtkwriter.addVertexData(adapt);
  vtkwriter.write("newassemblertestgeneo_basis_" + basis_type + "_part_unity_" + part_unity_type);
}


int main(int argc, char **argv)
{
  using Dune::PDELab::Backend::native;

  try{
    // initialize MPI, finalize is done automatically on exit
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc,argv);

    driver("geneo", "standard", helper);

    return 0;
  }

  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
    return 1;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
    return 1;
  }

}
