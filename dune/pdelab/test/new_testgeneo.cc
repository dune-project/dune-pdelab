#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/pdelab/boilerplate/pdelab.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionadapter.hh>
#include <dune/pdelab/localoperator/convectiondiffusionfem.hh>
#include <dune/pdelab/gridfunctionspace/vtk.hh>

#include <dune/pdelab/backend/istl/geneo/geneo.hh>

#include <dune/pdelab/backend/istl/geneo/novlp_operators.hh>

#include <dune/pdelab/test/gridexamples.hh>

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
    //return true;
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
    //return true; // FIXME
    if (entity.partitionType() == Dune::PartitionType::GhostEntity) // NOTE: Replaced by EntitySet excluder
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

template <typename Vector>
class EntitySetExcluder {
public:

  void setPartUnity(std::shared_ptr<Vector> partUnity) {
    partUnity_ = partUnity;
  }

  template<typename Entity>
  bool includeEntity(const Entity& entity) const {
    if (entity.partitionType() == Dune::PartitionType::GhostEntity)
      return false;
    if (partUnity_ == nullptr)
      return true;

   /* for (std::size_t i = 0; i < cache.size(); i++)
    {
      if (partUnity_[cache.containerIndex(i)[0]] > 0.0 &&
        partUnity_[cache.containerIndex(i)[0]] < 1.0)
        return true;
    }*/
    return true;
  }
private:
  //typedef LocalFunctionSpace<GFS, TrialSpaceTag> LFSU;
  std::shared_ptr<Vector> partUnity_ = nullptr;
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

  typedef typename GM::LeafGridView GV_;
  auto gv_ = grid->leafGridView();

  // Transfer partitioning from ParMETIS to our grid
  std::vector<unsigned> part(Dune::ParMetisGridPartitioner<GV_>::partition(gv_, helper));
  grid->loadBalance(part, 0);


  /*typedef Dune::YaspGrid<dim> GM;


  Dune::FieldVector<double,dim> L(1.0);
  std::array<int,dim> N;
  N.fill(cells);
  std::bitset<dim> B(false);

  auto grid = std::make_shared<GM>(L,N,B,overlap,Dune::MPIHelper::getCollectiveCommunication());
  typedef typename GM::LeafGridView GV_;
  auto gv_ = grid->leafGridView();*/

  std::cout << "grid setup: " << timer_detailed.elapsed() << std::endl; timer_detailed.reset();




  const int components = 1;
  using K = double;
  using Vector = Dune::BlockVector<Dune::FieldVector<K,components>>;
  using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<K,components,components>>;

  //using ES = Dune::PDELab::NonOverlappingEntitySet<GV>;
  using ESExcluder = EntitySetExcluder<Vector>;
  ESExcluder es_excluder;
  using GV = Dune::PDELab::OverlapEntitySet<GV_,Dune::Partitions::All, ESExcluder>;
  GV gv(gv_, es_excluder);
  //typedef GV_ GV;
  //GV_& gv = gv_;

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
  //if (gfs.gridView().comm().rank() == 3) verb=3;
  //verb = 3;

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

  //typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,NumberType,NumberType,NumberType,CC,CC> GO;
  //auto go = GO(gfs,cc,gfs,cc,lop,MBE(nonzeros));
  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,NumberType,NumberType,NumberType,CC,CC,GhostExcluder> GO;
  GhostExcluder ghost_excluder;
  auto go = GO(gfs,cc,gfs,cc,lop,MBE(nonzeros),ghost_excluder);


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
  int avg = nonzeros;// (dim==1) ? 3 : ((dim==2) ? 7 : 14);
  //int algebraic_overlap = 10;
  int algebraic_overlap = 1;

  if (verb > 2) {
    Dune::printmatrix(std::cout, native(A), "A", "");
  }

  Dune::NonoverlappingOverlapAdapter<GV, Vector, Matrix> adapter(gv, native(A), avg, algebraic_overlap);
  std::cout << "NonoverlappingOverlapAdapter: " << timer_detailed.elapsed() << std::endl; timer_detailed.reset();

  std::shared_ptr<Matrix> A_extended = adapter.extendMatrix(native(A));
  std::cout << "extend: " << timer_detailed.elapsed() << std::endl; timer_detailed.reset();
  if (verb > 2) {
    Dune::printmatrix(std::cout, native(A), "A", "");
    Dune::printmatrix(std::cout, *A_extended, "A_extended", "");
  }
  std::shared_ptr<Vector> part_unity = Dune::makePartitionOfUnity(adapter, *A_extended);
  std::cout << "part_unity: " << timer_detailed.elapsed() << std::endl; timer_detailed.reset();

  {
    V out(gfs, 0.0);
    Vector basis_restricted(adapter.getExtendedSize());
    adapter.restrictVector(*part_unity, basis_restricted);
    for (int i = 0; i < native(out).N(); i++)
      native(out)[i] = basis_restricted[i];
    Dune::VTKWriter<GV_> vtkwriter(gv_);//gfs.gridView());
    typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
    DGF xdgf(gfs,out);
    typedef Dune::PDELab::VTKGridFunctionAdapter<DGF> ADAPT;
    auto adapt = std::make_shared<ADAPT>(xdgf,"PoU");
    vtkwriter.addVertexData(adapt);
    vtkwriter.write("part_unity");
  }


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



  int nev = 2;
  //auto subdomainbasis = std::make_shared<Dune::PDELab::SubdomainBasis<Vector>>(*Dune::makePartitionOfUnity(adapter_A));


  M newmat(go_overlap);
  auto extended_matrices = adapter.lambdaMultiExtendMatrix(native(A_ovlp), native(A), [&](int i){

    std::shared_ptr<Vector> neighbor_part_unity = remotePartUnities.getVectorForRank(i); //.neighbor_basis[i];

    adapter.restrictVector(*neighbor_part_unity, *part_unity_restricted);

    InteriorExcluder<Vector> excluder(*part_unity_restricted);
    auto go_overlap = GO_OVLP(gfs,cc,gfs,cc,lop,MBE(nonzeros),excluder);

    // TODO: es_excluder.setPartUnity(part_unity_restricted);
    es_excluder.setPartUnity(part_unity_restricted);

    for (auto rIt=native(newmat).begin(); rIt!=native(newmat).end(); ++rIt)
      for (auto cIt=rIt->begin(); cIt!=rIt->end(); ++cIt) {
        *cIt = 0.0;
      }

    go_overlap.jacobian(x,newmat);

    return stackobject_to_shared_ptr(native(newmat));
  });

  std::shared_ptr<Matrix> A_ovlp_extended = extended_matrices.first;
  A_extended = extended_matrices.second;

  //std::shared_ptr<Matrix> A_ovlp_extended = adapter.multiExtendMatrix(native(A_ovlp), A_ovlp_for_neighbors);
  //A_extended = adapter.multiExtendMatrix(native(A), A_ovlp_for_neighbors); // NOTE: Need to redo extension... chicken and egg with partition of unity
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

  auto subdomainbasis = std::make_shared<Dune::PDELab::NewGenEOBasis<GV, Matrix, Vector>>(adapter, *A_extended, *A_ovlp_extended, *part_unity, -1.0, nev);
  //auto subdomainbasis = std::make_shared<Dune::PDELab::SubdomainBasis<Vector>>(*part_unity);
  std::cout << "eigenproblems: " << timer_detailed.elapsed() << std::endl; timer_detailed.reset();


  {
    V out(gfs, 0.0);
    Vector basis_restricted(adapter.getExtendedSize());
    adapter.restrictVector(*part_unity, basis_restricted);
    for (int i = 0; i < native(out).N(); i++)
      native(out)[i] = basis_restricted[i];
    Dune::VTKWriter<GV_> vtkwriter(gv_);
    typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
    DGF xdgf(gfs,out);
    typedef Dune::PDELab::VTKGridFunctionAdapter<DGF> ADAPT;
    auto adapt = std::make_shared<ADAPT>(xdgf,"PoU");
    vtkwriter.addVertexData(adapt);
    vtkwriter.write("newpart_unity_with_dir");
  }

  if(verb > 2) {
    Dune::printvector(std::cout, *part_unity, "part_unity with Dirichlet", "");
    Dune::printmatrix(std::cout, native(A), "A", "");
    Dune::printmatrix(std::cout, *A_extended, "A_extended", "");
    Dune::printmatrix(std::cout, native(A_ovlp), "A_ovlp", "");
    Dune::printmatrix(std::cout, *A_ovlp_extended, "A_ovlp_extended", "");
  }

  auto coarse_space = std::make_shared<Dune::PDELab::NewSubdomainProjectedCoarseSpace<GV, Matrix, Vector>>(adapter, gv, *A_extended, subdomainbasis, verbose);

  if(verb > 2) {
    Dune::printmatrix(std::cout, *coarse_space->get_coarse_system(), "coarse_system", "");
  }

  for (int i = 0; i < subdomainbasis->basis_size(); i++) {
    if (verb > 2)
      Dune::printvector(std::cout, *subdomainbasis->get_basis_vector(i), "basis vector " + std::to_string(i), "", 1, 10, 17);

    V out(gfs, 0.0);
    Vector basis_restricted(adapter.getExtendedSize());
    adapter.restrictVector(*subdomainbasis->get_basis_vector(i), basis_restricted);
    for (int i = 0; i < native(out).N(); i++)
      native(out)[i] = basis_restricted[i];
    Dune::VTKWriter<GV_> vtkwriter(gv_);
    typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
    DGF xdgf(gfs,out);
    typedef Dune::PDELab::VTKGridFunctionAdapter<DGF> ADAPT;
    auto adapt = std::make_shared<ADAPT>(xdgf,"basis vector");
    vtkwriter.addVertexData(adapt);
    vtkwriter.write("newtestgeneo_basis_" + std::to_string(i));
  }


  // Apply Dirichlet conditions on processor boundaries, needed for Schwarz method
  for (auto rIt=A_extended->begin(); rIt!=A_extended->end(); ++rIt)
    for (auto cIt=rIt->begin(); cIt!=rIt->end(); ++cIt)
    {
      if ((*part_unity)[rIt.index()] == .0) {
        *cIt = rIt.index() == cIt.index() ? 1.0 : 0.0;
      }
    }

  Dune::PDELab::ISTL::NewTwoLevelOverlappingAdditiveSchwarz<GV, Matrix, Vector> prec(adapter, *A_extended, coarse_space, true, verb);


  NonoverlappingOperator<GV, Matrix,Vector> linearOperator(gv,native(A));
  NonoverlappingScalarProduct<GV,Vector> scalarproduct(gv,native(x));
  Dune::CGSolver<Vector> solver(linearOperator,scalarproduct,prec,1e-6,500,verbose);

  //Vector b_cpy(d);
  Vector v(native(x)); // TODO: init via size
  // FIXME: Dirichlet unterstützen
  Dune::InverseOperatorResult stat;
  solver.apply(native(v),native(d),stat);
  native(x) -= v;

  //native(x) -= x_cpy;

  // now solve defect equation A*v = d using a CG solver with our shiny preconditioner
  /*V v(gfs,0.0);
  auto solver_ref = std::make_shared<Dune::CGSolver<V> >(*popf,ospf,*prec,1E-6,1000,verb,true);
  Dune::InverseOperatorResult result;
  solver_ref->apply(v,d,result);
  x -= v;*/


  // Write solution to VTK
  Dune::VTKWriter<GV_> vtkwriter(gv_);//gfs.gridView());
  typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
  DGF xdgf(gfs,x);
  typedef Dune::PDELab::VTKGridFunctionAdapter<DGF> ADAPT;
  auto adapt = std::make_shared<ADAPT>(xdgf,"solution");
  vtkwriter.addVertexData(adapt);
  vtkwriter.write("newtestgeneo_basis_" + basis_type + "_part_unity_" + part_unity_type);
}


int main(int argc, char **argv)
{
  using Dune::PDELab::Backend::native;

  try{
    // initialize MPI, finalize is done automatically on exit
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc,argv);

    driver("geneo", "standard", helper);
    //driver("geneo", "sarkis");
    //driver("lipton_babuska", "standard");

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
