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




template <typename Vector, typename GV>
class EntitySetExcluder {
public:

  typedef typename GV::template Codim<0>::Entity Entity;
  typedef typename GV::template Codim<1>::Entity Entity1;
  typedef typename GV::template Codim<2>::Entity Entity2;

  bool includeEntity(const Entity1& entity) const {
    return true;
  }
  bool includeEntity(const Entity2& entity) const {
    return true;
  }

  virtual bool includeEntity(const Entity& entity) const {
    if (entity.partitionType() == Dune::PartitionType::GhostEntity)
      return false;

    return true;
  }

};

template <typename Vector, typename GV, typename GFS>
class EntitySetPartUnityExcluder : public EntitySetExcluder<Vector, GV> {
public:

  typedef typename GV::template Codim<0>::Entity Entity;

  EntitySetPartUnityExcluder(const GFS& gfs, std::shared_ptr<Vector> partUnity) : lfs(gfs), partUnity_(partUnity) {}

  bool includeEntity(const Entity& entity) const override {

    if (entity.partitionType() == Dune::PartitionType::GhostEntity)
      return false;

    typedef Dune::PDELab::LFSIndexCache<LFS,Dune::PDELab::EmptyTransformation> LFSCache;
    LFSCache lfs_cache(lfs);
    lfs.bind( entity );
    lfs_cache.update();

    for (std::size_t i = 0; i < lfs_cache.size(); i++)
    {
      if ((*partUnity_)[lfs_cache.containerIndex(i).back()] > 0.0 &&
        (*partUnity_)[lfs_cache.containerIndex(i).back()] < 1.0) {
        return true;
        }
    }
    return false;
  }

private:
  typedef Dune::PDELab::LocalFunctionSpace<GFS, Dune::PDELab::TrialSpaceTag> LFS;
  mutable LFS lfs;
  std::shared_ptr<Vector> partUnity_ = nullptr;
};

template <typename Vector, typename GV>
class EntitySetGhostExcluder : public EntitySetExcluder<Vector, GV> {
public:

  typedef typename GV::template Codim<0>::Entity Entity;

  bool includeEntity(const Entity& entity) const override {
    return entity.partitionType() != Dune::PartitionType::GhostEntity;
  }
};


template <typename GO, typename ScalarMatrix, typename Matrix, typename ScalarVector, typename Vector>
class GenEOPreconditioner : public Dune::Preconditioner<Vector,Vector> {

  using V = typename GO::Domain;
  using M = typename GO::Jacobian;
  using GFS = typename GO::Traits::TrialGridFunctionSpace;
  using GV = typename GFS::Traits::GridView;


public:

  // define the category
  virtual Dune::SolverCategory::Category category() const override
  {
    return prec->category();
  }


  /*!
   *          \brief Prepare the preconditioner.
   *
   *          \copydoc Preconditioner::pre(X&,Y&)
   */
  virtual void pre (Vector& x, Vector& b) override
  {
    prec->pre(x, b);
  }

  /*!
   *          \brief Apply the precondioner.
   *
   *          \copydoc Preconditioner::apply(X&,const Y&)
   */
  virtual void apply (Vector& v, const Vector& d) override
  {
    prec->apply(v, d);
  }
  /*!
   *          \brief Clean up.
   *
   *          \copydoc Preconditioner::post(X&)
   */
  virtual void post (Vector& x) override
  {
    prec->post(x);
  }

  GenEOPreconditioner(const GO& go, const typename GO::Jacobian& A, int algebraic_overlap, int avg_nonzeros, const double eigenvalue_threshold, int& nev,
                      int nev_arpack = -1, const double shift = 0.001, int verbose = 0)
   : A_ovlp(go){

    V x(go.trialGridFunctionSpace(),0.0); // NOTE: We assume linear problems, so simply set x to zero here!

    const GFS& gfs = go.trialGridFunctionSpace();
    const GV& gv = gfs.gridView();


    using Dune::PDELab::Backend::native;


    Dune::NonoverlappingOverlapAdapter<GV, Vector, Matrix> adapter(gv, native(A), avg_nonzeros, algebraic_overlap);

    A_extended = adapter.extendMatrix(native(A));
    part_unity = Dune::makePartitionOfUnity(adapter, *A_extended);


    using Attribute = Dune::EPISAttribute;
    Dune::AllSet<Attribute> allAttribute;
    auto allinterface = std::shared_ptr<Dune::Interface>(new Dune::Interface());
    allinterface->build(*adapter.getRemoteIndices(),allAttribute,allAttribute); // all to all communication
    auto communicator = std::shared_ptr<Dune::BufferedCommunicator>(new Dune::BufferedCommunicator());
    communicator->build<Vector>(*allinterface);

    Dune::PDELab::MultiVectorBundle<GV, Vector, Matrix> remotePartUnities(adapter);
    remotePartUnities.localVector_ = part_unity;
    communicator->forward<Dune::PDELab::MultiGatherScatter<Dune::PDELab::MultiVectorBundle<GV, Vector, Matrix>>>(remotePartUnities,remotePartUnities); // make function known in other subdomains


    // Assemble fine grid matrix defined only on overlap region
    auto part_unity_restricted = std::make_shared<Vector>(native(A).N());
    adapter.restrictVector(*part_unity, *part_unity_restricted);



    auto es_pou_excluder = std::make_shared<EntitySetPartUnityExcluder<Vector, GV, GFS>> (gfs, part_unity_restricted);
    gfs.entitySet().setExcluder(es_pou_excluder);



    go.jacobian(x,A_ovlp);



    M newmat(go);
    // Provide neighbors with matrices assembled exclusively on respective overlap area
    auto extended_matrices = adapter.lambdaMultiExtendMatrix(native(A_ovlp), native(A), [&](int i){
      std::shared_ptr<Vector> neighbor_part_unity = remotePartUnities.getVectorForRank(i);

      adapter.restrictVector(*neighbor_part_unity, *part_unity_restricted);

      std::shared_ptr<EntitySetExcluder<Vector, GV>> es_local_pou_excluder = std::make_shared<EntitySetPartUnityExcluder<Vector, GV, GFS>> (gfs, part_unity_restricted);
      gfs.entitySet().setExcluder(es_local_pou_excluder);

      for (auto rIt=native(newmat).begin(); rIt!=native(newmat).end(); ++rIt)
        for (auto cIt=rIt->begin(); cIt!=rIt->end(); ++cIt) {
          *cIt = 0.0;
        }

      go.jacobian(x,newmat);

      return stackobject_to_shared_ptr(native(newmat));
    });

    A_ovlp_extended = extended_matrices.first;
    A_extended = extended_matrices.second;


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

    auto subdomainbasis = std::make_shared<Dune::PDELab::NewGenEOBasis<GV, Matrix, Vector>>(adapter, *A_extended, *A_ovlp_extended, *part_unity, eigenvalue_threshold, nev, nev_arpack, shift);

    auto coarse_space = std::make_shared<Dune::PDELab::NewSubdomainProjectedCoarseSpace<GV, Matrix, Vector>>(adapter, gv, *A_extended, subdomainbasis, verbose);


    // Apply Dirichlet conditions on processor boundaries, needed for Schwarz method
    for (auto rIt=A_extended->begin(); rIt!=A_extended->end(); ++rIt)
      for (auto cIt=rIt->begin(); cIt!=rIt->end(); ++cIt)
      {
        if ((*part_unity)[rIt.index()] == .0) {
          *cIt = rIt.index() == cIt.index() ? 1.0 : 0.0;
        }
      }

    prec = std::make_shared<Dune::PDELab::ISTL::NewTwoLevelOverlappingAdditiveSchwarz<GV, ScalarMatrix, Matrix, ScalarVector, Vector>>(adapter, *A_extended, coarse_space, true, verbose);

  }

private:

  M A_ovlp;
  std::shared_ptr<Matrix> A_ovlp_extended;
  std::shared_ptr<Matrix> A_extended;
  std::shared_ptr<Vector> part_unity;

  std::shared_ptr<Dune::PDELab::ISTL::NewTwoLevelOverlappingAdditiveSchwarz<typename GO::Traits::TrialGridFunctionSpace::Traits::GridView, ScalarMatrix, Matrix, ScalarVector, Vector>> prec = nullptr;

};


void driver(std::string basis_type, std::string part_unity_type, Dune::MPIHelper& helper) {


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
  typedef typename GM::LeafGridView GV;
  auto gv = grid->leafGridView();*/



  const int components = 1;
  using K = double;
  using Vector = Dune::BlockVector<Dune::FieldVector<K,components>>;
  using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<K,components,components>>;

  using ESExcluder = EntitySetExcluder<Vector, GV>;
  auto ghost_excluder = std::make_shared<EntitySetGhostExcluder<Vector, GV>>();


  using ES = Dune::PDELab::OverlapEntitySet<GV,Dune::Partitions::All, ESExcluder>;
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
  go.residual(x,d); // NOTE: Need operator without Dirichlet on proc boundaries for rhs setup!

  // Choose an eigenvalue threshold according to Spillane et al., 2014.
  // This particular choice is a heuristic working very well for Darcy problems.
  // Theoretically, any value delivers robustness; practically, you may need to
  // choose another value to achieve a good balance between condition bound and
  // global basis size
  //double eigenvalue_threshold = (double)overlap / (cells + overlap);
  double eigenvalue_threshold = -1;


  const int algebraic_overlap = 1;


  int nev = 2;
  auto prec = std::make_shared<GenEOPreconditioner<GO, Matrix, Matrix, Vector, Vector>>(go, A, algebraic_overlap, nonzeros, eigenvalue_threshold, nev, -1, 0.001, verbose);//, eigenvalue_threshold, 2, -1, .001, verbose);



  using Dune::PDELab::Backend::native;

  NonoverlappingOperator<ES, Matrix,Vector> linearOperator(es,native(A));
  NonoverlappingScalarProduct<ES,Vector> scalarproduct(es,native(x));
  Dune::CGSolver<Vector> solver(linearOperator,scalarproduct,*prec,1e-6,500,verbose);

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
  Dune::VTKWriter<GV> vtkwriter(gv);//gfs.gridView());
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
