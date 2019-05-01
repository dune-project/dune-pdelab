#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/common/parametertree.hh>
Dune::ParameterTree configuration;

#include <dune/pdelab/boilerplate/pdelab.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionadapter.hh>
#include <dune/pdelab/localoperator/convectiondiffusionfem.hh>
#include <dune/pdelab/gridfunctionspace/vtk.hh>

#include <dune/pdelab/backend/istl/geneo/geneo.hh>

/*
 * Defining a Darcy problem with alternating layers of permeability and a high contrast
 */
template<typename GV, typename RF>
class GenericEllipticProblem
{
  typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;

public:
  typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

  GenericEllipticProblem() {
    perm1 = 1e0;
    perm2 = configuration.get<double>("contrast");
    layer_thickness = 1.0 / (double)configuration.get<int>("layers");
    layer_model = configuration.get<bool>("layer_model");
  }

  //! tensor diffusion constant per cell? return false if you want more than one evaluation of A per cell.
  static constexpr bool permeabilityIsConstantPerCell()
  {
    return true;
  }

  //! tensor diffusion coefficient
  typename Traits::PermTensorType
  A (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    if (layer_model) {
      typename Traits::DomainType xglobal = e.geometry().global(x);

      RF coeff = (int)std::floor(xglobal[1] / layer_thickness) % 2 == 0 ? perm1 : perm2;

      typename Traits::PermTensorType I;
      I[0][0] = coeff;
      I[0][1] = 0;
      I[1][0] = 0;
      I[1][1] = coeff;
      return I;
    } else {

      typename Traits::DomainType xglobal = e.geometry().global(x);
      int num1 = floor(8*xglobal[0]);
      int num2 = floor(8*xglobal[1]);
      RF coeff = 0.0;
      if ( (num1 % 2 == 0) && (num2 % 2 == 0) )
        coeff = 0.1*perm2*(num2+1.0);
      else
        coeff = 1.0;

      RF duct1_b = 0.9*xglobal[0];
      RF duct1_t = duct1_b + 0.1;
      if ( (xglobal[1]>duct1_b) && (xglobal[1]<duct1_t))
        coeff = coeff + 0.2 * perm2;

      RF duct2_b = -1.0*xglobal[0] + 0.5;
      RF duct2_t = duct2_b + 0.1;
      if ( (xglobal[1]>duct2_b) && (xglobal[1]<duct2_t))
        coeff = coeff + 0.5 * perm2;

      RF duct3_b = 4.0*xglobal[0] - 2.0 ;
      RF duct3_t = duct3_b + 0.1;
      if ( (xglobal[1]>duct3_b) && (xglobal[1]<duct3_t))
        coeff = coeff + 0.2 * perm2;

      RF eps = 1.0;
      RF th = 0;
      typename Traits::PermTensorType I;
      th=th*M_PI/180.0;
      I[0][0]=coeff*(pow(cos(th),2.0) + pow(sin(th),2.0)*eps);
      I[0][1]=coeff*sin(2.0*th)*(eps-1.0)/2.0;
      I[1][0]=coeff*sin(2.0*th)*(eps-1.0)/2.0;
      I[1][1]=coeff*(pow(sin(th),2.0) + pow(cos(th),2.0)*eps);
      return I;
    }
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
    return 0.0;
  }

  BCType
  bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    typename Traits::DomainType xglobal = is.geometry().global(x);
    if (!(xglobal[1]<1E-5 || xglobal[1]>1.0-1E-5))
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
    else
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
  }

  //! Dirichlet boundary condition value
  typename Traits::RangeFieldType
  g (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::DomainType xglobal = e.geometry().global(x);
    if (xglobal[1] > 1.0-1E-5)
      return 1.0;
    else
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
private:
  RF perm1, perm2, layer_thickness;
  bool layer_model;
};

void driver(std::string basis_type, std::string part_unity_type) {

  Dune::Timer timer_full;

  int cells = configuration.get<int>("cells");
  int overlap = configuration.get<int>("overlap");

  // define parameters
  const unsigned int dim = 2;
  const unsigned int degree = 1;
  const std::size_t nonzeros = std::pow(2*degree+1,dim);
  typedef double NumberType;

  // build a grid
  typedef Dune::YaspGrid<dim> GM;

  Dune::FieldVector<double,dim> L(1.0);
  std::array<int,dim> N = {cells, cells};
  std::bitset<dim> B(false);

  typedef Dune::YaspFixedSizePartitioner<dim> YP;
  std::array<int,dim> yasppartitions;

  std::shared_ptr<GM> grid = nullptr;
  if (configuration.get<bool>("extend_domain_with_procs")) {
    yasppartitions = {Dune::MPIHelper::getCollectiveCommunication().size() / 2, 2};
    L[0] *= Dune::MPIHelper::getCollectiveCommunication().size();
    auto yp = new YP(yasppartitions);
    grid = std::make_shared<GM>(L,N,B,overlap,Dune::MPIHelper::getCollectiveCommunication(),yp);
  } else {
    grid = std::make_shared<GM>(L,N,B,overlap,Dune::MPIHelper::getCollectiveCommunication());
  }

  // make problem parameters
  typedef GenericEllipticProblem<GM::LevelGridView,NumberType> Problem;
  Problem problem;
  typedef Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<Problem> BCType;
  BCType bctype(grid->levelGridView(grid->maxLevel()),problem);


  // make a finite element space
  typedef typename GM::LevelGridView GV;

  auto gv = grid->levelGridView(grid->maxLevel());


  typedef typename GV::Grid::ctype DF;
  // instantiate finite element maps
  typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,double,1> FEM;
  FEM fem(gv);

  // make function space, with overlapping constraints as usual for overlapping Schwarz methods
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,
                                          Dune::PDELab::OverlappingConformingDirichletConstraints,
                                          Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed,1>
                                          > GFS;
  GFS gfs(gv,fem);


  // and a second function space with no constraints on processor boundaries, needed for the GenEO eigenproblem
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,
                                          Dune::PDELab::ConformingDirichletConstraints,
                                          Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed,1>
                                          > GFS_EXTERIOR;
  GFS_EXTERIOR gfs_exterior(gv,fem);


  // make a degree of freedom vector on fine grid and initialize it with interpolation of Dirichlet condition
  typedef Dune::PDELab::Backend::Vector<GFS,NumberType> V;
  V x(gfs,0.0);

  // also create a wrapper using the same data based on the GFS without processor boundary constraints
  typedef Dune::PDELab::Backend::Vector<GFS_EXTERIOR,NumberType> V_EXTERIOR;
  V_EXTERIOR x_exterior(gfs_exterior, Dune::PDELab::Backend::unattached_container());
  x_exterior.attach(x.storage());

  // Extract domain boundary constraints from problem definition, apply trace to solution vector
  typedef Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<Problem> G;
  G g(grid->levelGridView(grid->maxLevel()),problem);
  Dune::PDELab::interpolate(g,gfs,x);


  // Set up constraints containers
  //  - with boundary constraints and processor constraints as usual
  //  - with boundary constraints, but without processor constraints for the eigenprobleme
  //  - with only Neumann boundary constraints, but with processor constraints as needed for the partition of unity
  typedef typename GFS::template ConstraintsContainer<NumberType>::Type CC;
  auto cc = CC();
  typedef typename GFS_EXTERIOR::template ConstraintsContainer<NumberType>::Type CC_EXTERIOR;
  auto cc_exterior = CC_EXTERIOR();
  auto cc_bnd_neu_int_dir = CC();

  // assemble constraints
  Dune::PDELab::constraints(bctype,gfs,cc);
  Dune::PDELab::constraints(bctype,gfs_exterior,cc_exterior);

  Dune::PDELab::NoDirichletConstraintsParameters pnbc;
  Dune::PDELab::constraints(pnbc,gfs,cc_bnd_neu_int_dir);

  // set initial guess
  V x0(gfs,0.0);
  Dune::PDELab::copy_nonconstrained_dofs(cc,x0,x);


  // LocalOperator for given problem
  typedef Dune::PDELab::ConvectionDiffusionFEM<Problem,FEM> LOP;
  LOP lop(problem);

  // LocalOperator wrapper zeroing out subdomains' interiors in order to set up overlap matrix
  typedef Dune::PDELab::LocalOperatorOvlpRegion<LOP, GFS> LOP_OVLP;
  LOP_OVLP lop_ovlp(lop, gfs);

  typedef Dune::PDELab::ISTL::BCRSMatrixBackend<> MBE;


  // Construct GridOperators from LocalOperators
  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,NumberType,NumberType,NumberType,CC,CC> GO;
  auto go = GO(gfs,cc,gfs,cc,lop,MBE(nonzeros));

  typedef Dune::PDELab::GridOperator<GFS_EXTERIOR,GFS_EXTERIOR,LOP,MBE,NumberType,NumberType,NumberType,CC,CC> GO_EXTERIOR;
  auto go_exterior = GO_EXTERIOR(gfs_exterior,cc_exterior,gfs_exterior,cc_exterior,lop,MBE(nonzeros));

  typedef Dune::PDELab::GridOperator<GFS_EXTERIOR,GFS_EXTERIOR,LOP_OVLP,MBE,NumberType,NumberType,NumberType,CC,CC> GO_OVLP;
  auto go_overlap = GO_OVLP(gfs_exterior,cc_exterior,gfs_exterior,cc_exterior,lop_ovlp,MBE(nonzeros));

  // set up and assemble right hand side w.r.t. l(v)-a(u_g,v)
  V d(gfs,0.0);
  go.residual(x,d);


  // types
  typedef GO::Jacobian M;
  typedef GO_EXTERIOR::Jacobian M_EXTERIOR;

  // fine grid objects
  M AF(go);
  go.jacobian(x,AF);
  typedef Dune::PDELab::OverlappingOperator<CC,M,V,V> POP;
  auto popf = std::make_shared<POP>(cc,AF);
  typedef Dune::PDELab::ISTL::ParallelHelper<GFS> PIH;
  PIH pihf(gfs);
  typedef Dune::PDELab::OverlappingScalarProduct<GFS,V> OSP;
  OSP ospf(gfs,pihf);

  // Assemble fine grid matrix defined without processor constraints
  M_EXTERIOR AF_exterior(go_exterior);
  go_exterior.jacobian(x_exterior,AF_exterior);

  // Assemble fine grid matrix defined only on overlap region
  M_EXTERIOR AF_ovlp(go_overlap);
  go_overlap.jacobian(x_exterior,AF_ovlp);


  // Choose an eigenvalue threshold according to Spillane et al., 2014.
  // This particular choice is a heuristic working very well for Darcy problems.
  // Theoretically, any value delivers robustness; practically, you may need to
  // choose another value to achieve a good balance between condition bound and
  // global basis size
  double eigenvalue_threshold = (double)overlap / (cells + overlap);
  if (!configuration.get<bool>("use_threshold"))
    eigenvalue_threshold = -1;

  int verb=0;
  if (gfs.gridView().comm().rank()==0) verb=2;

  // Create a local function space needed for the Sarkis partition of unity only
  typedef Dune::PDELab::LocalFunctionSpace<GFS, Dune::PDELab::AnySpaceTag> LFS;
  LFS lfs(gfs);


  // Generate a partition of unity
  std::shared_ptr<V> part_unity;
  if (part_unity_type == "standard")
    part_unity = std::make_shared<V>(standardPartitionOfUnity<V>(gfs, cc_bnd_neu_int_dir));
  else if (part_unity_type == "sarkis")
    part_unity = std::make_shared<V>(sarkisPartitionOfUnity<V>(gfs, lfs, cc_bnd_neu_int_dir, cells, cells, overlap, yasppartitions[0], yasppartitions[1]));
  else
    DUNE_THROW(Dune::Exception, "Unkown selection in test driver!");


  // Choose how many eigenvalues to compute
  int nev = configuration.get<int>("nev");
  int nev_arpack = configuration.get<int>("nev_arpack");

  if (verb > 0) std::cout << "Basis setup" << std::endl;
  Dune::Timer timer_setup;

  // Construct per-subdomain basis functions
  std::shared_ptr<Dune::PDELab::SubdomainBasis<V> > subdomain_basis;
  if (basis_type == "geneo")
    subdomain_basis = std::make_shared<Dune::PDELab::GenEOBasis<GFS,M_EXTERIOR,V,1> >(gfs, AF_exterior, AF_ovlp, eigenvalue_threshold, *part_unity, nev, nev_arpack, 0.001, false, verb);
  //else if (basis_type == "lipton_babuska")
  //  subdomain_basis = std::make_shared<Dune::PDELab::LiptonBabuskaBasis<GFS,M_EXTERIOR,V,V,1> >(gfs, AF_exterior, AF_ovlp, -1, *part_unity, nev, nev_arpack);
  else if (basis_type == "part_unity") // We can't test this one, it does not lead to sufficient error reduction. Let's instantiate it anyway for test's sake.
    subdomain_basis = std::make_shared<Dune::PDELab::SubdomainBasis<V> >(*part_unity);
  else
    DUNE_THROW(Dune::Exception, "Unkown selection in test driver!");

  if (verb > 0) std::cout << "Basis setup finished: G=" << timer_setup.elapsed() << std::endl;

  // Fuse per-subdomain basis functions to a global coarse space
  auto coarse_space = std::make_shared<Dune::PDELab::SubdomainProjectedCoarseSpace<GFS,M_EXTERIOR,V,PIH> >(gfs, AF_exterior, subdomain_basis, pihf);


  // Plug coarse basis into actual preconditioner
  auto prec = std::make_shared<Dune::PDELab::ISTL::TwoLevelOverlappingAdditiveSchwarz<GFS,M,V,V>>(gfs, AF, coarse_space, true, verb);


  // now solve defect equation A*v = d using a CG solver with our shiny preconditioner
  V v(gfs,0.0);
  auto solver_ref = std::make_shared<Dune::CGSolver<V> >(*popf,ospf,*prec,1E-6,1000,verb,configuration.get<bool>("condition_estimate"));
  Dune::InverseOperatorResult result;
  solver_ref->apply(v,d,result);
  x -= v;

  if (verb > 0) std::cout << "CG solver finished: S=" << result.elapsed << std::endl;

  if (verb > 0) std::cout << "Solver finished: F=" << timer_full.elapsed() << std::endl;

  // Write solution to VTK
  Dune::VTKWriter<GV> vtkwriter(gfs.gridView());
  typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
  DGF xdgf(gfs,x);
  typedef Dune::PDELab::VTKGridFunctionAdapter<DGF> ADAPT;
  auto adapt = std::make_shared<ADAPT>(xdgf,"solution");
  vtkwriter.addVertexData(adapt);
  vtkwriter.write("testgeneo_basis_" + basis_type + "_part_unity_" + part_unity_type);
}


int main(int argc, char **argv)
{
  Dune::ParameterTreeParser parser;
  parser.readINITree("config.ini",configuration);
  parser.readOptions(argc, argv, configuration);

  using Dune::PDELab::Backend::native;

  try{
    // initialize MPI, finalize is done automatically on exit
    Dune::MPIHelper::instance(argc,argv);

    driver("geneo", "standard");
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
