#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/pdelab/boilerplate/pdelab.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionadapter.hh>
#include <dune/pdelab/localoperator/convectiondiffusionfem.hh>
#include <dune/pdelab/gridfunctionspace/vtk.hh>

#include <dune/pdelab/preconditioner/geneo/geneo.hh>

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
    RF perm2 = 1e6;
    RF layer_thickness = 1.0 / 20.0;

    RF coeff = (int)std::floor(xglobal[1] / layer_thickness) % 2 == 0 ? perm1 : perm2;

    typename Traits::PermTensorType I;
    I[0][0] = coeff;
    I[0][1] = 0;
    I[1][0] = 0;
    I[1][1] = coeff;
    return I;
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
    if (!(xglobal[1]<1E-6 || xglobal[1]>1.0-1E-6))
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
    else
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
  }

  //! Dirichlet boundary condition value
  typename Traits::RangeFieldType
  g (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::DomainType xglobal = e.geometry().global(x);
    if (xglobal[1] > 1.0-1E-6)
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
};

void driver(std::string basis_type, std::string part_unity_type) {

  int cells = 100;
  int overlap = 1;

  // define parameters
  const unsigned int dim = 2;
  const unsigned int degree = 1;
  const std::size_t nonzeros = std::pow(2*degree+1,dim);
  const Dune::GeometryType::BasicType elemtype = Dune::GeometryType::cube;
  const Dune::PDELab::MeshType meshtype = Dune::PDELab::MeshType::conforming;
  const Dune::SolverCategory::Category solvertype = Dune::SolverCategory::overlapping;
  typedef double NumberType;

  // build a grid
  typedef Dune::YaspGrid<dim> GM;

  Dune::FieldVector<double,dim> L(1.0);
  std::array<int,dim> N = {cells, cells};
  std::bitset<dim> B(false);

  typedef Dune::YaspFixedSizePartitioner<dim> YP;
  std::array<int,dim> yasppartitions = {2, 1};
  auto yp = new YP(yasppartitions);

  auto grid = std::make_shared<GM>(L,N,B,overlap,Dune::MPIHelper::getCollectiveCommunication(),yp);


  // make problem parameters
  typedef GenericEllipticProblem<GM::LevelGridView,NumberType> Problem;
  Problem problem;
  typedef Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<Problem> BCType;
  BCType bctype(grid->levelGridView(grid->maxLevel()),problem);


  // make a finite element space
  typedef typename GM::LevelGridView GV;
  typedef typename GM::ctype ctype;
  static const int dimworld = GM::dimensionworld;


  typedef Dune::PDELab::CGFEMBase<GV,ctype,NumberType,degree,dim,elemtype> FEMB;
  typedef Dune::PDELab::CGCONBase<GM,degree,elemtype,meshtype,solvertype,BCType> CONB;

  typedef typename FEMB::FEM FEM;
  typedef typename CONB::CON CON;
  typedef Dune::PDELab::ISTL::VectorBackend<> VBE;

  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
  auto gv = grid->levelGridView(grid->maxLevel());
  auto fem = FEMB(gv);
  auto con = CONB(*grid, bctype);
  auto gfs = GFS(gv, fem.getFEM(), con.getCON()); gfs.name("Solution");
  con.postGFSHook(gfs);

  // make a degree of freedom vector on fine grid and initialize it with interpolation of Dirichlet condition
  typedef Dune::PDELab::Backend::Vector<GFS,NumberType> V;
  V x(gfs,0.0);
  typedef Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<Problem> G;
  G g(grid->levelGridView(grid->maxLevel()),problem);
  Dune::PDELab::interpolate(g,gfs,x);

  // Set up constraints containers
  typedef typename GFS::template ConstraintsContainer<NumberType>::Type CC;
  auto cc = CC();
  auto cc_exterior = CC();
  auto cc_bnd_neu_int_dir = CC();

  // assemble constraints
  Dune::PDELab::constraints(bctype,gfs,cc);
  Dune::PDELab::constraints_exterior(bctype,gfs,cc_exterior);

  Dune::PDELab::NoDirichletConstraintsParameters pnbc;
  Dune::PDELab::constraints(pnbc,gfs,cc_bnd_neu_int_dir);

  // set initial guess
  V x0(gfs,0.0);
  {
    Dune::PDELab::copy_nonconstrained_dofs(cc,x0,x);
    con.make_consistent(gfs,x);
  }

  // assembler for finite element problem
  typedef Dune::PDELab::ConvectionDiffusionFEM<Problem,FEM> LOP;
  LOP lop(problem);

  typedef Dune::PDELab::ISTL::BCRSMatrixBackend<> MBE;
  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,NumberType,NumberType,NumberType,CC,CC> GO;
  auto go = GO(gfs,cc,gfs,cc,lop,MBE(nonzeros));

  // set up and assemble right hand side w.r.t. l(v)-a(u_g,v)
  V d(gfs,0.0);
  go.residual(x,d);


  // types
  typedef GO::Jacobian M;

  // fine grid objects
  M AF(go);
  go.jacobian(x,AF); // assemble fine grid matrix
  typedef Dune::PDELab::OverlappingOperator<CC,M,V,V> POP;
  auto popf = std::make_shared<POP>(cc,AF);
  typedef Dune::PDELab::ISTL::ParallelHelper<GFS> PIH;
  PIH pihf(gfs);
  typedef Dune::PDELab::OverlappingScalarProduct<GFS,V> OSP;
  OSP ospf(gfs,pihf);

  // Compute system
  auto go_exterior = GO(gfs,cc_exterior,gfs,cc_exterior,lop,MBE(nonzeros));
  M AF_exterior(go_exterior);
  go_exterior.jacobian(x,AF_exterior); // assemble fine grid matrix

  // Compute overlap system using a wrapper for the LocalOperator
  typedef Dune::PDELab::LocalOperatorOvlpRegion<LOP, GFS> LOP_OVLP;
  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP_OVLP,MBE,NumberType,NumberType,NumberType,CC,CC> GO_OVLP;
  LOP_OVLP lop_ovlp(lop, gfs);
  auto go_overlap = GO_OVLP(gfs,cc_exterior,gfs,cc_exterior,lop_ovlp,MBE(nonzeros));
  M AF_ovlp(go_overlap);
  go_overlap.jacobian(x,AF_ovlp); // assemble fine grid matrix

  // Choose an eigenvalue threshold according to Spillane et al., 2014
  double eigenvalue_threshold = (double)overlap / (cells + overlap);

  int verb=0;
  if (gfs.gridView().comm().rank()==0) verb=2;
  typedef Dune::PDELab::LocalFunctionSpace<GFS, Dune::PDELab::AnySpaceTag> LFS;
  LFS lfs(gfs);

  std::shared_ptr<V> part_unity;
  if (part_unity_type == "standard")
    part_unity = standardPartitionOfUnity<V>(gfs, lfs, cc_bnd_neu_int_dir);
  else if (part_unity_type == "sarkis")
    part_unity = sarkisPartitionOfUnity<V>(gfs, lfs, cc_bnd_neu_int_dir, cells, cells, overlap, yasppartitions[0], yasppartitions[1]);
  else
    DUNE_THROW(Dune::Exception, "Unkown selection in test driver!");

  int nev = 10;
  int nev_arpack = 10;

  std::shared_ptr<Dune::PDELab::SubdomainBasis<V> > subdomain_basis;
  if (basis_type == "geneo")
    subdomain_basis = std::make_shared<Dune::PDELab::GenEOBasis<GFS,M,V,1> >(gfs, AF_exterior, AF_ovlp, eigenvalue_threshold, *part_unity, nev, nev_arpack);
  else if (basis_type == "lipton_babuska")
    subdomain_basis = std::make_shared<Dune::PDELab::LiptonBabuskaBasis<GFS,M,V,V,1> >(gfs, AF_exterior, AF_ovlp, -1, *part_unity, nev, nev_arpack);
  else if (basis_type == "part_unity") // We can't test this one, it does not lead to sufficient error reduction. Let's instantiate it anyway for test's sake.
    subdomain_basis = std::make_shared<Dune::PDELab::OneDimensionalSubdomainBasis<V> >(*part_unity);
  else
    DUNE_THROW(Dune::Exception, "Unkown selection in test driver!");


  auto partunityspace = std::make_shared<Dune::PDELab::SubdomainProjectedCoarseSpace<GFS,M,V,1> >(gfs, AF_exterior, subdomain_basis, verb);
  auto prec = std::make_shared<Dune::PDELab::ISTL::TwoLevelOverlappingAdditiveSchwarz<GFS,M,V,V>>(gfs, AF, partunityspace);


  // now solve defect equation A*v = d
  V v(gfs,0.0);
  auto solver_ref = std::make_shared<Dune::CGSolver<V> >(*popf,ospf,*prec,1E-6,1000,verb,true);
  Dune::InverseOperatorResult result;
  solver_ref->apply(v,d,result);
  x -= v;


  // Write solution to VTK
  Dune::VTKWriter<GV> vtkwriter(gfs.gridView());
  typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
  DGF xdgf(gfs,x);
  typedef Dune::PDELab::VTKGridFunctionAdapter<DGF> ADAPT;
  auto adapt = std::make_shared<ADAPT>(xdgf,"solution");
  vtkwriter.addVertexData(adapt);
  vtkwriter.write("basis_" + basis_type + "_part_unity_" + part_unity_type);
}


int main(int argc, char **argv)
{
  using Dune::PDELab::Backend::native;

  try{
    // initialize MPI, finalize is done automatically on exit
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc,argv);

    driver("geneo", "standard");
    driver("geneo", "sarkis");
    driver("lipton_babuska", "standard");

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
