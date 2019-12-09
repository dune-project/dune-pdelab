#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/common/parametertree.hh>
Dune::ParameterTree configuration;

#include <dune/common/timer.hh>

#include <dune/pdelab/boilerplate/pdelab.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionadapter.hh>
#include <dune/pdelab/localoperator/convectiondiffusionfem.hh>
#include <dune/pdelab/gridfunctionspace/vtk.hh>

#include <dune/pdelab/backend/istl/geneo/geneo.hh>

#include <cholmod.h>

#include <stdio.h>

/*
 * Defining a Darcy problem with alternating layers of permeability and a high contrast
 */
template<typename GV, typename RF>
class GenericEllipticProblem
{
  typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;

public:
  typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

  GenericEllipticProblem()
  : layers(configuration.get<int>("layers")),
    contrast(configuration.get<double>("contrast")) {}

  //! tensor diffusion coefficient
  typename Traits::PermTensorType
  A (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::DomainType xglobal = e.geometry().global(x);

    RF perm1 = 1e0;
    RF perm2 = contrast;
    RF layer_thickness = 1.0 / (double)layers;

    RF coeff = (int)std::floor(xglobal[1]/ layer_thickness) % 2 == 0 ? perm1 : perm2;

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
    typename Traits::DomainType xglobal = e.geometry().global(x);
    return - std::exp(-std::sqrt(std::pow(xglobal[0] - 0.25, 2) + std::pow(xglobal[1] - 0.25, 2)))
           + std::exp(-std::sqrt(std::pow(xglobal[0] - 0.75, 2) + std::pow(xglobal[1] - 0.75, 2)));
    //return 0.0;
  }

  BCType
  bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    typename Traits::DomainType xglobal = is.geometry().global(x);
    /*if (!(xglobal[0]<1E-10 || xglobal[0]>1.0-1E-10))
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
    else
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;*/
    /*if (!((xglobal[1]<0.5 && xglobal[0]>1.0-1E-10) || (xglobal[1]>0.50 && xglobal[0]<1E-10)))
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
    else
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;*/
    //if (//(xglobal[1] < 1E-10 && xglobal[0] > 0.75) || (xglobal[1] > 1.0-1E-10 && xglobal[0] < 0.25) ||
        //xglobal[0] < 1E-10 || xglobal[0] > 1.0-1E-10)
    if (xglobal[1] > 1.0-1E-10 || xglobal[1] < 1E-10)
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
    else
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
  }

  //! Dirichlet boundary condition value
  typename Traits::RangeFieldType
  g (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::DomainType xglobal = e.geometry().global(x);
    /*if (xglobal[1] < 1E-10 && xglobal[0] > 0.75)
      return std::cos(xglobal[0] * 8* 3.141592);
    else if (xglobal[1] > 1.0-1E-10 && xglobal[0] < 0.25)
      return std::sin(xglobal[0] * 8 * 3.141592);
    else if (xglobal[0] < 1E-10)
      return std::sin(xglobal[1] * 2 * 3.141592) + 1.0;
    else if (xglobal[0] > 1.0 - 1E-10)
      return std::cos(xglobal[1] * 2 * 3.141592) - 1.0;
    else
      return 0.0;*/
    if (xglobal[1] > 1.0-1E-10)
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

  int layers;
  double contrast;

};




template<typename T, typename FiniteElementMap>
class DTNRhsLocalOperator :
  public Dune::PDELab::NumericalJacobianApplyVolume<DTNRhsLocalOperator<T,FiniteElementMap> >,
  public Dune::PDELab::NumericalJacobianApplyBoundary<DTNRhsLocalOperator<T,FiniteElementMap> >,
  public Dune::PDELab::NumericalJacobianVolume<DTNRhsLocalOperator<T,FiniteElementMap> >,
  public Dune::PDELab::NumericalJacobianBoundary<DTNRhsLocalOperator<T,FiniteElementMap> >,
  public Dune::PDELab::FullVolumePattern,
  public Dune::PDELab::LocalOperatorDefaultFlags,
  public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<typename T::Traits::RangeFieldType>
{
public:
  // pattern assembly flags
  enum { doPatternVolume = true };
  //enum { doPatternBoundary = false };

  // residual assembly flags
  //enum { doAlphaVolume = false };
  enum { doAlphaBoundary = true };

  DTNRhsLocalOperator (T& param_, int intorderadd_=0)
    : param(param_), intorderadd(intorderadd_)
  {
  }

  // boundary integral
  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_boundary (const IG& ig,
                        const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                        R& r_s) const
  {
    // Define types
    using RF = typename LFSV::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::RangeFieldType;
    using size_type = typename LFSV::Traits::SizeType;

    // Reference to the inside cell
    const auto& cell_inside = ig.inside();

    // Get geometry
    auto geo = ig.geometry();

    // Get geometry of intersection in local coordinates of cell_inside
    auto geo_in_inside = ig.geometryInInside();

    // evaluate boundary condition type
    auto ref_el = referenceElement(geo_in_inside);
    auto local_face_center = ref_el.position(0,0);
    auto intersection = ig.intersection();
    auto bctype = param.bctype(intersection,local_face_center);

    // skip rest if we are on Dirichlet boundary
    if (bctype==Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet) return;

    // loop over quadrature points and integrate normal flux
    auto intorder = intorderadd+2*lfsu_s.finiteElement().localBasis().order();
    for (const auto& ip : quadratureRule(geo,intorder))
      {
        // position of quadrature point in local coordinates of element
        auto local = geo_in_inside.global(ip.position());

        // evaluate shape functions (assume Galerkin method)
        auto& phi = cache.evaluateFunction(local,lfsu_s.finiteElement().localBasis());
        //auto& phitest = cache.evaluateFunction(local,lfsv_s.finiteElement().localBasis());

        RF u=0.0;
        for (size_type i=0; i<lfsu_s.size(); i++)
          u += x_s(lfsu_s,i)*phi[i];

        auto factor = ip.weight()*geo.integrationElement(ip.position());
        for (size_type i=0; i<lfsu_s.size(); i++)
          r_s.accumulate(lfsu_s,i,u*phi[i]*factor);

      }
  }


  //! set time in parameter class
  void setTime (double t)
  {
    param.setTime(t);
  }

private:
  T& param;
  int intorderadd;
  using LocalBasisType = typename FiniteElementMap::Traits::FiniteElementType::Traits::LocalBasisType;
  Dune::PDELab::LocalBasisCache<LocalBasisType> cache;
};




void driver(std::string basis_type, std::string part_unity_type) {
  using Dune::PDELab::Backend::native;

  Dune::Timer timer_assembly;


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
  std::array<int,dim> yasppartitions = {configuration.get<int>("subdomainsx"), configuration.get<int>("subdomainsy")};
  auto yp = new YP(yasppartitions);

  auto grid = std::make_shared<GM>(L,N,B,overlap,Dune::MPIHelper::getCollectiveCommunication(),yp);


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
  typedef DTNRhsLocalOperator<Problem,FEM> DTN_LOP;
  LOP lop(problem);
  DTN_LOP dtn_lop(problem);

  // LocalOperator wrapper zeroing out subdomains' interiors in order to set up overlap matrix
  typedef Dune::PDELab::LocalOperatorOvlpRegion<LOP, GFS> LOP_OVLP;
  LOP_OVLP lop_ovlp(lop, gfs);

  typedef Dune::PDELab::ISTL::BCRSMatrixBackend<> MBE;


  // Construct GridOperators from LocalOperators
  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,NumberType,NumberType,NumberType,CC,CC> GO;
  auto go = GO(gfs,cc,gfs,cc,lop,MBE(nonzeros));

  typedef Dune::PDELab::GridOperator<GFS,GFS,DTN_LOP,MBE,NumberType,NumberType,NumberType,CC,CC> DTN_GO;
  auto dtn_go = DTN_GO(gfs,cc,gfs,cc,dtn_lop,MBE(nonzeros));

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
  Dune::storeMatrixMarket(native(AF), std::to_string(gfs.gridView().comm().rank()) + "FineMatrix.mm");

  M M_dtn(dtn_go);
  dtn_go.jacobian(x,M_dtn);
  Dune::storeMatrixMarket(native(M_dtn), std::to_string(gfs.gridView().comm().rank()) + "DTNFineMatrix.mm");

  typedef Dune::PDELab::OverlappingOperator<CC,M,V,V> POP;
  auto popf = std::make_shared<POP>(cc,AF);
  typedef Dune::PDELab::ISTL::ParallelHelper<GFS> PIH;
  PIH pihf(gfs);
  typedef Dune::PDELab::OverlappingScalarProduct<GFS,V> OSP;
  OSP ospf(gfs,pihf);

  // Assemble fine grid matrix defined without processor constraints
  M_EXTERIOR AF_exterior(go_exterior);
  go_exterior.jacobian(x_exterior,AF_exterior);
  Dune::storeMatrixMarket(native(AF), std::to_string(gfs.gridView().comm().rank()) + "FineMatrixExterior.mm");

  // Assemble fine grid matrix defined only on overlap region
  M_EXTERIOR AF_ovlp(go_overlap);
  go_overlap.jacobian(x_exterior,AF_ovlp);
  Dune::storeMatrixMarket(native(AF_ovlp), std::to_string(gfs.gridView().comm().rank()) + "OverlapMatrix.mm");


  // Choose an eigenvalue threshold according to Spillane et al., 2014.
  // This particular choice is a heuristic working very well for Darcy problems.
  // Theoretically, any value delivers robustness; practically, you may need to
  // choose another value to achieve a good balance between condition bound and
  // global basis size
  double eigenvalue_threshold = (double)overlap / (cells + overlap);

  int verb=0;
  if (gfs.gridView().comm().rank()==0) verb=2;

  if (verb > 0)
    std::cout << "Running basis " << basis_type << " with " << part_unity_type << " partition of unity" << std::endl;

  // Create a local function space needed for the Sarkis partition of unity only
  typedef Dune::PDELab::LocalFunctionSpace<GFS, Dune::PDELab::AnySpaceTag> LFS;
  LFS lfs(gfs);

  MPI_Barrier (MPI_COMM_WORLD);
  if (verb > 0) std::cout << "Assembly: " << timer_assembly.elapsed() << std::endl;
  Dune::Timer timer_full_solve;
  Dune::Timer timer_part_unity;

  // Generate a partition of unity
  std::shared_ptr<V> part_unity;
  if (part_unity_type == "standard")
    part_unity = std::make_shared<V>(standardPartitionOfUnity<V>(gfs, cc_bnd_neu_int_dir));
  else if (part_unity_type == "sarkis")
    part_unity = std::make_shared<V>(sarkisPartitionOfUnity<V>(gfs, lfs, cc_bnd_neu_int_dir, cells, cells, overlap, yasppartitions[0], yasppartitions[1]));
  else
    DUNE_THROW(Dune::Exception, "Unkown selection in test driver!");

  /*M ovlp_mat(AF_ovlp);
  for (auto row_iter = native(ovlp_mat).begin(); row_iter != native(ovlp_mat).end(); row_iter++) {
    for (auto col_iter = row_iter->begin(); col_iter != row_iter->end(); col_iter++) {
      *col_iter *= native(*part_unity)[row_iter.index()] * native(*part_unity)[col_iter.index()];
    }
  }*/
  M_EXTERIOR ovlp_mat(AF_ovlp);
  for (auto row_iter = native(ovlp_mat).begin(); row_iter != native(ovlp_mat).end(); row_iter++) {
    for (auto col_iter = row_iter->begin(); col_iter != row_iter->end(); col_iter++) {
      *col_iter *= native(*part_unity)[row_iter.index()] * native(*part_unity)[col_iter.index()];
    }
  }
  Dune::storeMatrixMarket(native(ovlp_mat), std::to_string(gfs.gridView().comm().rank()) + "PartUnityOverlapMatrixPartUnity.mm");


  M_EXTERIOR ext_mat_part_unity(AF_exterior);
  for (auto row_iter = native(ext_mat_part_unity).begin(); row_iter != native(ext_mat_part_unity).end(); row_iter++) {
    for (auto col_iter = row_iter->begin(); col_iter != row_iter->end(); col_iter++) {
      *col_iter *= native(*part_unity)[row_iter.index()] * native(*part_unity)[col_iter.index()];
    }
  }
  Dune::storeMatrixMarket(native(ext_mat_part_unity), std::to_string(gfs.gridView().comm().rank()) + "PartUnityExteriorMatrixPartUnity.mm");



  MPI_Barrier (MPI_COMM_WORLD);
  if (verb > 0) std::cout << "Part unity: " << timer_part_unity.elapsed() << std::endl;

  // Choose how many eigenvalues to compute
  int nev = configuration.get<int>("nev");
  //if (basis_type == "rndgeneo") nev += 10;
  int nev_arpack = std::ceil(double(nev) * configuration.get<double>("arpack_factor")); // FIXME: This still correct? // configuration.get<int>("nev_arpack");
  eigenvalue_threshold = -1.0;

  Dune::Timer timer_basis;
  // Construct per-subdomain basis functions
  std::shared_ptr<Dune::PDELab::SubdomainBasis<V> > subdomain_basis;
  if (basis_type == "geneo")
    subdomain_basis = std::make_shared<Dune::PDELab::GenEOBasis<GFS,M_EXTERIOR,V,1> >(gfs, AF_exterior, AF_ovlp, eigenvalue_threshold, *part_unity, nev, nev_arpack, 0.001, false, verb, 0.0);
  else if (basis_type == "geneo_1e-3")
    subdomain_basis = std::make_shared<Dune::PDELab::GenEOBasis<GFS,M_EXTERIOR,V,1> >(gfs, AF_exterior, AF_ovlp, eigenvalue_threshold, *part_unity, nev, nev_arpack, 0.001, false, verb, 1e-3);
  else if (basis_type == "geneo_1e-6")
    subdomain_basis = std::make_shared<Dune::PDELab::GenEOBasis<GFS,M_EXTERIOR,V,1> >(gfs, AF_exterior, AF_ovlp, eigenvalue_threshold, *part_unity, nev, nev_arpack, 0.001, false, verb, 1e-6);
  else if (basis_type == "rndgeneo")
    subdomain_basis = std::make_shared<Dune::PDELab::RndGenEOBasis<GFS,M_EXTERIOR,V,1> >(gfs, AF_exterior, AF_ovlp, eigenvalue_threshold, *part_unity, nev, verb);
  else if (basis_type == "fastrndgeneo")
    subdomain_basis = std::make_shared<Dune::PDELab::FastRndGenEOBasis<GFS,M_EXTERIOR,V,1> >(gfs, AF_exterior, AF_ovlp, eigenvalue_threshold, *part_unity, nev, verb);
  else if (basis_type == "fastrndgeneoadaptive")
    subdomain_basis = std::make_shared<Dune::PDELab::FastRndGenEOBasisAdaptive<GFS,M_EXTERIOR,V,1> >(gfs, AF_exterior, AF_ovlp, eigenvalue_threshold, *part_unity, nev, verb);
  else if (basis_type == "fastrndgeneo2")
    subdomain_basis = std::make_shared<Dune::PDELab::FastRndGenEOBasis<GFS,M_EXTERIOR,V,1> >(gfs, AF_exterior, AF_ovlp, eigenvalue_threshold, *part_unity, nev, verb, 2);
  else if (basis_type == "lipton_babuska")
    subdomain_basis = std::make_shared<Dune::PDELab::LiptonBabuskaBasis<GFS,M_EXTERIOR,V,V,1> >(gfs, AF_exterior, AF_ovlp, -1, *part_unity, nev, nev_arpack);
  else if (basis_type == "part_unity") // We can't test this one, it does not lead to sufficient error reduction. Let's instantiate it anyway for test's sake.
    subdomain_basis = std::make_shared<Dune::PDELab::SubdomainBasis<V> >(*part_unity);
  else if (basis_type == "onelevel")
    subdomain_basis = std::make_shared<Dune::PDELab::SubdomainBasis<V> >(*part_unity);
  else
    DUNE_THROW(Dune::Exception, "Unkown selection in test driver!");
  MPI_Barrier (MPI_COMM_WORLD);
  if (verb > 0) std::cout << "Basis setup: " << timer_basis.elapsed() << std::endl;


  // Fuse per-subdomain basis functions to a global coarse space
  auto coarse_space = std::make_shared<Dune::PDELab::SubdomainProjectedCoarseSpace<GFS,M_EXTERIOR,V,PIH> >(gfs, AF_exterior, subdomain_basis, pihf);


  // Plug coarse basis into actual preconditioner
  auto prec = std::make_shared<Dune::PDELab::ISTL::TwoLevelOverlappingAdditiveSchwarz<GFS,M,M_EXTERIOR,V,V>>(gfs, AF, AF_exterior, coarse_space, *part_unity, false, basis_type != "onelevel", configuration.get<bool>("hybrid"), verb);

  //Dune::Richardson<V,V> richardson();

  Dune::Timer timer_solve;

  // now solve defect equation A*v = d using a CG solver with our shiny preconditioner
  V v(gfs,0.0);
  auto solver_ref = std::make_shared<Dune::CGSolver<V> >(*popf,ospf,*prec,1E-6,1000,verb,false);
  Dune::InverseOperatorResult result;
  solver_ref->apply(v,d,result);
  x -= v;

  MPI_Barrier (MPI_COMM_WORLD);
  if (verb > 0) std::cout << "pCG solve: " << timer_solve.elapsed() << std::endl;
  if (verb > 0) std::cout << "full solve: " << timer_full_solve.elapsed() << std::endl;

  // Write solution to VTK
  Dune::VTKWriter<GV> vtkwriter(gfs.gridView());
  typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
  DGF xdgf(gfs,x);
  typedef Dune::PDELab::VTKGridFunctionAdapter<DGF> ADAPT;
  auto adapt = std::make_shared<ADAPT>(xdgf,"solution");
  vtkwriter.addVertexData(adapt);
  vtkwriter.write("testgeneo_basis_" + basis_type + "_part_unity_" + part_unity_type);


  if  (configuration.get<bool>("write_basis", false)) {
    for (int i = 0; i < subdomain_basis->basis_size(); i++) {
      std::shared_ptr<V> vec = subdomain_basis->get_basis_vector(i);

      Dune::VTKWriter<GV> vtkwriter(gfs.gridView());
      typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
      DGF xdgf(gfs,*vec);
      typedef Dune::PDELab::VTKGridFunctionAdapter<DGF> ADAPT;
      auto adapt = std::make_shared<ADAPT>(xdgf,"basis");
      vtkwriter.addVertexData(adapt);
      vtkwriter.write("testgeneo_basis_" + basis_type + "_part_unity_" + part_unity_type
                      + "_r_" + std::to_string(gfs.gridView().comm().rank()) + "_vec_" + std::to_string(i));

    }
  }


}

void error_handler(int status, const char *file,
 int line, const char *message) {
std::cout<<"!!"<<std::endl;
}

int main(int argc, char **argv)
{
  Dune::ParameterTreeParser parser;
  parser.readINITree("configmatrixwriter.ini", configuration);
  parser.readOptions(argc, argv, configuration);

  using Dune::PDELab::Backend::native;

  try{
    // initialize MPI, finalize is done automatically on exit
    Dune::MPIHelper::instance(argc,argv);

    //std::cout << "#################### geneo" << std::endl;
    //driver("geneo", "standard");
    //std::cout << "#################### rndgeneo" << std::endl;
    driver(configuration.get<std::string>("method"), configuration.get<std::string>("part_unity"));
    //std::cout << "#################### babuska" << std::endl;
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
