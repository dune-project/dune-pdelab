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
class GenericEllipticProblem2
{
  typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;

public:
  typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

  GenericEllipticProblem2()
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


    int layer = (int)std::floor(xglobal[1]/ layer_thickness);
    bool low_perm = layer % 2 == 0;
    RF coeff;

    //if (layer % 3 == 0)
      coeff = low_perm ? perm1 : perm2;
    /*else if (layer % 3 == 1)
      coeff = low_perm && xglobal[0] > .1 ? perm1 : perm2;
    else if (layer % 3 == 2)
      coeff = low_perm && xglobal[0] < .9 ? perm1 : perm2;*/

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
    //return 1.0;//std::exp(-std::sqrt(std::pow(xglobal[0] - .5, 2) + std::pow(xglobal[1] - .5, 2)));
    return (- std::exp(-std::sqrt(std::pow(xglobal[0] - 0.25, 2) + std::pow(xglobal[1] - 0.25, 2)))
           + std::exp(-std::sqrt(std::pow(xglobal[0] - 0.75, 2) + std::pow(xglobal[1] - 0.75, 2)))) * 10;
    //return 1.0;
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
    /*typename Traits::DomainType xglobal = e.geometry().global(x);

    RF perm1 = 1e0;
    RF perm2 = contrast;
    RF layer_thickness = 1.0 / (double)layers;


    int layer = (int)std::floor(xglobal[1]/ layer_thickness);
    bool low_perm = layer % 2 == 0;
    RF coeff;

    if (layer % 3 == 0)
      coeff = low_perm ? perm1 : perm2;
    else if (layer % 3 == 1)
      coeff = low_perm && xglobal[0] > .1 ? perm1 : perm2;
    else if (layer % 3 == 2)
      coeff = low_perm && xglobal[0] < .9 ? perm1 : perm2;

    typename Traits::PermTensorType I;
    I[0][0] = coeff;
    I[0][1] = 0;
    I[1][0] = 0;
    I[1][1] = coeff;
    return I;*/

    typename Traits::DomainType xglobal = e.geometry().global(x);
    int num1 = floor(8*xglobal[0]);
    int num2 = floor(8*xglobal[1]);
    RF coeff = 0.0;
    if ( (num1 % 2 == 0) && (num2 % 2 == 0) )
      coeff = 0.1*contrast*(num2+1.0);
    else
      coeff = 1.0;

    RF duct1_b = 0.9*xglobal[0];
    RF duct1_t = duct1_b + 0.1;
    if ( (xglobal[1]>duct1_b) && (xglobal[1]<duct1_t))
      coeff = coeff + 0.2 * contrast;

    RF duct2_b = -1.0*xglobal[0] + 0.5;
    RF duct2_t = duct2_b + 0.1;
    if ( (xglobal[1]>duct2_b) && (xglobal[1]<duct2_t))
      coeff = coeff + 0.5 * contrast;

    RF duct3_b = 4.0*xglobal[0] - 2.0 ;
    RF duct3_t = duct3_b + 0.1;
    if ( (xglobal[1]>duct3_b) && (xglobal[1]<duct3_t))
      coeff = coeff + 0.2 * contrast;

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
    //return 1.0;//std::exp(-std::sqrt(std::pow(xglobal[0] - .5, 2) + std::pow(xglobal[1] - .5, 2)));
    //return (- std::exp(-std::sqrt(std::pow(xglobal[0] - 0.25, 2) + std::pow(xglobal[1] - 0.25, 2)))
    //+ std::exp(-std::sqrt(std::pow(xglobal[0] - 0.75, 2) + std::pow(xglobal[1] - 0.75, 2)))) * 100;
    return 0.0;
  }

  BCType
  bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    typename Traits::DomainType xglobal = is.geometry().global(x);
    /*if (!(xglobal[0]<1E-10 || xglobal[0]>1.0-1E-10))
     *    return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
     *  else
     *    return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;*/
    /*if (!((xglobal[1]<0.5 && xglobal[0]>1.0-1E-10) || (xglobal[1]>0.50 && xglobal[0]<1E-10)))
     *    return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
     *  else
     *    return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;*/
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
     *    return std::cos(xglobal[0] * 8* 3.141592);
     *  else if (xglobal[1] > 1.0-1E-10 && xglobal[0] < 0.25)
     *    return std::sin(xglobal[0] * 8 * 3.141592);
     *  else if (xglobal[0] < 1E-10)
     *    return std::sin(xglobal[1] * 2 * 3.141592) + 1.0;
     *  else if (xglobal[0] > 1.0 - 1E-10)
     *    return std::cos(xglobal[1] * 2 * 3.141592) - 1.0;
     *  else
     *    return 0.0;*/
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

void driver(std::string basis_type, std::string part_unity_type) {

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
  V x_orig(x);

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
  else if (basis_type == "geneo_1e-1")
    subdomain_basis = std::make_shared<Dune::PDELab::GenEOBasis<GFS,M_EXTERIOR,V,1> >(gfs, AF_exterior, AF_ovlp, eigenvalue_threshold, *part_unity, nev, nev_arpack, 0.001, false, verb, 1e-1);
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
  int basis_size_without_augmentation = subdomain_basis->local_basis.size();
  if (configuration.get<bool>("coarse_only") && configuration.get<bool>("dirichlet_interpol_in_coarse_only")) {
    //subdomain_basis->local_basis.push_back(Dune::stackobject_to_shared_ptr(x_orig));
    using Dune::PDELab::Backend::native;
    bool nonzero = false;
    for (auto& val : native(x_orig)) {
      if (val != .0) {
        nonzero = true;
        break;
      }
    }
    if (nonzero) {
      subdomain_basis->local_basis.push_back(std::make_shared<V>(x_orig));
    }
  }
  auto coarse_space = std::make_shared<Dune::PDELab::SubdomainProjectedCoarseSpace<GFS,M_EXTERIOR,V,PIH> >(gfs, AF_exterior, subdomain_basis, pihf);

  // Plug coarse basis into actual preconditioner
  bool restricted = configuration.get<bool>("hybrid"); // NOTE: hybrid and restricted mode coupled here!!
  auto prec = std::make_shared<Dune::PDELab::ISTL::TwoLevelOverlappingAdditiveSchwarz<GFS,M,M_EXTERIOR,V,V>>(gfs, AF, AF_exterior, coarse_space, *part_unity, restricted, basis_type != "onelevel", configuration.get<bool>("hybrid"), verb);

  std::shared_ptr<Dune::PDELab::GenEOBasis<GFS,M_EXTERIOR,V,1> > subdomain_basis_full;
  std::shared_ptr<Dune::PDELab::SubdomainProjectedCoarseSpace<GFS,M_EXTERIOR,V,PIH> > coarse_space_full;
  std::shared_ptr<Dune::PDELab::ISTL::TwoLevelOverlappingAdditiveSchwarz<GFS,M,M_EXTERIOR,V,V> > prec_full;
  if (configuration.get<bool>("coarse_only")) {
    subdomain_basis_full = std::make_shared<Dune::PDELab::GenEOBasis<GFS,M_EXTERIOR,V,1> >(gfs, AF_exterior, AF_ovlp, eigenvalue_threshold, *part_unity, nev, nev_arpack, 0.001, false, verb, 0.0);
    coarse_space_full = std::make_shared<Dune::PDELab::SubdomainProjectedCoarseSpace<GFS,M_EXTERIOR,V,PIH> >(gfs, AF_exterior, subdomain_basis_full, pihf);
    prec_full = std::make_shared<Dune::PDELab::ISTL::TwoLevelOverlappingAdditiveSchwarz<GFS,M,M_EXTERIOR,V,V>>(gfs, AF, AF_exterior, coarse_space_full, *part_unity, false, basis_type != "onelevel", configuration.get<bool>("hybrid"), verb);
  }
  //auto prec = std::make_shared<Dune::PDELab::ISTL::TwoLevelOverlappingAdditiveSchwarz<GFS,M,M_EXTERIOR,V,V>>(gfs, AF, AF_exterior, coarse_space, *part_unity, false, basis_type != "onelevel", configuration.get<bool>("hybrid"), verb);

  //Dune::Richardson<V,V> richardson();

while(true) {
  V x_orig2 = x_orig;
  V x2 = x;

  Dune::Timer timer_solve;

  // now solve defect equation A*v = d using a CG solver with our shiny preconditioner
  V v(gfs,0.0);
  if (configuration.get<bool>("coarse_only")) {
    typedef Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > COARSE_M;
    Dune::UMFPack<COARSE_M> coarse_solver(*coarse_space->get_coarse_system());



    typedef Dune::BlockVector<Dune::FieldVector<double,1> > COARSE_V;
    COARSE_V coarse_defect(coarse_space->basis_size(),coarse_space->basis_size());
    V d2 = d;
    coarse_space->restrict (d2, coarse_defect);
    // Solve coarse system
    COARSE_V v0(coarse_space->basis_size(),coarse_space->basis_size());

    Dune::InverseOperatorResult result;
    coarse_solver.apply(v0, coarse_defect, result);
    // Prolongate coarse solution on local domain
    V prolongated(gfs, 0.0);
    coarse_space->prolongate(v0, prolongated);

    Dune::PDELab::AddDataHandle<GFS,V> prolongated_addh(gfs,prolongated);
    gfs.gridView().communicate(prolongated_addh,Dune::All_All_Interface,Dune::ForwardCommunication);

    x2 -= prolongated;
  } else {
    if (configuration.get<bool>("hybrid")) {
      auto solver_ref = std::make_shared<Dune::RestartedGMResSolver<V> >(*popf,ospf,*prec,1E-6,30,1000,verb);
      Dune::InverseOperatorResult result;
      solver_ref->apply(v,d,result);
      x2 -= v;
    } else {
      auto solver_ref = std::make_shared<Dune::CGSolver<V> >(*popf,ospf,*prec,1E-6,1000,verb,false);
      Dune::InverseOperatorResult result;
      solver_ref->apply(v,d,result);
      x2 -= v;
    }
  }

  MPI_Barrier (MPI_COMM_WORLD);
  if (verb > 0) std::cout << "pCG solve: " << timer_solve.elapsed() << std::endl;
  if (verb > 0) std::cout << "full solve: " << timer_full_solve.elapsed() << std::endl;

  std::string vtk_prefix = "testgeneo_basis_" + basis_type + "_part_unity_" + part_unity_type + "_dir_" + (configuration.get<bool>("dirichlet_interpol_in_coarse_only") ? "true" : "false");

  if (configuration.get<bool>("coarse_only")) {
    auto solver_ref = std::make_shared<Dune::CGSolver<V> >(*popf,ospf,*prec_full,1E-6,1000,verb,false);
    Dune::InverseOperatorResult result;
    V d2 = d;
    solver_ref->apply(v,d2,result);
    x_orig2 -= v;

    typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
    DGF dgf_ref(gfs, x_orig2);
    DGF dgf_soln(gfs, x2);

    {
    typedef Dune::PDELab::DifferenceAdapter<DGF,DGF> Difference;
    Difference error_diff(dgf_ref,dgf_soln);
    Dune::VTKWriter<GV> vtkwriter(gfs.gridView());
    typedef Dune::PDELab::VTKGridFunctionAdapter<Difference> ADAPT;
    auto adapt = std::make_shared<ADAPT>(error_diff,"error");
    vtkwriter.addVertexData(adapt);
    vtkwriter.write(vtk_prefix + "_error" + std::to_string(basis_size_without_augmentation));
    }

    typedef Dune::PDELab::DifferenceSquaredAdapter<DGF,DGF> DifferenceSquared;
    DifferenceSquared differencesquared(dgf_ref,dgf_soln);
    typename DifferenceSquared::Traits::RangeType error_norm_squared(0.0);
    Dune::PDELab::integrateGridFunction (differencesquared, error_norm_squared, 10);

    error_norm_squared = gfs.gridView().comm().sum(error_norm_squared);

    if (gfs.gridView().comm().rank() == 0)
      std::cout << "Approx error norm AE: " << std::sqrt(error_norm_squared) << std::endl;
  }



  // Write solution to VTK
  Dune::VTKWriter<GV> vtkwriter(gfs.gridView());
  typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
  DGF xdgf(gfs,x2);
  typedef Dune::PDELab::VTKGridFunctionAdapter<DGF> ADAPT;
  auto adapt = std::make_shared<ADAPT>(xdgf,"solution");
  vtkwriter.addVertexData(adapt);
  DGF xdgf_orig(gfs,x_orig2);
  if (configuration.get<bool>("coarse_only")) {
    auto adapt_orig = std::make_shared<ADAPT>(xdgf_orig,"ref");
    vtkwriter.addVertexData(adapt_orig);
  }
  vtkwriter.write(vtk_prefix);

  if  (configuration.get<bool>("write_basis", false)) {
    for (int i = 0; i < subdomain_basis->basis_size(); i++) {
      std::shared_ptr<V> vec = subdomain_basis->get_basis_vector(i);

      Dune::VTKWriter<GV> vtkwriter(gfs.gridView());
      typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
      DGF xdgf(gfs,*vec);
      typedef Dune::PDELab::VTKGridFunctionAdapter<DGF> ADAPT;
      auto adapt = std::make_shared<ADAPT>(xdgf,"basis");
      vtkwriter.addVertexData(adapt);
      vtkwriter.write(vtk_prefix
                      + "_r_" + std::to_string(gfs.gridView().comm().rank()) + "_vec_" + std::to_string(i));

    }
  }

  if (!configuration.get<bool>("coarse_reduction_test"))
    break;
  if (subdomain_basis->local_basis.size() == 1)
    break;
  subdomain_basis->local_basis.pop_back();
  coarse_space = std::make_shared<Dune::PDELab::SubdomainProjectedCoarseSpace<GFS,M_EXTERIOR,V,PIH> >(gfs, AF_exterior, subdomain_basis, pihf);
  prec = std::make_shared<Dune::PDELab::ISTL::TwoLevelOverlappingAdditiveSchwarz<GFS,M,M_EXTERIOR,V,V>>(gfs, AF, AF_exterior, coarse_space, *part_unity, false, basis_type != "onelevel", configuration.get<bool>("hybrid"), verb);
}


}

void error_handler(int status, const char *file,
 int line, const char *message) {
std::cout<<"!!"<<std::endl;
}

int main(int argc, char **argv)
{
  Dune::ParameterTreeParser parser;
  parser.readINITree("config.ini", configuration);
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
