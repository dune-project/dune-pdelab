#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/pdelab.hh>
#include <dune/pdelab/backend/istl/geneo/OfflineOnline/geneobasisOnline.hh>
#include <dune/pdelab/backend/istl/geneo/OfflineOnline/OnlineTools.hh>
#include <dune/pdelab/backend/istl/geneo/OfflineOnline/SubDomainGmshReader.hh>

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
    RF perm2 = 1e5; // FIXME we want high contrast
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

void driver(std::string path_to_storage, std::vector<int> targeted, Dune::MPIHelper& helper) {

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

  // Dune::gridlevellist<GRID>(*grid, 0, "");

  // ~~~~~~~~~~~~~~~~~~
//  Type definitions
  // ~~~~~~~~~~~~~~~~~~
  const int components = 1;
  using K = double;
  using Vector = Dune::BlockVector<Dune::FieldVector<K,components>>;
  using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<K,components,components>>;

  typedef Dune::BlockVector<Dune::FieldVector<K, 1>> CoarseVector;
  typedef Dune::BCRSMatrix<Dune::FieldMatrix<K, 1, 1>> CoarseMatrix;

  typedef Dune::BlockVector<Dune::FieldVector<int, 1>> vector1i;

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

  // ~~~~~~~~~~~~~~~~~~
//  Solving process begin here: First some parameters
  // ~~~~~~~~~~~~~~~~~~
  double eigenvalue_threshold = -1;
  // const int algebraic_overlap = 0;
  int nev = 4;
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
  vector1i lb = FromOffline<vector1i>(path_to_storage, "localBasisSizes");
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
  vector1i offlineDoF2GI = FromOffline<vector1i>(path_to_storage, "GI", targeted[0]);
  std::vector<int> DofOffline_to_DofOnline = offlineDoF2GI2gmsh2onlineDoF<vector1i>(targeted[0], gmsh2dune, offlineDoF2GI, path_to_storage);

  // ~~~~~~~~~~~~~~~~~~
  // Load PoU
  // ~~~~~~~~~~~~~~~~~~
  Vector PoU = FromOffline<Vector>(path_to_storage, "PoU", targeted[0]);
  V newPoU(gfs, 0.0);
  Vector nPoU(v_size);
  // int cnt=0;
  for (int i=0; i<v_size; i++){
    native(newPoU)[DofOffline_to_DofOnline[i]] = PoU[i];
    nPoU[DofOffline_to_DofOnline[i]] = PoU[i];
    // cnt+=1;
  }

  /* Write a field in vtu */
  Dune::VTKWriter<GV> vtkwriter(gv);
  DGF xdgf(gfs,newPoU);
  auto adapt = std::make_shared<ADAPT>(xdgf,"PoU");
  vtkwriter.addVertexData(adapt);
  vtkwriter.write("PoU",Dune::VTK::ascii);

  // ~~~~~~~~~~~~~~~~~~
//  Subdomain basis computation or loading
  // ~~~~~~~~~~~~~~~~~~

  int basis_size = local_basis_sizes[targeted[0]];

  std::shared_ptr<Dune::PDELab::SubdomainBasis<Vector>> online_subdomainbasis;
  online_subdomainbasis = std::make_shared<Dune::PDELab::GenEOBasisOnline<GO, Matrix, Vector>>(native(A), nPoU, eigenvalue_threshold, nev, nev_arpack);


  for (int i=0; i<basis_size; i++){
    /* Plot EV */
    Dune::VTKWriter<GV> vtkwriterEV(gv);
    V EV(gfs,0.0);
    native(EV) = *online_subdomainbasis->get_basis_vector(i);
    // for (int i=0; i<v_size; i++)
    //   std::cout << i << " :: " << native(EV)[i] << std::endl;
    /* Write a field in the vtu file */
    DGF xdgfEV(gfs,EV);
    auto adaptEV = std::make_shared<ADAPT>(xdgfEV,"EV");
    vtkwriterEV.addVertexData(adaptEV);
    vtkwriterEV.write("EV-"+std::to_string(i),Dune::VTK::ascii);
  }

  std::shared_ptr<Dune::PDELab::SubdomainBasis<Vector>> offline_subdomainbasis;
  offline_subdomainbasis = std::make_shared<Dune::PDELab::GenEOBasisFromFiles<GO, Matrix, Vector>>(path_to_storage, basis_size, targeted[0]);

  for (int i=0; i<basis_size; i++){
    /* Plot EV */
    Dune::VTKWriter<GV> vtkwriterEV(gv);
    V EV(gfs,0.0);

    for (int j=0; j<v_size; j++){
      native(EV)[DofOffline_to_DofOnline[j]] = (*offline_subdomainbasis->get_basis_vector(i))[j];
    }
    /* Write a field in the vtu file */
    DGF xdgfEV(gfs,EV);
    auto adaptEV = std::make_shared<ADAPT>(xdgfEV,"EV");
    vtkwriterEV.addVertexData(adaptEV);
    vtkwriterEV.write("offEV-"+std::to_string(i),Dune::VTK::ascii);
  }

}


int main(int argc, char **argv)
{
  using Dune::PDELab::Backend::native;

  // Offline folder
  std::string path_to_storage = "Offline/";

  // Define what subdomain need to be solved
  std::vector<int> targeted = {1}; // Subdomains that need a second solve

  // initialize MPI, finalize is done automatically on exit
  Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc,argv);

  driver(path_to_storage, targeted, helper);

  return 0;
}