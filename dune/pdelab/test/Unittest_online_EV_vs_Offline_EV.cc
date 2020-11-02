#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/pdelab.hh>
#include <dune/pdelab/backend/istl/geneo/nonoverlapping/geneobasisOnline.hh>
#include <dune/pdelab/backend/istl/geneo/nonoverlapping/geneobasisfromfiles.hh>

#include <dune/grid/common/gridinfo.hh> // visualize grid information -> could be remove at the end

// #include <dune/pdelab/test/testordering.cc>
#include <dune/pdelab/test/SubDomainGmshReader.hh>

template<typename GFS>
void check_ordering(const GFS& gfs)
{
    const typename GFS::Ordering& ordering = gfs.ordering();

    Dune::PDELab::LocalFunctionSpace<GFS> lfs(gfs);

    typedef typename GFS::Traits::GridView GV;

    for (typename GV::template Codim<0>::Iterator it = gfs.gridView().template begin<0>();
         it != gfs.gridView().template end<0>(); ++it)
    {
      lfs.bind(*it);

      std::vector<typename GFS::Ordering::Traits::DOFIndex> vdi(lfs.size());
      std::vector<typename GFS::Ordering::Traits::ContainerIndex> vci(lfs.size());
      std::array<std::size_t,Dune::TypeTree::TreeInfo<GFS>::leafCount> leaf_sizes;
      std::fill(begin(leaf_sizes),end(leaf_sizes),0);
      for (unsigned i = 0; i < lfs.size(); ++i)
      {
        vdi[i] = lfs.dofIndex(i);
      }
      Dune::PDELab::map_dof_indices_to_container_indices<
        typename decltype(vdi)::iterator,
        typename decltype(vci)::iterator,
        typename decltype(leaf_sizes)::iterator,
        Dune::TypeTree::TreeInfo<GFS>::depth,
        false
        > visitor(begin(vdi),begin(vci),begin(leaf_sizes));

      Dune::TypeTree::applyToTree(ordering,visitor);

      for (unsigned i = 0; i < lfs.size(); ++i)
      {
        const typename GFS::Ordering::Traits::DOFIndex& di = lfs.dofIndex(i);
        typename GFS::Ordering::Traits::ContainerIndex ci;
        ordering.mapIndex(di.view(),ci);
        std::cout << di << "    " << vci[i] << "    " << ci << std::endl;
      }
    }

    using V = Dune::PDELab::Backend::Vector<GFS,double>;
    V x(gfs);
    x = 0.0;
    std::cout << std::endl;
}

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

void driver(std::string basis_type, std::string part_unity_type, Dune::MPIHelper& helper) {

  // ~~~~~~~~~~~~~~~~~~
  // ReCreate offline hierarchy
  // ~~~~~~~~~~~~~~~~~~
  std::string path_to_storage = "Offline/";

  // ~~~~~~~~~~~~~~~~~~
  // Define what subdomain need to be solved
  // ~~~~~~~~~~~~~~~~~~
  std::vector<int> targeted = {0}; // Subdomains that need a second solve

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
  V d(gfs,0.0);
  go.residual(x,d); // The rhs is loaded from offline directly restricted in the coarse spaces

  // ~~~~~~~~~~~~~~~~~~
//  Solving process begin here: First some parameters
  // ~~~~~~~~~~~~~~~~~~
  double eigenvalue_threshold = -1;
  // const int algebraic_overlap = 0;
  int nev = 3;
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
  Vector lb;
  std::ifstream file_lb;
  std::string filename_lb = path_to_storage + "localBasisSizes.mm";
  file_lb.open(filename_lb.c_str(), std::ios::in);
  Dune::readMatrixMarket(lb,file_lb);
  file_lb.close();

  const int number_of_rank_used_offline = lb.size();

  std::vector<int> local_basis_sizes(number_of_rank_used_offline), local_offset(number_of_rank_used_offline+1);
  local_offset[0]=0;
  for (int i=0; i<lb.size();i++) {
    local_basis_sizes[i] = lb[i];
    local_offset[i+1] = lb[i]+local_offset[i];
  }

  // ~~~~~~~~~~~~~~~~~~
  // Load PoU
  // ~~~~~~~~~~~~~~~~~~
  Vector PoU(A.N());
  std::string filename_PoU = path_to_storage + std::to_string(targeted[0]) + "_PoU.mm";
  std::ifstream file_PoU;
  file_PoU.open(filename_PoU.c_str(), std::ios::in);
  Dune::readMatrixMarket(PoU,file_PoU);
  file_PoU.close();

  Vector online2GI(v_size);
  std::string filename_on2GI = path_to_storage + std::to_string(targeted[0]) + "_LocalToGlobalNode.mm";
  std::ifstream file_on2GI;
  file_on2GI.open(filename_on2GI.c_str(), std::ios::in);
  Dune::readMatrixMarket(online2GI,file_on2GI);
  file_on2GI.close();

  Vector offline2GI(v_size);
  std::string filename_off2GI = path_to_storage + std::to_string(targeted[0]) + "_GI.mm";
  std::ifstream file_off2GI;
  file_off2GI.open(filename_off2GI.c_str(), std::ios::in);
  Dune::readMatrixMarket(offline2GI,file_off2GI);
  file_off2GI.close();

  // /* First method to deal with online2GI and offline2GI */
  // std::vector<int> GI2gmsh;
  // auto result = std::max_element(online2GI.begin(), online2GI.end());
  // GI2gmsh.resize(online2GI[std::distance(online2GI.begin(), result)]);
  // for(int i=0; i<GI2gmsh.size(); i++)
  //   GI2gmsh[i]=-1;
  // for(int i=0; i<v_size; i++){
  //   GI2gmsh[online2GI[i]] = i;
  // }
  // /* End first method */

  /* Second method to deal with online2GI and offline2GI */
  std::vector<std::pair<int, int>> online(v_size), offline(v_size);
  for(int i=0; i<v_size; i++){
    online[i] = std::make_pair(online2GI[i],i);
    offline[i] = std::make_pair(offline2GI[i],i);
  }
  std::sort(online.begin(), online.end());
  std::sort(offline.begin(), offline.end());
  /* End second method */

  std::vector<int> DofOffline_to_DofOnline(v_size);
  // Change order from Dof offline to offline global
  for(int i=0; i<v_size; i++){ // go through dof offline ordering
    // DofOffline_to_DofOnline[i] = GI2gmsh[offline2GI[i]];
    DofOffline_to_DofOnline[offline[i].second] = online[i].second;
  }

  auto tmp = DofOffline_to_DofOnline;
  // We know that a re-ordering of nodes is done during the GmshReader operation
  for(int i=0; i<v_size; i++){ // go through dof offline ordering
    DofOffline_to_DofOnline[i] = gmsh2dune[tmp[i]];
  }

  // /* Visualise changes in terminal */
  // for(int i=0; i<v_size; i++){ // go through dof offline ordering
  //   std::cout << i << " -> " << DofOffline_to_DofOnline[i] << std::endl;
  // }

  V newPoU(gfs, 0.0);
  Vector nPoU(v_size);
  // int cnt=0;
  for (int i=0; i<v_size; i++){
    native(newPoU)[DofOffline_to_DofOnline[i]] = PoU[i];
    nPoU[DofOffline_to_DofOnline[i]] = PoU[i];
    // cnt+=1;
  }

  /* Write a field in the vtu file */
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

  // initialize MPI, finalize is done automatically on exit
  Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc,argv);

  driver("geneo", "standard", helper);

  return 0;
}


/////////////////////// KEEP ALL PREVIOUS TESTS FOR NOW

  // int cntE=0;
  // std::vector<int> nodes_from_elements;
  // nodes_from_elements.push_back(0);
  // for (const auto& e : elements(gv)){
  //   for (unsigned int i=0; i < e.subEntities(dim); i++){
  //     std::cout << " " << i << ":" << grid->levelIndexSet(0).subIndex(e,i,dim) << std::endl;
  //     auto it = std::find(nodes_from_elements.begin(), nodes_from_elements.end(), grid->levelIndexSet(0).subIndex(e,i,dim));
  //     if (it == nodes_from_elements.end()){
  //       nodes_from_elements.push_back(grid->levelIndexSet(0).subIndex(e,i,dim));
  //     }
  //   }
  //   cntE+=1;
  // }


  // std::cout << "v_size : " << v_size << std::endl;
  // std::cout << "this vector size : " << nodes_from_elements.size() << std::endl;
  // for(int i=0; i<nodes_from_elements.size(); i++)
  //   std::cout << nodes_from_elements[i] << std::endl;

  //   // // Create a vector from GFS ordering to node ordering
  // std::vector<std::pair<int, int>> element_nodes_ordering(v_size);

  // int cnt=0;
  // for(int i=0; i<nodes_from_elements.size(); i++) {
  //   element_nodes_ordering[cnt] = std::make_pair(nodes_from_elements[i],cnt);
  //   // std::cout << element_nodes_ordering[cnt].first << " : " << element_nodes_ordering[cnt].second << std::endl;
  //   cnt+=1;
  // }
  // std::sort(element_nodes_ordering.begin(), element_nodes_ordering.end());

  // V newPoU(gfs, 0.0);
  // cnt=0;
  // // for (const auto& vertex : vertices(gv)){
  // for(int i=0; i<v_size; i++){
  //   std::cout << cnt << " // " << element_nodes_ordering[cnt].second << std::endl;
  //   // std::cout << i << " // " << std::endl;

  //   native(newPoU)[cnt] = tmp[element_nodes_ordering[cnt].second];
  //   // native(newPoU)[cnt] = tmp[gfs2node[cnt].second];
  //   // native(newPoU)[cnt] = tmp[cnt];
  //   // std::cout << cnt << " // " << native(newPoU)[cnt] << std::endl;
  //   cnt+=1;
  // }

  // Dune::VTKWriter<GV> vtkwriter(gv);
  // typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
  // DGF xdgf(gfs,newPoU);
  // typedef Dune::PDELab::VTKGridFunctionAdapter<DGF> ADAPT;
  // auto adapt = std::make_shared<ADAPT>(xdgf,"solution");
  // vtkwriter.addVertexData(adapt);
  // vtkwriter.write("solution");

/////////////////////


  // check_ordering(gfs);
  // auto ordering = gfs.ordering();

  // // auto map_lfs_indices = ordering.map_lfs_indices();

  // std::cout << ordering.blockCount() << std::endl; // size of the ordering
  // // for(int i=0; i<ordering.size(); i++){
  // //   std::cout << ordering[i] << " : " << i << std::endl;
  // // }

  // const typename GFS::Ordering& ordering = gfs.ordering();
  // Dune::PDELab::LocalFunctionSpace<GFS> lfs(gfs);

/////////////////////


  // using LFS = Dune::PDELab::LocalFunctionSpace<GFS>;
  // using LFSCache = Dune::PDELab::LFSIndexCache<LFS>;

  // std::shared_ptr<V> data;
  // std::shared_ptr<const GFS> pgfs;
  // Dune::PDELab::LocalFunctionSpace<GFS> lfs(gfs);
  // LFSCache lfs_cache(lfs);

  // int iter=0;

  // for (const auto& e : elements(gv)){
  //   std::vector<int> nodes_from_elements;
  //   for (unsigned int i=0; i < e.subEntities(dim); i++){
  //     // std::cout << " " << i << ":" << grid->levelIndexSet(0).subIndex(e,i,dim) << std::endl;
  //     // auto it = std::find(nodes_from_elements.begin(), nodes_from_elements.end(), grid->levelIndexSet(0).subIndex(e,i,dim));
  //     // if (it == nodes_from_elements.end()){
  //     nodes_from_elements.push_back(grid->levelIndexSet(0).subIndex(e,i,dim));
  //     // }
  //   }

  //   // auto gt = e.type();
  //   // auto index = gfs.entitySet().indexSet().index(e);
  //   // std::cout << index << std::endl;
  //   // for (auto it = nodes_from_elements.begin(); it != nodes_from_elements.end(); ++it) {

  //   //   GFS::Ordering::Traits::DOFIndexAccessor::store(*it,gt,index,0);

  //   // }

  //   typename V::template LocalView<LFSCache> p_view(*data);
  //   lfs.bind(e);
  //   lfs_cache.update();
  //   // std::cout << iter << std::endl;

  //   std::vector<int> pw(lfs.size()-1);
  //   p_view.bind(lfs_cache);
  //   p_view.read(pw);

  //   // std::cout << lfs_cache.dofIndex(0) << ", " << lfs_cache.dofIndex(1) << ", " << lfs_cache.dofIndex(2) << ", " << lfs_cache.dofIndex(3) << std::endl;
  //   // for (int i=0; i<4; i++) {
  //   //   assert(lfs_cache.dofIndex(i).entityIndex()[1]==nodes_from_elements[i]);
  //   //   std::cout << lfs_cache.dofIndex(i).entityIndex()[1] << ", " << nodes_from_elements[i] << std::endl;
  //   // }
  //   // auto test = lfs_cache.dofIndex(0).entityIndex()[1];


  //   iter+=1;
  // }

/////////////////////

  // // Create a vector from GFS ordering to node ordering
  // std::vector<std::pair<int, int>> gfs2node(v_size);

  // cnt=0;
  // for (const auto& vertex : vertices(gv)){
  //   assert(cnt<v_size);
  //   gfs2node[cnt] = std::make_pair(gv.indexSet().index(vertex),cnt);
  //   std::cout << gfs2node[cnt].first << " :/: " << gfs2node[cnt].second << std::endl;
  //   cnt+=1;
  // }
  // std::sort(gfs2node.begin(), gfs2node.end());


  // V newPoU2(gfs, 0.0);
  // cnt=0;
  // for (const auto& vertex : vertices(gv)){
  // // for(int i=0; i<v_size; i++){
  //   // std::cout << gv.indexSet().index(vertex) << " // " << cnt << std::endl;
  //   // std::cout << i << " // " << std::endl;

  //   // native(newPoU2)[cnt] = tmp[element_nodes_ordering[cnt].second];
  //   native(newPoU2)[cnt] = tmp[element_nodes_ordering[gfs2node[cnt].second].second];
  //   // native(newPoU)[cnt] = tmp[cnt];
  //   // std::cout << cnt << " // " << native(newPoU)[cnt] << std::endl;
  //   cnt+=1;
  // }