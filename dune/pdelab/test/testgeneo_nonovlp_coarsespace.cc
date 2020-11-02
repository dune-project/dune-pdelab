#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/pdelab.hh>

#include <dune/grid/utility/globalindexset.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/common/rangegenerators.hh>
#include <dune/grid/common/gridinfo.hh> // visualize grid information -> could be remove at the end

#include <dune/grid/utility/parmetisgridpartitioner.hh>

#include <dune/pdelab/test/SubDomainGmshWriter.hh>

#if HAVE_PARMETIS
#include <parmetis.h>

#include <algorithm>
#include <vector>

// only enable for ParMETIS because the implementation uses functions that
// are not emulated by scotch
#ifdef PARMETIS_MAJOR_VERSION

/** \brief Returns a vector with a partition number for each element
 *
 * \author Benjamin Bykowski and adapted by Peter Bastian
 *
 * \param gv The grid view to be partitioned
 * \param mpihelper The MPIHelper object, needed to get the MPI communicator. This is needed by the function ParMETIS_V3_PartMeshKway and can unfortunately not be omitted
 * \param parts number of subdomains desired
 *
 * \return std::vector with one uint per All_Partition element.  For each element, the entry is the
 *    number of the partition the element is assigned to. This number is greater or equal zero and smaller as parts.
 *    No partitioning is done, only this vector is computed
 */
template<class GridView>
std::vector<unsigned> parmetis_partitioning (const GridView& gv, const Dune::MPIHelper& mpihelper, int parts) {

#if PARMETIS_MAJOR_VERSION > 3
  typedef idx_t idx_type;
  typedef ::real_t real_type;
#else
  typedef int idx_type;
  typedef float real_type;
#endif // PARMETIS_MAJOR_VERSION > 3

  const unsigned numElements = gv.size(0);

  std::vector<unsigned> part(numElements);

  // Setup parameters for ParMETIS
  idx_type wgtflag = 0;                                  // we don't use weights
  idx_type numflag = 0;                                  // we are using C-style arrays
  idx_type ncon = 1;                                     // number of balance constraints
  idx_type ncommonnodes = 2;                             // number of nodes elements must have in common to be considered adjacent to each other
  idx_type options[4] = {0, 0, 0, 0};                    // use default values for random seed, output and coupling
  idx_type edgecut;                                      // will store number of edges cut by partition
  idx_type nparts = parts;                               // number of partitions to create is a parameter
  std::vector<real_type> tpwgts(ncon*nparts, 1./nparts); // load per subdomain and weight (same load on every process)
  std::vector<real_type> ubvec(ncon, 1.05);              // weight tolerance (same weight tolerance for every weight there is)

  // The difference elmdist[i+1] - elmdist[i] is the number of nodes that are on process i
  std::vector<idx_type> elmdist(nparts+1);
  elmdist[0] = 0;
  std::fill(elmdist.begin()+1, elmdist.end(), gv.size(0)); // all elements are on process zero

  // Create and fill arrays "eptr", where eptr[i] is the number of vertices that belong to the i-th element, and
  // "eind" contains the vertex-numbers of the i-the element in eind[eptr[i]] to eind[eptr[i+1]-1]
  std::vector<idx_type> eptr, eind;
  int numVertices = 0;
  eptr.push_back(numVertices);

  for (const auto& element : elements(gv, Dune::Partitions::interior)) {
    const size_t curNumVertices = Dune::referenceElement<double,GridView::dimension>(element.type()).size(GridView::dimension);

    numVertices += curNumVertices;
    eptr.push_back(numVertices);

    for (size_t k = 0; k < curNumVertices; ++k)
      eind.push_back(gv.indexSet().subIndex(element, k, GridView::dimension));
  }

  // Partition mesh using ParMETIS
  if (0 == mpihelper.rank()) {
    MPI_Comm comm = Dune::MPIHelper::getLocalCommunicator();

#if PARMETIS_MAJOR_VERSION >= 4
    const int OK =
#endif
      ParMETIS_V3_PartMeshKway(elmdist.data(), eptr.data(), eind.data(), NULL, &wgtflag, &numflag,
                               &ncon, &ncommonnodes, &nparts, tpwgts.data(), ubvec.data(),
                               options, &edgecut, reinterpret_cast<idx_type*>(part.data()), &comm);

#if PARMETIS_MAJOR_VERSION >= 4
    if (OK != METIS_OK)
      DUNE_THROW(Dune::Exception, "ParMETIS returned an error code.");
#endif
  }

  return part;
}

#else // PARMETIS_MAJOR_VERSION
#warning "You seem to be using the ParMETIS emulation layer of scotch, which does not work with this file."
#endif

#else // HAVE_PARMETIS
#warning "PARMETIS was not found, please check your configuration"
#endif

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

template <typename V_>
struct AddGatherScatter
{
  static typename V_::value_type gather(const V_ &a, int i)
  {
      return a[i]; // I am sending my value
  }
  static void scatter(V_& a, typename V_::value_type v, int i)
  {
      a[i] += v; // add what I receive to my value
  }
};

// template <typename V_>
// struct GatherGlobalIndex
// {
//   static typename V_::value_type gather(const V_& a, int i)
//   {
//       return a[i]; // I am sending my value
//   }
//   static void scatter(V_& a, typename V_::value_type v, int i)
//   {
//     a[i] = v; // add what I receive to my value
//     // if(a[i]<2)
//     //   a[i] += v; // add what I receive to my value
//   }
// };

/**
 * \brief Takes a partition and extends all subdomains by a number of layers given by overlap
 *
 * \param gv The grid view to be operated on
 * \param partitions number of partitions (or subdomains) the grid has been partitioned into (using the function parmetis_partition)
 * \param partition partition information for each element
 * \param overlap overlap that should be added (zero is fine, then nothing is done)
 * \param mode determines how partitions are grown. Can have the value "vertex" or "element", default is "vertex"
 *
 * If mode has the value "element" the extension is done via element faces. Partitions are extended to neighbors of elements in overlap rounds.
 * If mode has the value "vertex" the extension is done via vertices of the grid. First, partitioning is converted from elements to
 * vertices. Then partitions are extended to neighbors of vertices in overlap rounds. Finally, partitioning
 * is converted back to elements.
 * If mode has neither the value "element" or "vertex" the original partitioning is returned, converted to a set for each element
 *
 * \return std::vector<std::set<unsigned>> for each element stores a set of subomains containing this element
 */
template<class GV>
std::vector<std::set<unsigned>> grow_subdomains (const GV& gv, unsigned partitions, std::vector<unsigned> partition, unsigned overlap, std::string mode="vertex")
{
  const int dim = GV::dimension; // extract dimension (codim of vertices)
  auto& indexset = gv.indexSet(); // to attach data to elements

  if (mode=="vertex") {
    std::vector<std::set<unsigned>> subdomainsv(indexset.size(dim)); // set of subdomains for each vertex

    // initialize subdomain list for each vertex by the partition
    for (const auto& e : elements(gv))
      for (unsigned int i=0; i<e.subEntities(dim); ++i)
        {
          auto v = e.template subEntity<dim>(i);
          subdomainsv[indexset.index(v)].insert(partition[indexset.index(e)]);
        }

    // in each round extend overlap by one
    for (int rounds=0; rounds<overlap; rounds++)
    {
      std::vector<std::set<unsigned>> old(subdomainsv); // copy current state
      for (const auto& e : elements(gv))
        {
          // build union of all partitions in all vertices of the element
          std::set<unsigned> unification;
          for (unsigned int i=0; i<e.subEntities(dim); ++i)
            for (const auto& j : old[indexset.index(e.template subEntity<dim>(i))])
              unification.insert(j);
          // now add union to all vertices (a clique)
          for (const auto& j : unification)
            for (unsigned int i=0; i<e.subEntities(dim); ++i)
              subdomainsv[indexset.index(e.template subEntity<dim>(i))].insert(j);
        }
    }

    // now convert again to elements: element is in subdomain if *all* vertices are in subdomain
    std::vector<std::set<unsigned>> subdomainse(indexset.size(0)); // set of subdomains for each element
    for (const auto& e : elements(gv))
    {
      std::set<unsigned> intersection(subdomainsv[indexset.index(e.template subEntity<dim>(0))]);
      for (unsigned int i=1; i<e.subEntities(dim); ++i)
        {
          std::set<unsigned> update;
          for (const auto& j : subdomainsv[indexset.index(e.template subEntity<dim>(i))])
            if (intersection.count(j)>0) update.insert(j);
                intersection = update;
        }
      subdomainse[indexset.index(e)] = intersection;
    }
    // and we are done
    return subdomainse;
  }

  // now the element mode
  std::vector<std::set<unsigned>> subdomains(indexset.size(0)); // set of subdomains for each element

  // initialize subdomain list for each element by the partition
  for (const auto& e : elements(gv))
    subdomains[indexset.index(e)].insert(partition[indexset.index(e)]);

  if (mode=="element")
    {
      // in each round extend overlap by one
      for (int rounds=0; rounds<overlap; rounds++){
        std::vector<std::set<unsigned>> old(subdomains); // copy current state
        for (const auto& e : elements(gv))
          for (const auto& is : intersections(gv,e))
            if (is.neighbor())
              for (const auto& i : old[indexset.index(is.outside())])
                subdomains[indexset.index(e)].insert(i);
      }
    }
  // and we are done
  return subdomains;
}

void driver(std::string basis_type, std::string part_unity_type, Dune::MPIHelper& helper) {

  // ~~~~~~~~~~~~~~~~~~
  // Create offline hierarchy
  // ~~~~~~~~~~~~~~~~~~
  std::string path_to_storage = "Offline/";

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
  Dune::GmshReader<GRID>::read(factory, "grids/24x24.msh", true, true);
  std::unique_ptr<GRID> grid (factory.createGrid());

  typedef typename GRID::LeafGridView GV;
  auto gv = grid->leafGridView();

  // auto indexset = gv.indexSet();

  const int algebraic_overlap = 5;


  // if (gv.comm().rank()==0){
  //   // Create a vector from GFS ordering to node ordering
  //   int v_size = grid->levelIndexSet(0).size(gv.indexSet().types(0)[1]);
  //   std::vector<std::pair<int, int>> gfs2node(v_size);

  //   int cnt=0;
  //   for (const auto& vertex : vertices(gv)){
  //     assert(cnt<v_size);
  //     gfs2node[cnt] = std::make_pair(gv.indexSet().index(vertex),cnt);
  //     std::cout << gfs2node[cnt].first << " : " << gfs2node[cnt].second << std::endl;
  //     cnt+=1;
  //   }
  // }

  // ~~~~~~~~~~~~~~~~~~
//  Associate to each global element its global ID
  // ~~~~~~~~~~~~~~~~~~
  // Dune::gridlevellist<GRID>(*grid, 0, "");

  // // Getting map global ID to global Index
  // int length;
  // if(gv.comm().rank()==0)
  //   length = grid->levelIndexSet(0).size(grid->levelIndexSet(0).types(0)[0]);
  //   // std::cout << "elements size: " << grid->levelIndexSet(0).size(grid->levelIndexSet(0).types(0)[0]) << std::endl;
  // gv.comm().broadcast(&length, 1, 0);
  // std::vector<long unsigned int> Element_GlobalIndex_to_GlobalID(length);
  // for (const auto& element : elements(levelGridView(*grid, 0)))
  //   Element_GlobalIndex_to_GlobalID[grid->levelIndexSet(0).index(element)] = grid->globalIdSet().template id<0>(element);

  // // if (gv.comm().barrier()==0)
  // for (int i=0; i<Element_GlobalIndex_to_GlobalID.size();i++)
  //   gv.comm().broadcast(&Element_GlobalIndex_to_GlobalID[i], 1, 0);

  // ~~~~~~~~~~~~~~~~~~
//  Associate to each global Vertex its global ID
  // ~~~~~~~~~~~~~~~~~~
  // Dune::gridlevellist<GRID>(*grid, 0, "");

  // Getting map global ID to global Index
  int length;
  if(gv.comm().rank()==0){
    length = grid->levelIndexSet(0).size(gv.indexSet().types(0)[1]);
    std::cout << "nodes size: " << length << std::endl;
  }
  gv.comm().broadcast(&length, 1, 0);
  std::vector<long unsigned int> Vertices_GlobalIndex_to_GlobalID(length);

  if(gv.comm().rank()==0)
    for (const auto& vertex : vertices(gv))
      Vertices_GlobalIndex_to_GlobalID[gv.indexSet().index(vertex)] = gv.grid().globalIdSet().id(vertex);


  // // if (gv.comm().barrier()==0)
  for (int i=0; i<Vertices_GlobalIndex_to_GlobalID.size();i++)
    gv.comm().broadcast(&Vertices_GlobalIndex_to_GlobalID[i], 1, 0);
    // std::cout << i << " : " << Vertices_GlobalIndex_to_GlobalID[i] << std::endl;


  // Transfer partitioning from ParMETIS to our grid
#if PARMETIS_MAJOR_VERSION
  // ~~~~~~~~~~~~~~~~~~
  // Save domain decomposition (meshes)
  // ~~~~~~~~~~~~~~~~~~
  // First recreate the overlapped subdomain
  unsigned subdomains = grid->comm().size();
  auto partition = parmetis_partitioning(gv,helper,subdomains);
  unsigned uoverlap = algebraic_overlap + 0u;
  std::string extensionmethod = "vertex";
  auto overlappingsubdomains = grow_subdomains(gv,subdomains,partition,uoverlap,extensionmethod);

  // Write subdomain as a msh file
  Dune::SubDomainGmshWriter<GV> sdwriter(gv, overlappingsubdomains, subdomains);
  sdwriter.write(path_to_storage, "subdomain.msh");

  grid->loadBalance(partition, 0);
#else
  grid->loadBalance();
#endif

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

  using Dune::PDELab::Backend::native;

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
  go.residual(x,d);


  // ~~~~~~~~~~~~~~~~~~
//  Solving process begin here: First some parameters
  // ~~~~~~~~~~~~~~~~~~
  double eigenvalue_threshold = -1;
  int nev = 3;
  int nev_arpack = nev;
  // double shift = 0.001;

  // ~~~~~~~~~~~~~~~~~~
//  ADAPTER
  // ~~~~~~~~~~~~~~~~~~
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

  auto rank = adapter.gridView().comm().rank();

  V vect(gfs, 0.0);
  adapter.restrictVector(*part_unity, native(vect));

  Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,Dune::refinementLevels(0));
  Dune::PDELab::vtk::DefaultFunctionNameGenerator fieldname;
  fieldname.prefix("PoU");
  Dune::PDELab::addSolutionToVTKWriter(vtkwriter,gfs,vect,fieldname);
  vtkwriter.write("PoU",Dune::VTK::ascii);

  // ~~~~~~~~~~~~~~~~~~
  // Save new2old / old2new index to be able to reproduce extension/reduction online
  // ~~~~~~~~~~~~~~~~~~

  auto new2old_localindex = adapter.get_new2old_localindex();
  auto old2new_localindex = adapter.get_old2new_localindex();

  Dune::BlockVector<Dune::FieldVector<int, 1>> n2o(new2old_localindex.size()), o2n(old2new_localindex.size());
  for (int i=0; i<new2old_localindex.size();i++) { n2o[i]=new2old_localindex[i]; }
  for (int i=0; i<old2new_localindex.size();i++) { o2n[i]=old2new_localindex[i]; }

  std::string filename_n2o = path_to_storage + std::to_string(rank) + "_restrict.mm";
  Dune::storeMatrixMarket(n2o, filename_n2o, 15);
  std::string filename_o2n = path_to_storage + std::to_string(rank) + "_extend.mm";
  Dune::storeMatrixMarket(o2n, filename_o2n, 15);



  // ~~~~~~~~~~~~~~~~~~
// Get local index to global id
  // ~~~~~~~~~~~~~~~~~~

  // auto vertices_globalID_by_subdomain = adapter.get_globalid();

  using GlobalId = typename GV::Grid::GlobalIdSet::IdType;
  using AttributedLocalIndex = Dune::ParallelLocalIndex<Attribute>;
  using ParallelIndexSet = Dune::ParallelIndexSet<GlobalId,AttributedLocalIndex,256>;
  auto lookup = Dune::GlobalLookupIndexSet<ParallelIndexSet>(*adapter.getEpis().parallelIndexSet(), adapter.getExtendedSize());

  // int cntfound=0;
  // int cntnot=0;
  Dune::BlockVector<Dune::FieldVector<int, 1>> global_indices(adapter.getExtendedSize());
  for (int incrSD=0; incrSD<global_indices.size(); incrSD++){
    auto it = std::find(Vertices_GlobalIndex_to_GlobalID.begin(), Vertices_GlobalIndex_to_GlobalID.end(), lookup.pair(incrSD)->global());
    if (it != Vertices_GlobalIndex_to_GlobalID.end()){
      // if (rank==0){
      //   cntfound+=1;
      //   std::cout << "found : " << *it << std::endl;
      // }
      global_indices[incrSD] = std::distance(Vertices_GlobalIndex_to_GlobalID.begin(), it);
    }
    // else {
    //   if (rank==0){
    //     cntnot+=1;
    //     std::cout << "not found : " << lookup.pair(incrSD)->global() << std::endl;
    //   }
    // }
  }

  std::string filename_GI = path_to_storage + std::to_string(rank) + "_GI.mm";
  Dune::storeMatrixMarket(global_indices, filename_GI, 15);


  // if(rank==0){
  //   for (int i=0; i<vertices_globalID_by_subdomain.size();i++){
  //     auto my_pair = lookup[vertices_globalID_by_subdomain[i]];
  //     std::cout << vertices_globalID_by_subdomain[i] << std::endl;
  //     std::cout << my_pair << std::endl;
  //     std::cout << global_indices[i] << std::endl;
  //   }
  // }
  // auto vertices_globalID_by_subdomain = adapter.get_globalid();

  // if (rank==1){
  //   std::cout << "local size : " << vertices_globalID_by_subdomain.size() << std::endl;
  //   std::cout << "local size with ovlp : " << adapter.getExtendedSize() << std::endl;
  // }

  // if (rank==1)
  //   for (int i=0; i<vertices_globalID_by_subdomain.size(); i++)
  //     std::cout << i << " : " << vertices_globalID_by_subdomain[i] << std::endl;

  // if (rank==1){
  //   std::cout << std::endl;
  //   std::cout << std::endl;
  // }

  // std::vector<long unsigned int> vertices_globalID_by_ovlp_subdomain(adapter.getExtendedSize());
  // for (int i=0; i<adapter.getExtendedSize(); i++)
  //   vertices_globalID_by_ovlp_subdomain[i] = 1;

  // std::cout << "new2old_localindex.size() : " << old2new_localindex.size() << std::endl;
  // for (int i=0; i<old2new_localindex.size(); i++)
  //   vertices_globalID_by_ovlp_subdomain[old2new_localindex[i]] = vertices_globalID_by_subdomain[i];

  // if (rank==1)
  //   for (int i=0; i<vertices_globalID_by_ovlp_subdomain.size(); i++)
  //     std::cout << i << " : " << vertices_globalID_by_ovlp_subdomain[i] << std::endl;

  // communicator->forward<GatherGlobalIndex<std::vector<long unsigned int>>>(vertices_globalID_by_ovlp_subdomain,vertices_globalID_by_ovlp_subdomain); // make function known in other subdomains

  // if (rank==1){
  //   std::cout << std::endl;
  //   std::cout << std::endl;
  // }

  // if (rank==1)
  //   for (int i=0; i<vertices_globalID_by_ovlp_subdomain.size(); i++)
  //     std::cout << i << " : " << vertices_globalID_by_ovlp_subdomain[i] << std::endl;

  // ~~~~~~~~~~~~~~~~~~
  // Associate local ovlp vectors to global indices
  // ~~~~~~~~~~~~~~~~~~

  // int cntnot=0;
  // int cntnot=0;
  // Dune::BlockVector<Dune::FieldVector<int, 1>> global_indices(adapter.getExtendedSize());
  // for (int incrSD=0; incrSD<global_indices.size(); incrSD++){
  //   auto it = std::find(Vertices_GlobalIndex_to_GlobalID.begin(), Vertices_GlobalIndex_to_GlobalID.end(), vertices_globalID_by_ovlp_subdomain[incrSD]);

  //   if (it != Vertices_GlobalIndex_to_GlobalID.end()){
  //     if (rank==0){
  //       cntfound+=1;
  //       std::cout << "found : " << *it << std::endl;
  //     }
  //     global_indices[incrSD] = std::distance(Vertices_GlobalIndex_to_GlobalID.begin(), it);
  //   } else {
  //     if (rank==0){
  //       cntnot+=1;
  //       std::cout << "not found : " << vertices_globalID_by_ovlp_subdomain[incrSD]<< std::endl;
  //     }
  //   }
  // }

  // if (rank==0){
  //   std::cout << "cntfound : " << cntfound << ", cntnot : " << cntnot << std::endl;
  // }
  // std::string filename_GI = path_to_storage + std::to_string(rank) + "_GI.mm";
  // Dune::storeMatrixMarket(global_indices, filename_GI, 15);

  // ~~~~~~~~~~~~~~~~~~
// Get local index to global id 1
  // ~~~~~~~~~~~~~~~~~~

  // auto vertices_globalID_by_subdomain = adapter.get_globalid();
  // // std::cout << "local size : " << vertices_globalID_by_subdomain.size() << std::endl;

  // // ~~~~~~~~~~~~~~~~~~
  // // Associate local vectors to global indices
  // // ~~~~~~~~~~~~~~~~~~

  // int cntfound=0;
  // int cntnot=0;
  // Dune::BlockVector<Dune::FieldVector<int, 1>> global_indices(vertices_globalID_by_subdomain.size());
  // for (int incrSD=0; incrSD<vertices_globalID_by_subdomain.size(); incrSD++){
  //   auto it = std::find(vertices_globalID_by_subdomain.begin(), vertices_globalID_by_subdomain.end(), vertices_globalID_by_subdomain[incrSD]);

  //   if (it != vertices_globalID_by_subdomain.end()){
  //     global_indices[incrSD] = std::distance(vertices_globalID_by_subdomain.begin(), it);
  //     if(rank==1)
  //       std::cout << global_indices[incrSD] << std::endl;
  //   }
  // }

  // Dune::BlockVector<Dune::FieldVector<int, 1>> ovlp_global_indices(adapter.getExtendedSize());
  // for (int i=0; i<adapter.getExtendedSize(); i++)
  //   ovlp_global_indices[i] = -1;

  // std::cout << "new2old_localindex.size() : " << new2old_localindex.size() << std::endl;
  // for (int i=0; i<new2old_localindex.size(); i++)
  //   ovlp_global_indices[i] = global_indices[new2old_localindex[i]];

  // if (rank==2)
  //   for (int i=0; i<ovlp_global_indices.size(); i++)
  //     std::cout << i << " : " << ovlp_global_indices[i] << std::endl;


  // communicator->forward<GatherGlobalIndex<Dune::BlockVector<Dune::FieldVector<int, 1>>>>(ovlp_global_indices,ovlp_global_indices); // make function known in other subdomains

  // if (rank==1)
  //   for (int i=0; i<ovlp_global_indices.size(); i++)
  //     std::cout << i << " : " << ovlp_global_indices[i] << std::endl;

  // ~~~~~~~~~~~~~~~~~~
// Save PoU
  // ~~~~~~~~~~~~~~~~~~
  Dune::BlockVector<Dune::FieldVector<K, 1>> PoU((*part_unity).size());
  for (int i=0; i<(*part_unity).size();i++)
    PoU[i]=(*part_unity)[i];
  std::string filename_PoU = path_to_storage + std::to_string(rank) + "_PoU.mm";
  Dune::storeMatrixMarket(PoU, filename_PoU, 15);


  // ~~~~~~~~~~~~~~~~~~
  // Save neighbor ranks vector
  // ~~~~~~~~~~~~~~~~~~

  std::vector<int> neighbor_ranks = adapter.findNeighboringRanks();

  Dune::BlockVector<Dune::FieldVector<int, 1>> NR(neighbor_ranks.size());
  for (int i=0; i<neighbor_ranks.size();i++) { NR[i]=neighbor_ranks[i]; }

  std::string filename_NR = path_to_storage + std::to_string(rank) + "_neighborRanks.mm";
  Dune::storeMatrixMarket(NR, filename_NR, 15);

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  Subdomain basis computations
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  std::shared_ptr<Dune::PDELab::SubdomainBasis<Vector>> subdomainbasis;
  subdomainbasis = std::make_shared<Dune::PDELab::NonoverlappingGenEOBasis<GO, Matrix, Vector>>(adapter, A_extended, A_extended, part_unity, eigenvalue_threshold, nev, nev_arpack);

  // ~~~~~~~~~~~~~~~~~~
  //  Particular solution
  // ~~~~~~~~~~~~~~~~~~
  // TODO: compute the particular solution just on certain subdomains not all

  // Extend the load vector d (right hand side) of the fine to space to its virtual overlap
  Vector b(adapter.getExtendedSize());

  //~ Add a solution to another rhs here fill with ones or
  // for(auto it = b.begin(); it!=b.end(); ++it){
  //   b[it.index()] += 1.0;
    // std::srand(std::time(0));
    // b[it.index()] = -1.0 + 2.0* (std::rand()+0.0) / (RAND_MAX + 1.0);
    // if(adapter.gridView().comm().rank()==0) std::cout<< b[it.index()] << std::endl;
  // }

  // Use of the real rhs
  adapter.extendVector(native(d), b);

  communicator->forward<AddGatherScatter<Vector>>(b,b); // make function known in other subdomains

  // ~~~~~~~~~~~~~~~~~~
  // Set up boundary condition for the computation of the particular solution
  // ~~~~~~~~~~~~~~~~~~
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

  // ~~~~~~~~~~~~~~~~~~
  // Solve the local fine space with ovlp system
  // ~~~~~~~~~~~~~~~~~~
  Vector ui(adapter.getExtendedSize());
  Dune::UMFPack<Matrix> subdomain_solver(A_dirichlet, false);
  Dune::InverseOperatorResult result1;
  subdomain_solver.apply(ui,b_cpy,result1);
  // subdomainbasis->append(ui);

  // ~~~~~~~~~~~~~~~~~~
  // Save local basis sizes
  // ~~~~~~~~~~~~~~~~~~
  // Communicate local coarse space dimensions
  Dune::BlockVector<Dune::FieldVector<int, 1>> local_basis_sizes(adapter.gridView().comm().size());
  int local_size = subdomainbasis->basis_size();
  adapter.gridView().comm().allgather(&local_size, 1, local_basis_sizes.data());
  if (rank==0) {
    std::string filename_ls = path_to_storage + "localBasisSizes.mm";
    Dune::storeMatrixMarket(local_basis_sizes, filename_ls, 15);
  }

  // ~~~~~~~~~~~~~~~~~~
  // Save local basis
  // ~~~~~~~~~~~~~~~~~~
  for (int basis_index = 0; basis_index < subdomainbasis->basis_size(); basis_index++) {
    std::string filename_EV = path_to_storage + std::to_string(rank) + "_EV_" + std::to_string(basis_index) + ".mm";
    Dune::storeMatrixMarket(*subdomainbasis->get_basis_vector(basis_index), filename_EV, 15);
  }

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  Coarse space construction and saving
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  auto coarse_space = std::make_shared<Dune::PDELab::NonoverlappingSubdomainProjectedCoarseSpace<GV, Matrix, Vector>>(adapter, gv, *A_extended, subdomainbasis, verbose);
  // Get the coarse matrix from the coarse space
  std::shared_ptr<CoarseMatrix> AH = coarse_space->get_coarse_system();
  // Save the coarse space
  std::string filename = path_to_storage + "OfflineAH.mm";
  Dune::storeMatrixMarket(*AH, filename, 15);

  // // Save the sizes of local basis
  // coarse_space->local_basis_sizes_to_file(path_to_storage + "local_basis_sizes");

  // ~~~~~~~~~~~~~~~~~~
//  Coarse space right hand side and saving
  // ~~~~~~~~~~~~~~~~~~
  // Use the correct rhs to solve the final system
  // adapter.extendVector(native(d), b); // if comment: this could already been declared and filled for the particular solution
  // Initializate a load vector (right hand side) in the coarse space : coarse_d = RH * b
  CoarseVector coarse_d(coarse_space->basis_size(), coarse_space->basis_size());
  coarse_space->restrict(b,coarse_d);
  // Save the coarse rhs
  if (rank==0){
    std::string filename_cb = path_to_storage + "OfflineCoarseb.mm";
    Dune::storeMatrixMarket(coarse_d, filename_cb, 15);
  }

  // ~~~~~~~~~~~~~~~~~~
// Solve the coarse space system  ::TODO:: MAKE THAT OPTIONAL for a normal offline go
  // ~~~~~~~~~~~~~~~~~~
  // Objective :  find x = RH^T * AH^-1 * RH * b
  // Use of UMFPack to solve the problem [AH * coarse_v = RH * b] instead of inversing AH :: objective is to have [coarse_v = AH^-1 *  RH * b]
  Dune::UMFPack<CoarseMatrix> coarse_solver(*AH, false); // Is there something better than UMFPACK?
  CoarseVector coarse_v(coarse_space->basis_size(), coarse_space->basis_size());
  Dune::InverseOperatorResult result;
  coarse_solver.apply(coarse_v,coarse_d,result);

  // ~~~~~~~~~~~~~~~~~~
// Prolongate the solution in order to have v_fine_vovlp on the fine space : v_fine_vovlp = RH^T * coarse_v = RH^T * AH^-1 * RH * b
  // ~~~~~~~~~~~~~~~~~~
  Vector v_fine_vovlp(adapter.getExtendedSize());
  coarse_space->prolongate(coarse_v, v_fine_vovlp);
  communicator->forward<AddGatherScatter<Vector>>(v_fine_vovlp, v_fine_vovlp);
  V v(gfs, 0.0);
  adapter.restrictVector(v_fine_vovlp, v);

  native(x) -= v;


  // Communicate local fine space dimensions to be able to initialize solution vector in the online phase
  Dune::BlockVector<Dune::FieldVector<int, 1>> local_nonovlp_sizes(adapter.gridView().comm().size());
  int local_nonovlp_size = native(v).size();
  adapter.gridView().comm().allgather(&local_nonovlp_size, 1, local_nonovlp_sizes.data());
  if (rank==0) {
    std::string filename_lns = path_to_storage + "localNovlpSizes.mm";
    Dune::storeMatrixMarket(local_basis_sizes, filename_lns, 15);
  }

  // Write solution to VTK
  // Dune::VTKWriter<GV> vtkwriter(gv);
  // typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
  // DGF xdgf(gfs,x);
  // typedef Dune::PDELab::VTKGridFunctionAdapter<DGF> ADAPT;
  // auto adapt = std::make_shared<ADAPT>(xdgf,"solution");
  // vtkwriter.addVertexData(adapt);
  // vtkwriter.write("solution");

  // ~~~~~~~~~~~~~~~~~~
// Visualise all the basis in vtk format
  // ~~~~~~~~~~~~~~~~~~
  // for (int basis_index = 0; basis_index < subdomainbasis->basis_size(); basis_index++) {
  //     Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,Dune::refinementLevels(0));
  //     V vect(gfs, 0.0);
  //     adapter.restrictVector(native(*subdomainbasis->get_basis_vector(basis_index)), native(vect));

  //     int rank = adapter.gridView().comm().rank();
  //     std::string filename = "BasisVector_"+std::to_string(basis_index);

  //     Dune::PDELab::vtk::DefaultFunctionNameGenerator fieldname;
  //     fieldname.prefix("EV");

  //     Dune::PDELab::addSolutionToVTKWriter(vtkwriter,gfs,vect,fieldname);
  //     vtkwriter.write(filename,Dune::VTK::ascii);
  // }
}


int main(int argc, char **argv)
{
  using Dune::PDELab::Backend::native;

  // initialize MPI, finalize is done automatically on exit
  Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc,argv);

  driver("geneo", "standard", helper);

  return 0;
}
