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

// #if HAVE_UG
// #include <dune/grid/uggrid/uggridfactory.hh>

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
    static void scatter(V_ &a, typename V_::value_type v, int i)
    {
        a[i] += v; // add what I receive to my value
    }
};

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
          for (const auto& j : subdomainsv[indexset.index(e.template subEntity<dim>(i))]) {
            if (intersection.count(j)>0)
              update.insert(j);
            if (j>2)
              std::cout << "intersection: " << j << std::endl;
          }
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

  if (mode=="element") {
    // in each round extend overlap by one
    for (int rounds=0; rounds<overlap; rounds++) {
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

  // std::vector<unsigned> partition(grid->leafGridView().size(0));
  // std::cout << grid->leafGridView().size(0) << std::endl;
  // for (int i=0; i< partition.size(); i++) {
  //   if (i < partition.size()/2 ) {partition[i]=0;}
  //   else {partition[i]=1;}
  // }

  unsigned subdomains = 3; // grid->comm().size();
  unsigned uoverlap = 1u;
  std::string extensionmethod = "vertex";

  auto partition = parmetis_partitioning(gv,helper,subdomains);

  // ~~~~~~~~~~~~~~~~~~
  // Save domain decomposition
  // ~~~~~~~~~~~~~~~~~~
  // First recreate the overlapped subdomain
  auto overlappingsubdomains = grow_subdomains(gv,subdomains,partition,uoverlap,extensionmethod);

  // test overlap generation
  auto& indexset = gv.indexSet(); // to attach data to elements
  std::cout << indexset.size(0) << " elements" << std::endl;
  std::cout << indexset.size(dim) << " vertices" << std::endl;
  std::cout << subdomains << " subdomains" << std::endl;
  int vp = 0;
  std::vector<unsigned> adomain(indexset.size(0));
  std::vector<unsigned> k0(indexset.size(0));
  size_t k_0 = 0;
  unsigned count_elements=0;
  for (const auto& e : elements(gv))
  {
    adomain[indexset.index(e)] = 0;
    if (overlappingsubdomains[indexset.index(e)].count(vp)>0) adomain[indexset.index(e)] = 1;
    if (partition[indexset.index(e)]==vp) adomain[indexset.index(e)] = 2;
    k0[indexset.index(e)] = overlappingsubdomains[indexset.index(e)].size();
    k_0 = std::max(k_0,overlappingsubdomains[indexset.index(e)].size());
    count_elements += overlappingsubdomains[indexset.index(e)].size();
  }
  std::cout << k_0 << " k_0" << std::endl;
  std::cout << indexset.size(0) << " elements, " << count_elements << " in all subdomains, ratio=" << ((double)count_elements)/indexset.size(0) << std::endl;

  // Write subdomain as a msh file
  Dune::SubDomainGmshWriter<GV> sdwriter(gv, overlappingsubdomains, subdomains);
  sdwriter.write("", "subdomain.msh");


  grid->loadBalance(partition, 0);

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


  // Write solution to VTK
  Dune::VTKWriter<GV> vtkwriter(gv);
  typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
  DGF xdgf(gfs,x);
  typedef Dune::PDELab::VTKGridFunctionAdapter<DGF> ADAPT;
  auto adapt = std::make_shared<ADAPT>(xdgf,"solution");
  vtkwriter.addVertexData(adapt);
  vtkwriter.write("vtksol");

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

// #endif


int main(int argc, char **argv)
{
  using Dune::PDELab::Backend::native;

  // initialize MPI, finalize is done automatically on exit
  Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc,argv);

  driver("geneo", "standard", helper);

  return 0;
}
