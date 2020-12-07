#ifndef DUNE_PDELAB_OFFLINE_TOOLS_HH
#define DUNE_PDELAB_OFFLINE_TOOLS_HH

template <typename V_>
struct AddGatherScatter{
  static typename V_::value_type gather(const V_ &a, int i)
  {
      return a[i]; // I am sending my value
  }
  static void scatter(V_& a, typename V_::value_type v, int i)
  {
      a[i] += v; // add what I receive to my value
  }
};

template<typename GRID>
std::vector<long unsigned int> extractNodalGID(std::unique_ptr<GRID>& grid, int verbose=0) {
  auto gv = grid->leafGridView();

  int length;
  if(gv.comm().rank()==0)
    length = grid->levelIndexSet(0).size(GRID::LeafGridView::dimension); // 0-> elements ; GRID::LeafGridView::dimension-> nodes

  gv.comm().broadcast(&length, 1, 0);

  if(verbose>0)
    std::cout << "Rank:" << gv.comm().rank() << " -> NodalGID's size: " << length << std::endl;

  std::vector<long unsigned int> NodalGID(length);

  if(gv.comm().rank()==0)
    for (const auto& vertex : vertices(gv))
      NodalGID[gv.indexSet().index(vertex)] = gv.grid().globalIdSet().id(vertex);

  for (int i=0; i<NodalGID.size();i++)
    gv.comm().broadcast(&NodalGID[i], 1, 0);

  return NodalGID;
}

template<typename GRID>
std::vector<long unsigned int> extractElementGID(std::unique_ptr<GRID>& grid, int verbose=0) { // Not used but maybe later for non-linear element
  auto gv = grid->leafGridView();
  // Getting map global ID to global Index
  int length;
  if(gv.comm().rank()==0)
    length = grid->levelIndexSet(0).size(0); // 0-> elements ; GRID::LeafGridView::dimension-> nodes

  gv.comm().broadcast(&length, 1, 0);

  if(verbose>0)
    std::cout << "Rank:" << gv.comm().rank() << " -> NodalGID's size: " << length << std::endl;

  std::vector<long unsigned int> ElementGID(length);
  for (const auto& element : elements(levelGridView(*grid, 0)))
    ElementGID[grid->levelIndexSet(0).index(element)] = grid->globalIdSet().template id<0>(element);

  for (int i=0; i<ElementGID.size();i++)
    gv.comm().broadcast(&ElementGID[i], 1, 0);

  return ElementGID;
}

template<typename vector1i, class GridView, class Matrix, class Vector>
vector1i extractGlobalIndices(Dune::NonoverlappingOverlapAdapter<GridView, Vector, Matrix>& adapter, std::vector<long unsigned int>& NodalGID, int verbose=0) {

  using Attribute = Dune::EPISAttribute;
  using GlobalId = typename GridView::Grid::GlobalIdSet::IdType; // long unsigned int
  using AttributedLocalIndex = Dune::ParallelLocalIndex<Attribute>;
  using ParallelIndexSet = Dune::ParallelIndexSet<GlobalId,AttributedLocalIndex,256>;
  auto lookup = Dune::GlobalLookupIndexSet<ParallelIndexSet>(*adapter.getEpis().parallelIndexSet(), adapter.getExtendedSize());

  int cntnot=0;
  vector1i global_indices(adapter.getExtendedSize());
  for (int incrSD=0; incrSD<global_indices.size(); incrSD++){
    auto it = std::find(NodalGID.begin(), NodalGID.end(), lookup.pair(incrSD)->global());
    if (it != NodalGID.end()){
      global_indices[incrSD] = std::distance(NodalGID.begin(), it);
    }
    else {
      cntnot+=1;
      if(verbose>0)
        std::cout << "Rank:" << adapter.gridView().comm().rank() << " not found : " << lookup.pair(incrSD)->global() << std::endl;
    }
  }
  assert(cntnot==0);
  return global_indices;
}

template <typename GridView, typename Vector, typename Matrix>
class ParticularSolution {

  public:

  ParticularSolution(Dune::NonoverlappingOverlapAdapter<GridView, Vector, Matrix>& adapter, std::shared_ptr<Matrix> A, const Vector& part_unity)
  : adapter_(adapter),
    A_(A)
  {

    b_.resize(adapter_.getExtendedSize());

  }

  void constantRHS(double value=1.0) {
    for(auto it = b_.begin(); it!=b_.end(); ++it)
      b_[it.index()] += value;

  }

  void randomRHS(double value=1.0) {
    for(auto it = b_.begin(); it!=b_.end(); ++it){
      std::srand(std::time(0));
      b_[it.index()] = value*(-1.0 + 2.0* (std::rand()+0.0) / (RAND_MAX + 1.0));
    }

  }

  void exactRHS(Vector d) {
    using Attribute = Dune::EPISAttribute;
    Dune::AllSet<Attribute> allAttribute;
    auto allinterface = std::shared_ptr<Dune::Interface>(new Dune::Interface());
    allinterface->build(*adapter_.getRemoteIndices(),allAttribute,allAttribute); // all to all communication
    // build up buffered communicator allowing communication over a dof vector
    auto communicator = std::shared_ptr<Dune::BufferedCommunicator>(new Dune::BufferedCommunicator());
    communicator->build<Vector>(*allinterface);

    adapter_.extendVector(d, b_);
    communicator->forward<AddGatherScatter<Vector>>(b_,b_); // make function known in other subdomains

  }

  void solveAndAppend(Dune::PDELab::SubdomainBasis<Vector>& subdomainbasis) {
    Vector ui(adapter_.getExtendedSize());
    Dune::UMFPack<Matrix> subdomain_solver(*A_, false);
    Dune::InverseOperatorResult result1;
    subdomain_solver.apply(ui,b_,result1);
    subdomainbasis.append(ui);
  }

  private:
  Dune::NonoverlappingOverlapAdapter<GridView, Vector, Matrix> adapter_;
  std::shared_ptr<Matrix> A_ = nullptr;
  Vector b_;
};



#endif // DUNE_PDELAB_OFFLINE_TOOLS_HH
