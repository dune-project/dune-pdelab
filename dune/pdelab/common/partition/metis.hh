#ifndef DUNE_PDELAB_COMMON_PARTITION_METIS_HH
#define DUNE_PDELAB_COMMON_PARTITION_METIS_HH


#include <dune/pdelab/common/partition/halo.hh>
#include <dune/pdelab/common/partition/coloring.hh>
// #include <dune/pdelab/common/trace.hh>

#include <dune/grid/concepts/entity.hh>
#include <dune/grid/concepts/gridview.hh>
#include <dune/grid/common/partitionset.hh>

#include <dune/common/iteratorfacades.hh>

#if HAVE_METIS
#include <metis.h>
#endif

#include <vector>
#include <memory>

namespace Dune::PDELab::inline Experimental::EntitySetPartition {

namespace Impl {

/**
 * @brief Container of entity seeds
 * @details This container stores entity seeds and converts them into rvalue entities on traversal.
 * Models a range of entities
 * @todo Improve documentation
 *
 * @tparam GridView  grid view to transform entity seeds into entities
 * @tparam Entity    Entity to "store" in the container
 */
template<Dune::Concept::GridView GridView, Dune::Concept::Entity Entity>
class EntitySeedContainer {
  using EntitySeed = typename Entity::EntitySeed;

  class ConstIterator : public RandomAccessIteratorFacade<ConstIterator, Entity, Entity> {
    using typename RandomAccessIteratorFacade<ConstIterator, Entity, Entity>::DifferenceType;
  public:

    explicit ConstIterator(typename std::vector<EntitySeed>::const_iterator seed_it, const GridView& grid_view)
      : _seed_it{std::move(seed_it)}
      , _grid_view{&grid_view}
    {}

    explicit ConstIterator() {}

    bool equals(const ConstIterator &other) const
    {
      return _seed_it == other._seed_it;
    }

    Entity dereference() const {
      return _grid_view->grid().entity(*_seed_it);
    }

    void increment() {
      ++_seed_it;
    }

    void decrement(){
      --_seed_it;
    }

    void advance(DifferenceType n) {
      std::advance(_seed_it, n);
    }

    DifferenceType distanceTo(const ConstIterator& other) const {
      return std::distance(_seed_it, other._seed_it);
    }

  private:
    typename std::vector<EntitySeed>::const_iterator _seed_it;
    GridView const * _grid_view;
  };

public:

  explicit EntitySeedContainer(const GridView& grid_view)
    : _grid_view{grid_view}
  {}

  [[nodiscard]] ConstIterator begin() const {
    return ConstIterator{_seeds.begin(), _grid_view};
  }

  [[nodiscard]] ConstIterator end() const {
    return ConstIterator{_seeds.end(), _grid_view};
  }

  std::vector<EntitySeed>& seeds() {
    return _seeds;
  }

private:
  std::vector<EntitySeed> _seeds;
  GridView _grid_view;
};


/**
 * @brief Partition set with entities split on several patches using METIS
 * @details This is a partition on an grid view where entities seeds of the grid
 * view are split using METIS and assigned to different patches of one label set.
 *
 * @tparam GV     Grid view to partition
 * @tparam codim  Codimension to partition
 * @tparam pit    Parallel partition to operate on
 */
template<Dune::Concept::GridView GV, int codim = 0, PartitionIteratorType pit = PartitionIteratorType::All_Partition>
class MetisMixin {
public:
  //! Uderlying grid view
  using GridView = GV;
  //! Entity being partitioned
  using Entity = typename GridView::Codim<codim>::Entity;
  //! Range of entities grouped by a patch
  using PatchSet = EntitySeedContainer<GridView, Entity>;
  //! Range of patches grouped by a label
  using LabelSet = std::vector<PatchSet>;
  //! Range of labels
  using PartitionSet = std::array<LabelSet,1>;

  /**
   * @brief Construct a Simple partition
   *
   * @param grid_view  Grid view to operate on
   * @param patches Number of patches to split the entity set
   */
  explicit MetisMixin(const GridView& grid_view, std::size_t patches)
   : _partition_set{std::make_shared<PartitionSet>()}
   , _grid_view{grid_view}
  {
    update(patches);
  }

  //! Uderlying grid view
  [[nodiscard]] GridView gridView() const noexcept { return _grid_view; }

  //! Range of the partition set
  [[nodiscard]] const PartitionSet& range() const noexcept { return *_partition_set; }

protected:
  /**
   * @brief Sets a new set of patches
   * @details The number of patches will be set up anew
   *
   * @param patches  Number of patches in the set
   */
  void update(std::size_t patches) {
    // TRACE_EVENT("dune", "EntitySet::updatePartition");
    if (patches == 0)
      DUNE_THROW(InvalidStateException, "There cannot be zero partitions");
#if HAVE_METIS
    const std::size_t dimension = GridView::dimension;
    // setup METIS parameters
    idx_t parts = patches;                              // number of partitions
    idx_t ncommonnodes = 2;                             // number of nodes elements must have in common to be considered adjacent to each other
    idx_t edgecut = std::numeric_limits<idx_t>::max();  // will store number of edges cut by partition
    idx_t options[METIS_NOPTIONS];                      // use default values for random seed, output and coupling
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_NUMBERING] = 0;

    std::vector<idx_t> cells, nodes, element_part, node_part;
    std::vector<typename Entity::EntitySeed> seeds;
    cells.reserve(_grid_view.size(0));
    nodes.reserve(_grid_view.size(dimension));

    // create graph of elements and vertices
    int vertices = 0;
    cells.push_back(vertices);
    for (const auto& element : entities(_grid_view, Dune::Codim<codim>{}, Dune::partitionSet<pit>())) {
      const auto& ref_element = referenceElement<double, dimension>(element.type());
      edgecut = std::min<idx_t>(edgecut, ref_element.size(1));
      vertices += ref_element.size(dimension);
      cells.push_back(vertices);
      seeds.emplace_back(element.seed());

      for (int k = 0; k != ref_element.size(dimension); ++k)
        nodes.push_back(_grid_view.indexSet().subIndex(element, k, dimension));
    }

    idx_t element_count = cells.size()-1;
    idx_t node_count = nodes.size();
    element_part.assign(element_count, 0);
    node_part.assign(node_count, 0);

    // actual partition of elements
    auto result = METIS_PartMeshDual(
      &element_count,
      &node_count,
      cells.data(),
      nodes.data(),
      nullptr,
      nullptr,
      &ncommonnodes,
      &parts,
      nullptr,
      options,
      &edgecut,
      element_part.data(),
      node_part.data()
    );

    if (result != METIS_OK)
      DUNE_THROW(InvalidStateException, "Metis could not partition the grid");

    LabelSet label_set(patches, PatchSet{_grid_view});

    // set entity seeds to each assigned patch
    for (std::size_t i = 0; i != element_part.size(); ++i)
      label_set[element_part[i]].seeds().emplace_back(std::move(seeds[i]));

    (*_partition_set)[0] = std::move(label_set);
#else
    DUNE_THROW(InvalidStateException, "You are trying to use METIS partitioner, but METIS was not found!");
#endif
  }

private:
  std::shared_ptr<PartitionSet> _partition_set;
  GridView _grid_view;
};

} //namespace Impl


/**
 * @brief Partition set with entities split on several patches using METIS
 * @details This is a partition on an grid view where entities seeds of the grid
 * view are split using METIS and assigned to different patches of one label set.
 * The patches may be neigboring each other, thus, inducing a shared memory region on
 * several entities of the set.
 * This class models the EntitySetPartition concept.
 *
 * @tparam GV     Grid view to partition
 * @tparam codim  Codimension to partition
 * @tparam pit    Parallel partition to operate on
 */
template<Dune::Concept::GridView GV, int codim = 0, PartitionIteratorType pit = PartitionIteratorType::All_Partition>
struct Metis
  : public Impl::MetisMixin<GV, codim, pit>
  , public Impl::UncoloredOverlapMixin<GV>
{

  /**
   * @brief Construct a Metis uncolored entiy set partition
   *
   * @param grid_view   The grid view to operate on
   * @param patches     Number of patches to have in total
   * @param halo        Distance another entity in the same label set is considered connected
   */
  explicit Metis(const GV& grid_view, std::size_t patches, std::size_t halo)
   : Impl::MetisMixin<GV, codim, pit>{grid_view, patches}
   , Impl::UncoloredOverlapMixin<GV>::UncoloredOverlapMixin{grid_view}
  {
    Impl::UncoloredOverlapMixin<GV>::updatePartitionOverlap(*this, halo);
  }

  /**
   * @brief Update partition with new patches
   *
   * @param patches     Number of patches to have in total
   * @param halo        Distance another entity in the same label set is considered connected
   */
  void update(std::size_t patches, std::size_t halo) {
    Impl::MetisMixin<GV, codim, pit>::update(patches);
    Impl::UncoloredOverlapMixin<GV>::updatePartitionOverlap(*this, halo);
  }
};

/**
 * @brief Colored partition set with entities split on several patches using METIS
 * @details This is a partition on an grid view where entities seeds of the grid
 * view are split using METIS and assigned to different patches.
 * Labels are colored so the memory region of every entity in the set is private.
 * This class models the EntitySetPartition concept.
 *
 * @tparam GV     Grid view to partition
 * @tparam codim  Codimension to partition
 * @tparam pit    Parallel partition to operate on
 */
template<Dune::Concept::GridView GridView, int codim = 0, PartitionIteratorType pit = PartitionIteratorType::All_Partition>
struct MetisColored
  : public Impl::ColoringAdaptor<Impl::MetisMixin<GridView, codim, pit>>
{

  /**
   * @brief Construct a Metis colored entiy set partition
   *
   * @param grid_view   The grid view to operate on
   * @param patches     Number of patches to have in total
   * @param halo        Distance another entity in the same label set is considered connected
   */
  explicit MetisColored(const GridView& grid_view, std::size_t patches, std::size_t halo)
   : Impl::ColoringAdaptor<Impl::MetisMixin<GridView, codim, pit>>{Impl::MetisMixin<GridView, codim, pit>{grid_view, patches}, halo}
  {}

  /**
   * @brief Update partition with new patches
   *
   * @param patches     Number of patches to have in total
   * @param halo        Distance another entity in the same label set is considered connected
   */
  void update(std::size_t patches, std::size_t halo) {
    Impl::ColoringAdaptor<Impl::MetisMixin<GridView, codim, pit>>::update(Impl::MetisMixin<GridView, codim, pit>{this->gridView(), patches}, halo);
  }
};


} // namespace Dune::PDELab::EntitySetPartition

#endif // DUNE_PDELAB_COMMON_PARTITION_METIS_HH
