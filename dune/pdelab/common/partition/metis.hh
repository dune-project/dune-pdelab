#ifndef DUNE_PDELAB_COMMON_PARTITION_METIS_HH
#define DUNE_PDELAB_COMMON_PARTITION_METIS_HH


#include <dune/pdelab/common/partition/halo/region.hh>
#include <dune/pdelab/common/partition/halo/colored.hh>
#include <dune/pdelab/common/partition/halo/uncolored.hh>

#include <dune/grid/concepts/entity.hh>
#include <dune/grid/concepts/gridview.hh>
#include <dune/grid/common/partitionset.hh>

#include <dune/common/iteratorfacades.hh>
#include <dune/common/metis.hh>

#if !HAVE_METIS
#error This file should only be included if metis is available
#endif

#include <vector>
#include <memory>

namespace Dune::PDELab::EntitySetPartition {

namespace Impl {

/**
 * @brief Container of entity seeds
 * @details This container stores entity seeds and converts them into rvalue entities on traversal.
 * Models a range of entities
 * @todo Improve documentation
 *
 * @tparam GridView  entity set to transform entity seeds into entities
 * @tparam Entity    Entity to "store" in the container
 */
template<Dune::Concept::GridView EntitySet, Dune::Concept::Entity Entity>
class EntitySeedContainer {
  using EntitySeed = typename Entity::EntitySeed;

  // Iterator over entity seeds
  class ConstIterator : public RandomAccessIteratorFacade<ConstIterator, Entity, Entity> {
    using typename RandomAccessIteratorFacade<ConstIterator, Entity, Entity>::DifferenceType;
  public:

    // Construct an iterator from a seed iterator and a entity set
    explicit ConstIterator(typename std::vector<EntitySeed>::const_iterator seed_it, const EntitySet& entity_set)
      : _entity_set{&entity_set}
      , _seed_it{std::move(seed_it)}
    {}

    // Default constructor
    explicit ConstIterator() {}

    // Iterator concept methods
    bool equals(const ConstIterator &other) const
    {
      return _seed_it == other._seed_it;
    }

    // Dereference to get an entity
    Entity dereference() const {
      return _entity_set->grid().entity(*_seed_it);
    }

    // Advance the iterator
    void increment() {
      ++_seed_it;
    }

    // Decrement the iterator
    void decrement(){
      --_seed_it;
    }

    // Advance the iterator by n positions
    void advance(DifferenceType n) {
      std::advance(_seed_it, n);
    }

    // Get the distance to another iterator
    DifferenceType distanceTo(const ConstIterator& other) const {
      return std::distance(_seed_it, other._seed_it);
    }

  private:
    EntitySet const * _entity_set;
    typename std::vector<EntitySeed>::const_iterator _seed_it;
  };

public:

  // Construct an empty container for entity seeds
  explicit EntitySeedContainer(const EntitySet& entity_set, std::vector<EntitySeed>&& seeds)
    : _entity_set{entity_set}
    , _seeds{std::move(seeds)}
  {}

  // Begin of the entity seed container
  [[nodiscard]] ConstIterator begin() const {
    return ConstIterator{_seeds.begin(), _entity_set};
  }

  // End of the entity seed container
  [[nodiscard]] ConstIterator end() const {
    return ConstIterator{_seeds.end(), _entity_set};
  }

  [[nodiscard]] std::size_t size() const {
    return _seeds.size();
  }

private:
  EntitySet _entity_set;
  std::vector<EntitySeed> _seeds;
};


/**
 * @brief Partition set with entities split on several patches using METIS
 * @details This is a partition on an entity set where entities seeds of the grid
 * view are split using METIS and assigned to different patches of one label set.
 *
 * @tparam ES     Grid view to partition
 */
template<Dune::Concept::GridView ES>
class MetisMixin {
public:
  //! Uderlying entity set
  using EntitySet = ES;
  //! Entity being partitioned
  using Element = typename EntitySet::template Codim<0>::Entity;
  //! Range of entities grouped by a patch
  using PatchSet = EntitySeedContainer<EntitySet, Element>;
  //! Range of patches grouped by a label
  using LabelSet = std::vector<PatchSet>;
  //! Range of labels
  using PartitionSet = std::array<LabelSet,1>;

  /**
   * @brief Construct a Simple partition
   *
   * @param entity_set  Grid view to operate on
   * @param patches Number of patches to split the entity set
   */
  explicit MetisMixin(const EntitySet& entity_set, std::size_t patches)
   : _entity_set{entity_set}
   , _patches{patches}
  {
    update(_entity_set);
  }

  //! Uderlying entity set
  [[nodiscard]] EntitySet entitySet() const noexcept { return _entity_set; }

  //! Begin of the partition set
  [[nodiscard]] auto begin() const {
    return _partition_set->begin();
  }

  //! End of the partition set
  [[nodiscard]] auto end() const {
    return _partition_set->end();
  }

protected:
  /**
   * @brief Sets a new set of patches
   * @details The number of patches will be set up anew
   *
   * @param patches  Number of patches in the set
   */
  void update(EntitySet entity_set) {

    _entity_set = entity_set;
    if (_patches == 0) {
      std::size_t entities_per_patch = 100;
      switch(EntitySet::dimension) {
        case 1: entities_per_patch = 125; break;
        case 2: entities_per_patch = 250; break;
        default: entities_per_patch = 500; break;
      }
      _patches = std::max<std::size_t>(1,entity_set.size(0) / entities_per_patch);
    }
    if (_patches == 1) {
      // all entities in one patch
      _partition_set = std::make_shared<PartitionSet>();

      // set entity seeds to each assigned patch
      std::vector<typename Element::EntitySeed> seeds;
      for (const auto& element : elements(_entity_set))
        seeds.emplace_back(element.seed());

      _partition_set = std::make_shared<PartitionSet>();
      (*_partition_set)[0] = LabelSet{PatchSet{_entity_set, std::move(seeds)}};
      return;
    }


    const std::size_t dimension = EntitySet::dimension;
    using idx_t = Dune::Metis::idx_t;
    // setup METIS parameters
    idx_t parts = _patches;                             // number of partitions
    idx_t ncommonnodes = 2;                             // number of nodes elements must have in common to be considered adjacent to each other
    idx_t edgecut = std::numeric_limits<idx_t>::max();  // will store number of edges cut by partition

    std::vector<idx_t> cells, nodes, element_part, node_part;
    std::vector<typename Element::EntitySeed> seeds;
    cells.reserve(_entity_set.size(0));
    nodes.reserve(_entity_set.size(dimension));

    // create graph of elements and vertices
    int vertices = 0;
    cells.push_back(vertices);
    for (const auto& element : elements(_entity_set)) {
      const auto& ref_element = referenceElement<double, dimension>(element.type());
      edgecut = std::min<idx_t>(edgecut, ref_element.size(1));
      vertices += ref_element.size(dimension);
      cells.push_back(vertices);
      seeds.emplace_back(element.seed());

      for (int k = 0; k != ref_element.size(dimension); ++k)
        nodes.push_back(_entity_set.indexSet().subIndex(element, k, dimension));
    }

    idx_t element_count = cells.size()-1;
    idx_t node_count = nodes.size();
    element_part.assign(element_count, 0);
    node_part.assign(node_count, 0);

    // actual partition of elements
#if (!HAVE_SCOTCH_METIS) || (HAVE_SCOTCH_METIS && SCOTCH_VERSION >= 7)
    // METIS 5+ supports dual graph partitioning of meshes, but in case we use scotch, it's only implemented in 7+
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
      nullptr,
      &edgecut,
      element_part.data(),
      node_part.data()
    );

    if (result != METIS_OK)
      DUNE_THROW(PartitionError, "Metis could not partition the grid");
#else
    DUNE_THROW(PartitionError, "Metis version does not support mesh partitioning, please update to metis 5+ or scotch 7+");
#endif

    // assign entity seeds to each patch
    std::vector<std::vector<typename Element::EntitySeed>> seed_patches(_patches);
    for (std::size_t i = 0; i != element_part.size(); ++i)
      seed_patches[element_part[i]].emplace_back(std::move(seeds[i]));

    // move seeds into patch sets
    _partition_set = std::make_shared<PartitionSet>();
    for (std::size_t i = 0; i != _patches; ++i)
      (*_partition_set)[0].emplace_back(_entity_set, std::move(seed_patches[i]));
  }

private:
  std::shared_ptr<PartitionSet> _partition_set;
  EntitySet _entity_set;
  std::size_t _patches;
};

} //namespace Impl


/**
 * @brief Partition set with entities split on several patches using METIS
 * @details This is a partition on an entity set where entities seeds of the grid
 * view are split using METIS and assigned to different patches of one label set.
 * The patches are neighboring each other, thus, inducing a halo overlap in some entities of the partition.
 *
 * @tparam ES     Grid view to partition
 */
template<Dune::Concept::GridView EntitySet>
struct Metis
  : public Impl::MetisMixin<EntitySet>
  , public Impl::UncoloredHaloMixin<EntitySet>
{
  /**
   * @brief Construct a Metis uncolored entity set partition
   *
   * @param entity_set      The entity set to operate on
   * @param patches         Number of patches to have in total (0 means automatic choice)
   * @param halo_distance   Neighbor distance where another entity in the same label set is considered connected
   */
  explicit Metis(const EntitySet& entity_set, std::size_t patches = 0, std::size_t halo_distance = all_interior_halo_region)
   : Impl::MetisMixin<EntitySet>{entity_set, patches}
   , Impl::UncoloredHaloMixin<EntitySet>::UncoloredHaloMixin{halo_distance}
  {
    update(entity_set);
  }

  //! Update the partition with a new entity set
  void update(const EntitySet& entity_set) {
    Impl::MetisMixin<EntitySet>::update(entity_set);
    Impl::UncoloredHaloMixin<EntitySet>::updateHalo(*this);
  }
};

/**
 * @brief Colored partition set with entities split on several patches using METIS
 * @details This is a partition on an entity set where entities seeds of the grid
 * view are split using METIS and assigned to different patches.
 * Labels are colored so the halo region of every entity in the set is interior.
 *
 * @tparam ES     Grid view to partition
 */
template<Dune::Concept::GridView EntitySet>
struct MetisColored
  : public Impl::ColoredHaloAdaptor<Impl::MetisMixin<EntitySet>>
{

  /**
   * @brief Construct a Metis colored entiy set partition
   *
   * @param entity_set      The entity set to operate on
   * @param patches         Number of patches to have in total
   * @param halo_distance   Distance another entity in the same label set is considered connected
   */
  explicit MetisColored(const EntitySet& entity_set, std::size_t patches = 0, std::size_t halo_distance = all_interior_halo_region)
    : Impl::ColoredHaloAdaptor<Impl::MetisMixin<EntitySet>>{Impl::MetisMixin<EntitySet>{entity_set, patches}, halo_distance}
    , _patches{patches}
  {}

  //! Update the partition with a new entity set
  void update(const EntitySet& entity_set) {
    Impl::MetisMixin<EntitySet> base(entity_set, _patches);
    Impl::ColoredHaloAdaptor<Impl::MetisMixin<EntitySet>>::update(std::move(base));
  }

private:
  std::size_t _patches;
};


} // namespace Dune::PDELab::EntitySetPartition

#endif // DUNE_PDELAB_COMMON_PARTITION_METIS_HH
