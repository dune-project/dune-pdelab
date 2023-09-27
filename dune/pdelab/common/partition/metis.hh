#ifndef DUNE_PDELAB_COMMON_PARTITION_METIS_HH
#define DUNE_PDELAB_COMMON_PARTITION_METIS_HH


#include <dune/pdelab/common/partition/halo.hh>
#include <dune/pdelab/common/partition/coloring.hh>
#include <dune/pdelab/common/trace.hh>

#include <dune/grid/concepts/entity.hh>
#include <dune/grid/concepts/gridview.hh>
#include <dune/grid/common/partitionset.hh>

#include <dune/common/iteratorfacades.hh>

#if !HAVE_METIS
#error This file should only be included if metis is available
#endif

#include <metis.h>

#include <vector>
#include <memory>

namespace Dune::PDELab::inline Experimental::EntitySetPartitioner {

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
template<Dune::Concept::GridView EntitySet, Dune::Concept::Entity Entity>
class EntitySeedContainer {
  using EntitySeed = typename Entity::EntitySeed;

  class ConstIterator : public RandomAccessIteratorFacade<ConstIterator, Entity, Entity> {
    using typename RandomAccessIteratorFacade<ConstIterator, Entity, Entity>::DifferenceType;
  public:

    explicit ConstIterator(typename std::vector<EntitySeed>::const_iterator seed_it, const EntitySet& entity_set)
      : _seed_it{std::move(seed_it)}
      , _entity_set{&entity_set}
    {}

    explicit ConstIterator() {}

    bool equals(const ConstIterator &other) const
    {
      return _seed_it == other._seed_it;
    }

    Entity dereference() const {
      return _entity_set->grid().entity(*_seed_it);
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
    EntitySet const * _entity_set;
  };

public:

  explicit EntitySeedContainer(const EntitySet& entity_set)
    : _entity_set{entity_set}
  {}

  [[nodiscard]] ConstIterator begin() const {
    return ConstIterator{_seeds.begin(), _entity_set};
  }

  [[nodiscard]] ConstIterator end() const {
    return ConstIterator{_seeds.end(), _entity_set};
  }

  std::vector<EntitySeed>& seeds() {
    return _seeds;
  }

private:
  std::vector<EntitySeed> _seeds;
  EntitySet _entity_set;
};


/**
 * @brief Partition set with entities split on several patches using METIS
 * @details This is a partition on an grid view where entities seeds of the grid
 * view are split using METIS and assigned to different patches of one label set.
 *
 * @tparam ES     Grid view to partition
 * @tparam codim  Codimension to partition
 * @tparam pit    Parallel partition to operate on
 */
template<Dune::Concept::GridView ES, int codim = 0, PartitionIteratorType pit = PartitionIteratorType::All_Partition>
class MetisMixin {
public:
  //! Uderlying grid view
  using EntitySet = ES;
  //! Entity being partitioned
  using Entity = typename EntitySet::template Codim<codim>::Entity;
  //! Range of entities grouped by a patch
  using PatchSet = EntitySeedContainer<EntitySet, Entity>;
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
   : _partition_set{std::make_shared<PartitionSet>()}
   , _entity_set{entity_set}
  {
    update(patches);
  }

  //! Uderlying grid view
  [[nodiscard]] EntitySet entitySet() const noexcept { return _entity_set; }

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
    TRACE_EVENT("dune", "EntitySet::updatePartition");
    if (patches == 0)
      DUNE_THROW(InvalidStateException, "There cannot be zero partitions");
    const std::size_t dimension = EntitySet::dimension;
    // setup METIS parameters
    idx_t parts = patches;                              // number of partitions
    idx_t ncommonnodes = 2;                             // number of nodes elements must have in common to be considered adjacent to each other
    idx_t edgecut = std::numeric_limits<idx_t>::max();  // will store number of edges cut by partition
    idx_t options[METIS_NOPTIONS];                      // use default values for random seed, output and coupling
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_NUMBERING] = 0;

    std::vector<idx_t> cells, nodes, element_part, node_part;
    std::vector<typename Entity::EntitySeed> seeds;
    cells.reserve(_entity_set.size(0));
    nodes.reserve(_entity_set.size(dimension));

    // create graph of elements and vertices
    int vertices = 0;
    cells.push_back(vertices);
    for (const auto& element : entities(_entity_set, Dune::Codim<codim>{}, Dune::partitionSet<pit>())) {
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

    LabelSet label_set(patches, PatchSet{_entity_set});

    // set entity seeds to each assigned patch
    for (std::size_t i = 0; i != element_part.size(); ++i)
      label_set[element_part[i]].seeds().emplace_back(std::move(seeds[i]));

    (*_partition_set)[0] = std::move(label_set);
  }

private:
  std::shared_ptr<PartitionSet> _partition_set;
  EntitySet _entity_set;
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
 * @tparam ES     Grid view to partition
 * @tparam codim  Codimension to partition
 * @tparam pit    Parallel partition to operate on
 */
template<Dune::Concept::GridView ES, int codim = 0, PartitionIteratorType pit = PartitionIteratorType::All_Partition>
struct Metis
  : public Impl::MetisMixin<ES, codim, pit>
  , public Impl::UncoloredHaloMixin<ES>
{

  /**
   * @brief Construct a Metis uncolored entiy set partition
   *
   * @param entity_set   The grid view to operate on
   * @param patches     Number of patches to have in total
   * @param halo        Distance another entity in the same label set is considered connected
   */
  explicit Metis(const ES& entity_set, std::size_t patches, std::size_t halo)
   : Impl::MetisMixin<ES, codim, pit>{entity_set, patches}
   , Impl::UncoloredHaloMixin<ES>::UncoloredHaloMixin{entity_set}
  {
    Impl::UncoloredHaloMixin<ES>::updatePartitionHalo(*this, halo);
  }

  /**
   * @brief Update partition with new patches
   *
   * @param patches     Number of patches to have in total
   * @param halo        Distance another entity in the same label set is considered connected
   */
  void update(std::size_t patches, std::size_t halo) {
    Impl::MetisMixin<ES, codim, pit>::update(patches);
    Impl::UncoloredHaloMixin<ES>::updatePartitionHalo(*this, halo);
  }
};

/**
 * @brief Colored partition set with entities split on several patches using METIS
 * @details This is a partition on an grid view where entities seeds of the grid
 * view are split using METIS and assigned to different patches.
 * Labels are colored so the memory region of every entity in the set is private.
 * This class models the EntitySetPartition concept.
 *
 * @tparam ES     Grid view to partition
 * @tparam codim  Codimension to partition
 * @tparam pit    Parallel partition to operate on
 */
template<Dune::Concept::GridView EntitySet, int codim = 0, PartitionIteratorType pit = PartitionIteratorType::All_Partition>
struct MetisColored
  : public Impl::ColoringAdaptor<Impl::MetisMixin<EntitySet, codim, pit>>
{

  /**
   * @brief Construct a Metis colored entiy set partition
   *
   * @param entity_set   The grid view to operate on
   * @param patches     Number of patches to have in total
   * @param halo        Distance another entity in the same label set is considered connected
   */
  explicit MetisColored(const EntitySet& entity_set, std::size_t patches, std::size_t halo)
   : Impl::ColoringAdaptor<Impl::MetisMixin<EntitySet, codim, pit>>{Impl::MetisMixin<EntitySet, codim, pit>{entity_set, patches}, halo}
  {}

  /**
   * @brief Update partition with new patches
   *
   * @param patches     Number of patches to have in total
   * @param halo        Distance another entity in the same label set is considered connected
   */
  void update(std::size_t patches, std::size_t halo) {
    Impl::ColoringAdaptor<Impl::MetisMixin<EntitySet, codim, pit>>::update(Impl::MetisMixin<EntitySet, codim, pit>{this->entitySet(), patches}, halo);
  }
};


} // namespace Dune::PDELab::EntitySetPartitioner

#endif // DUNE_PDELAB_COMMON_PARTITION_METIS_HH
