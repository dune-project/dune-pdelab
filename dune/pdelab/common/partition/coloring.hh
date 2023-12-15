#ifndef DUNE_PDELAB_COMMON_PARTITION_COLORING_HH
#define DUNE_PDELAB_COMMON_PARTITION_COLORING_HH

#include <dune/pdelab/common/partition/halo.hh>
#include <dune/pdelab/concepts/entityset_partition.hh>

#include <dune/grid/common/mcmgmapper.hh>

#include <set>
#include <vector>

namespace Dune::PDELab::inline Experimental::EntitySetPartitioner::Impl {

/**
 * @brief Adapts a partition with one label set into a colored partition
 * @details This class first creates a conectiviy graph based on the halo
 * parameter, and then uses a DSatur algorithm to give a color to each patch.
 * Finally, each patch is moved to its respective (colored) label set.
 * @todo Improve graph data structure to an adjacency list
 * @tparam BasePartition  An entity set partition with one label set
 */
template<class BasePartition>
class ColoringAdaptor : public ColoredHaloMixin {
public:
  //! Uderlying grid view
  using EntitySet = typename BasePartition::EntitySet;
  //! Entity being partitioned
  using Entity = typename BasePartition::Entity;
  //! Range of entities grouped by a patch
  using PatchSet = typename BasePartition::PatchSet;
  //! Range of patches grouped by a label
  using LabelSet = std::vector<PatchSet>;
  //! Range of labels
  using PartitionSet = std::vector<LabelSet>;

protected:
  void update(BasePartition&& base, std::size_t halo) {
    // pre-condition: there is only one label on the base partition
    // post-condition: patches in the same color do not have entities in their halo
   Dune::MultipleCodimMultipleGeomTypeMapper<EntitySet> mapper{base.entitySet(), [](GeometryType gt, int dimgrid) { return true; }};

    std::size_t no_id = std::numeric_limits<std::size_t>::max();
    // assign patch id to each entity
    std::vector<std::size_t> patch_id(mapper.size(), no_id);

    std::size_t label = 0;
    std::size_t patch = 0;
    for (const auto& label_set : base.range()) {
      for (const auto& patch_set : label_set) {
        for (const auto& entity : patch_set) {
          patch_id[mapper.index(entity)] = patch;
        }
        ++patch;
      }
      ++label;
    }

    // check pre-condition
    if (label != 1)
      DUNE_THROW(InvalidStateException, "Partition should only contain one label at setup");

    auto get_patch_id = [&](const Entity& entity) -> std::size_t {
      return patch_id[mapper.index(entity)];
    };

    std::vector<std::set<std::size_t>> connectivity(patch);

    auto add_entity_link = [&](const Entity& entity_in, const Entity& entity_out) {
      connectivity[get_patch_id(entity_in)].insert(get_patch_id(entity_out));
      connectivity[get_patch_id(entity_out)].insert(get_patch_id(entity_in));
    };

    std::function<void(const Entity&, const Entity&, std::size_t)> add_link;
    add_link = [&](const Entity& entity_in, const Entity& entity_out, std::size_t current_halo) {
      add_entity_link(entity_in, entity_out);
      if ((current_halo--) != 0)
        for (const auto& intersection : intersections(base.entitySet(), entity_out))
          if (intersection.neighbor())
            add_link(entity_in, intersection.outside(), current_halo);
    };

    // make patch conectivity graph according to halo
    for (const auto& label_set : base.range())
      for (const auto& patch_set : label_set)
        for (const auto& entity : patch_set)
          add_link(entity, entity, halo+1);

    // perfom actual coloring
    auto [colors, color_count] = make_colors(connectivity);

    // debug: print vertex index, color, and neighbors
    // for (std::size_t i = 0; i != connectivity.size(); ++i) {
    //   std::cout << i << ", " << colors[i] << ": [";
    //   for (auto j : connectivity[i])
    //     std::cout << j << ", ";
    //   std::cout << "]" << std::endl;
    // }

    // assing colored partitions
    _partition_set->assign(color_count, LabelSet{});
    patch = 0;
    for (auto& label_set : base.range())
      for (auto& patch_set : label_set)
        (*_partition_set)[colors[patch++]].emplace_back(std::move(patch_set));
  }

  // use DSatur to color graph (https://en.wikipedia.org/wiki/DSatur#Pseudocode)
  auto make_colors(const std::vector<std::set<std::size_t>>& graph, bool verify = false) {
    const std::size_t uncolored = std::numeric_limits<std::size_t>::max();
    std::vector<std::size_t> colors(graph.size(), uncolored);
    std::vector<std::size_t> saturation(graph.size(), 0);

    // stores the vertex id of the saturation in a binary heap
    std::vector<std::size_t> sat_heap;

    // saturation definition
    auto sat_comp = [&](auto lhs, auto rhs) {
      auto comp = saturation[lhs] <=> saturation[rhs];
      if (comp == 0)
        return graph[lhs].size() > graph[rhs].size();
      else
        return comp < 0;
    };

    { // find starting point: maximum saturation
      std::size_t max_deg = 0;
      std::size_t saturated_vertex = 0;
      for (std::size_t vertex = 0; vertex != graph.size(); ++vertex)
        if (max_deg < graph[vertex].size()) {
          max_deg = graph[vertex].size();
          saturated_vertex = vertex;
        }
      sat_heap.push_back(saturated_vertex);
    }

    std::vector<std::size_t> neighbor_colors;
    std::size_t max_color = 0;
    while (not sat_heap.empty()) {
      std::pop_heap(sat_heap.begin(), sat_heap.end(), sat_comp);
      std::size_t saturated_vertex = sat_heap.back();
      sat_heap.pop_back();

      neighbor_colors.clear();
      for(auto neighbor : graph[saturated_vertex]) {
        if (saturated_vertex == neighbor) continue; // skip self
        neighbor_colors.push_back(colors[neighbor]);
        saturation[neighbor] += 1;
        if (colors[neighbor] == uncolored) {
          sat_heap.push_back(neighbor);
          std::push_heap(sat_heap.begin(), sat_heap.end(), sat_comp);
        }
      }

      // choose lowest color not used
      std::sort(neighbor_colors.begin(), neighbor_colors.end());
      std::size_t color_candidate = 0;
      for(auto neighbor_color : neighbor_colors)
        if (color_candidate == neighbor_color)
          color_candidate = neighbor_color+1;

      colors[saturated_vertex] = color_candidate;
      max_color = std::max(max_color, color_candidate);
    }

    // check post-condition: vertices do not neighbor any vertex with the same color
    if (verify)
      for (std::size_t vertex = 0; vertex != graph.size(); ++vertex)
        for(auto neighbor : graph[vertex])
          if (vertex != neighbor and colors[vertex] == colors[neighbor])
            DUNE_THROW(MathError, "Coloring is incorrect!");

    return std::make_tuple(colors, max_color+1);
  }

public:

  /**
 * @brief Construct a colored partition using base partition patches
 *
 * @param base_partition  Base partition with patches to color
 * @param halo            Distance another entity in the same label set is considered connected
 */
  explicit ColoringAdaptor(BasePartition&& base_partition, std::size_t halo)
    : _entity_set{base_partition.entitySet()}
    , _partition_set{std::make_shared<PartitionSet>()}
  {
    update(std::move(base_partition), halo);
  }

  //! Uderlying grid view
  [[nodiscard]] EntitySet entitySet() const { return _entity_set; }

  //! Range of the partition set
  [[nodiscard]] const PartitionSet& range() const noexcept { return *_partition_set; }

private:
  EntitySet _entity_set;
  std::shared_ptr<PartitionSet> _partition_set;
};

} // namespace Dune::PDELab::EntitySetPartitioner::Impl

#endif // DUNE_PDELAB_COMMON_PARTITION_COLORING_HH
