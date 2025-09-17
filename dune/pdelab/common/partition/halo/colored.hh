#ifndef DUNE_PDELAB_COMMON_PARTITION_HALO_COLORED_HH
#define DUNE_PDELAB_COMMON_PARTITION_HALO_COLORED_HH

#include <dune/pdelab/common/partition/halo/region.hh>

#include <dune/istl/matrixindexset.hh>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/concepts/entity.hh>

#include <set>
#include <span>
#include <vector>


namespace Dune::PDELab::EntitySetPartition {

  //! Exception for coloring errors
struct ColoringError : public PartitionError {};

namespace Impl {

/**
 * @brief Adapts a partition with patches into a colored partition
 * @details This class first creates a conectiviy graph based on the halo distance
 * and then uses a DSatur algorithm to give a color to each patch.
 * Finally, each patch is moved to its respective (colored) label set.
 * @tparam BasePartition  A partition of an entity set
 */
template<class BasePartition>
class ColoredHaloAdaptor {
public:
  //! Uderlying entity set
  using EntitySet = typename BasePartition::EntitySet;
  //! Entity being partitioned
  using Element = typename BasePartition::Element;
  //! Range of entities grouped by a patch
  using PatchSet = typename BasePartition::PatchSet;
  //! Range of patches grouped by a label
  using LabelSet = std::vector<PatchSet>;
  //! Range of labels
  using PartitionSet = std::vector<LabelSet>;


  /**
   * @brief Construct a colored partition using base partition patches
   *
   * @param base               Base partition to operate on
   * @param halo_distance      Distance another entity in the same label set is considered connected
   */
  explicit ColoredHaloAdaptor(BasePartition&& base, std::size_t halo_distance)
    : _entity_set{base.entitySet()}
    , _halo_distance{halo_distance}
  {
    update(std::move(base));
  }

  //! Halo region of the partition for an entity
  [[nodiscard]] constexpr static auto haloRegion(const Dune::Concept::Entity auto& entity) {
    return interior_halo_region;
  }

  //! Underlying entity set
  [[nodiscard]] EntitySet entitySet() const noexcept { return _entity_set; }

  //! begin of the partition set range
  [[nodiscard]] auto begin() const { return _partition_set->begin(); }

  //! end of the partition set range
  [[nodiscard]] auto end() const { return _partition_set->end(); }

private:

  // A simple graph representation using CSR format
  class Graph {
  public:
    Graph(std::vector<std::size_t>&& row_ptr, std::vector<std::size_t>&& col_idx)
      : row_ptr_(std::move(row_ptr)), col_idx_(std::move(col_idx)) {}

    auto operator[](std::size_t vertex) const {
      return std::span(col_idx_.data() + row_ptr_[vertex], col_idx_.data() + row_ptr_[vertex + 1]);
    }

    auto size() const { return row_ptr_.size() - 1; }
  private:
    std::vector<std::size_t> row_ptr_;
    std::vector<std::size_t> col_idx_;
  };

protected:

  void update(BasePartition&& base) {
    // post-condition: patches in the same color do not have entities in their halo
    _entity_set = base.entitySet();

    if (_halo_distance == all_overlap_halo_region) {
      DUNE_THROW(PartitionError, "ColoredHaloAdaptor does not support all-overlap halo region.");
    } else if (_halo_distance == all_interior_halo_region) {
      // all patches must be in one label
      _partition_set = std::make_shared<PartitionSet>(1);
      for (auto& label_set : base)
        for (auto& patch_set : label_set)
          (*_partition_set)[0].emplace_back(std::move(patch_set));
      return;
    }

    // build connectivity graph between patches
    Graph graph = makeGraph(base, _halo_distance);

    // perfom actual coloring
    auto [colors, color_count] = makeColors(graph);

    // assing colored partitions
    _partition_set = std::make_shared<PartitionSet>(color_count);
    std::size_t patch = 0;
    for (auto& label_set : base)
      for (auto& patch_set : label_set)
        (*_partition_set)[colors[patch++]].emplace_back(std::move(patch_set));
  }

private:

  // build connectivity graph between patches based on shared entities in their halo
  Graph makeGraph(const BasePartition& base, std::size_t halo_distance) {

    std::size_t patch = 0;
    for (auto& label_set : base)
      patch += std::distance(label_set.begin(), label_set.end());

    auto all_entities_layout =  [](GeometryType gt, int dimgrid) { return true; };
    Dune::MultipleCodimMultipleGeomTypeMapper<EntitySet> all_mapper(base.entitySet(), all_entities_layout);

    // build connectivity between patches based on shared entities in their halo
    Dune::MatrixIndexSet patch_links(patch, patch);
    Dune::MatrixIndexSet entity_owners(all_mapper.size(), patch);

    // if an entity is owned by different patches, it is in the overlap and a link between patches is made
    auto mark_entity_overlap = [&](const Element& entity, std::size_t id) {
      Hybrid::forEach(std::make_index_sequence<Element::dimension+1>{}, [&](auto codim){
        for (const auto& sub_entity : subEntities(entity, Dune::Codim<codim>{})) {
          auto entity_index = all_mapper.index(sub_entity);
          std::visit([&](const auto& owners){
            for (auto owner : owners) {
              patch_links.add(owner,id);
              patch_links.add(id,owner);
            }
          }, entity_owners.columnIndices(entity_index));
          entity_owners.add(entity_index, id);
        }
      });
    };

    // propagate the overlap on neighboring entities recursively until the halo distance is reached
    std::function<void(const Element&,std::size_t,std::size_t)> mark_overlap;
    mark_overlap = [&](const Element& entity, std::size_t id, std::size_t halo_distance) {
      mark_entity_overlap(entity, id);
      if (halo_distance != 0) {
        for (const auto& intersection : intersections(base.entitySet(), entity))
          if (intersection.neighbor())
            mark_overlap(intersection.outside(), id, halo_distance-1);
      }
    };

    patch = 0;
    for (const auto& label_set : base) {
      for (const auto& patch_set : label_set) {
        for (const auto& entity : patch_set) {
          mark_overlap(entity, patch, _halo_distance);
        }
        ++patch;
      }
    }

    // build graph from connectivity
    std::vector<std::size_t> row_ptr(patch_links.rows() + 1, 0);
    std::vector<std::size_t> col_idx;
    col_idx.reserve(patch_links.size());

    for (std::size_t i = 0; i != patch_links.rows(); ++i) {
      std::visit([&](const auto& neighbors){
        row_ptr[i+1] = row_ptr[i] + neighbors.size();
        col_idx.insert(col_idx.end(), neighbors.begin(), neighbors.end());
      }, patch_links.columnIndices(i));
    }
    return Graph{std::move(row_ptr), std::move(col_idx)};
  }

  // use DSatur to color graph (https://en.wikipedia.org/wiki/DSatur#Pseudocode)
  auto makeColors(const auto& graph) {
    constexpr std::size_t uncolored = std::numeric_limits<std::size_t>::max();
    constexpr std::size_t max_degree = std::numeric_limits<unsigned char>::max();
    std::vector<std::size_t> color(graph.size(), uncolored);
    std::vector<unsigned char> degree(graph.size(), 0);
    Dune::MatrixIndexSet adjacent_colors(graph.size(), graph.size());
    std::vector<bool> used(graph.size(), false);
    std::size_t max_color = 0;

    // // Struct to store information
    std::set<std::array<std::size_t,3>, std::greater<>> queue;

    // initialize the priority queue with saturation degree 0 and degree of each vertex
    for (std::size_t vertex = 0; vertex != graph.size(); ++vertex) {
      if (graph[vertex].size() >= max_degree)
        DUNE_THROW(ColoringError,
          "Vertex " << vertex << " has degree " << graph[vertex].size()
          << " which exceeds the maximum supported degree of " << max_degree << " for coloring.");
      degree[vertex] = graph[vertex].size();
      queue.emplace(std::array<std::size_t, 3>{ 0, degree[vertex], vertex });
    }

    // process all vertices in the queue
    while (not queue.empty()) {
      // choose the vertex with highest saturation degree, breaking ties with its degree
      // and remove it from the priority queue
      auto it = queue.begin();
      std::size_t node = (*it)[2];
      queue.erase(it);

      // mark used colors by neighbors
      for (std::size_t neighbor : graph[node])
        if (color[neighbor] != uncolored)
          used[color[neighbor]] = true;
      // identify the lowest feasible colour for this node
      std::size_t this_color = 0;
      while (this_color != used.size() && used[this_color])
        ++this_color;
      // reset used colors by neighbors
      for (std::size_t neighbor : graph[node])
        if (color[neighbor] != uncolored)
          used[color[neighbor]] = false;

      // assign color to node
      color[node] = this_color;
      max_color = std::max(max_color, this_color);

      // update priority queue: neighbor saturation and degrees change with the coloring of this node
      for (std::size_t neighbor : graph[node]) {
        if (color[neighbor] == uncolored) {
          auto it = queue.find( std::array<std::size_t, 3>{ adjacent_colors.rowsize(neighbor), degree[neighbor], neighbor });
          auto suggest_it = std::next(it);
          queue.erase(it);
          adjacent_colors.add(neighbor, this_color);
          degree[neighbor]--;
          queue.emplace_hint( suggest_it, std::array<std::size_t, 3>{ adjacent_colors.rowsize(neighbor), degree[neighbor], neighbor });
        }
      }
    }

    // check post-condition: vertices do not neighbor any vertex with the same color
    for (std::size_t vertex = 0; vertex != graph.size(); ++vertex)
      for(auto neighbor : graph[vertex]) {
        if (vertex != neighbor and color[vertex] == color[neighbor])
          DUNE_THROW(ColoringError, "Coloring is incorrect!");
        if (color[vertex] == uncolored)
          DUNE_THROW(ColoringError, "Vertex " << vertex << " is uncolored!");
        if (color[vertex] >= max_color+1)
          DUNE_THROW(ColoringError, "Vertex " << vertex << " is an unexpected color!");
      }

    return std::make_tuple(std::move(color), max_color+1);
  }

private:
  EntitySet _entity_set;
  std::shared_ptr<PartitionSet> _partition_set;
  std::size_t _halo_distance;
};

} // namespace Impl
} // namespace Dune::PDELab::EntitySetPartition

#endif // DUNE_PDELAB_COMMON_PARTITION_HALO_COLORED_HH
