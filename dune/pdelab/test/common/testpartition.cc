#include <dune/pdelab/common/partition/iteratorsplit.hh>
#include <dune/pdelab/common/partition/identity.hh>
#if HAVE_METIS
#include <dune/pdelab/common/partition/metis.hh>
#endif

#include <dune/grid/io/file/vtk.hh>

#include <dune/grid/uggrid.hh>

#include <dune/common/test/testsuite.hh>
#include <dune/common/classname.hh>

#include "../gridexamples.hh"

#include <mutex>


using namespace Dune::PDELab::EntitySetPartition;

template<class EntitySet> std::string partitionName(IteratorSplit<EntitySet>)        { return "IteratorSplit"; }
template<class EntitySet> std::string partitionName(IteratorSplitColored<EntitySet>) { return "IteratorSplitColored"; }
#if HAVE_METIS
template<class EntitySet> std::string partitionName(Metis<EntitySet>)                { return "Metis"; }
template<class EntitySet> std::string partitionName(MetisColored<EntitySet>)         { return "MetisColored"; }
#endif
template<class EntitySet> std::string partitionName(Identity<EntitySet>)             { return "Identity"; }

std::string haloName(std::size_t halo_distance) {
  if (halo_distance == all_interior_halo_region)
    return "i";
  else if (halo_distance == all_overlap_halo_region)
    return "o";
  else
    return std::to_string(halo_distance);
}

template<class Partition>
void testPartitionOverlapRegion(Dune::TestSuite& test, Partition partition, std::size_t halo_distance) {

   if (halo_distance == Dune::PDELab::EntitySetPartition::all_interior_halo_region or halo_distance == Dune::PDELab::EntitySetPartition::all_overlap_halo_region) {
     for (const auto& label_set : partition) {
       for (const auto& patch_set : label_set) {
         for (const auto& entity : patch_set) {
            if (halo_distance == Dune::PDELab::EntitySetPartition::all_overlap_halo_region)
              test.check(partition.haloRegion(entity) == Dune::PDELab::EntitySetPartition::HaloRegion::Overlap,
                        "All entities should be in the overlap region");
            else
              test.check(partition.haloRegion(entity) == Dune::PDELab::EntitySetPartition::HaloRegion::Interior,
                        "All entities should be in the interior region");
         }
       }
     }
     return;
   }

    using Element = typename Partition::Element;
    using EntitySet = typename Partition::EntitySet;
    using mo = std::memory_order;


    auto all_entities_layout = [](Dune::GeometryType gt, int dimgrid) { return true; };
    Dune::MultipleCodimMultipleGeomTypeMapper<EntitySet> all_mapper(partition.entitySet(), all_entities_layout);
    Dune::MultipleCodimMultipleGeomTypeMapper<EntitySet> element_mapper(partition.entitySet(), Dune::mcmgElementLayout());

    // temporary storage where every entry can be updated concurrently
    // note that proxy objects in std::vector<bool> are not thread-safe
    std::vector<std::atomic<std::size_t>> entity_owner;
    std::vector<std::atomic<bool>> entity_in_overlap;

    // mark every entity with a unique owner
    auto mark_owner = [&](const Element& entity, std::size_t id, auto) {
      Dune::Hybrid::forEach(std::make_index_sequence<Element::dimension+1>{}, [&](auto codim){
        for (const auto& sub_entity : subEntities(entity, Dune::Codim<codim>{})) {
          auto entity_index = all_mapper.index(sub_entity);
          entity_owner[entity_index].store(id, mo::relaxed);
        }
      });
    };

    // if an entity is owned by different patches, it is in the overlap
    auto mark_entity_overlap = [&](const Element& entity, std::size_t id) {
      Dune::Hybrid::forEach(std::make_index_sequence<Element::dimension+1>{}, [&](auto codim){
        for (const auto& sub_entity : subEntities(entity, Dune::Codim<codim>{})) {
          auto entity_index = all_mapper.index(sub_entity);
          auto owner = entity_owner[entity_index].load(mo::relaxed);
          if (owner != id)
            entity_in_overlap[entity_index].store(true, mo::relaxed);
        }
      });
    };

    // propagate the overlap on neighboring entities recursively until the halo distance is reached
    std::function<void(const Element&,std::size_t,std::size_t)> mark_overlap;
    mark_overlap = [&](const Element& entity, std::size_t id, std::size_t halo_distance) {
      mark_entity_overlap(entity, id);
      if (halo_distance != 0) {
        for (const auto& intersection : intersections(partition.entitySet(), entity))
          if (intersection.neighbor())
            mark_overlap(intersection.outside(), id, halo_distance-1);
      }
    };

    [[maybe_unused]] std::mutex testmutex;
    // if any of the sub-entities is in the overlap, the entity is in the halo and should be private
    auto check_overlap = [&](const Element& entity, auto, auto) {
      Dune::Hybrid::forEach(std::make_index_sequence<Element::dimension+1>{}, [&](auto codim){
        for (const auto& sub_entity : subEntities(entity, Dune::Codim<codim>{})) {
          auto sub_entity_index = all_mapper.index(sub_entity);
          bool in_overlap = entity_in_overlap[sub_entity_index].load(mo::relaxed);
          if (in_overlap) {
#if __cpp_lib_execution >= 201603L
            std::lock_guard lock{testmutex};
#endif
            test.check(partition.haloRegion(entity) == Dune::PDELab::EntitySetPartition::HaloRegion::Overlap, "Entity in overlap should be private");
          }
        }
      });
    };

    // helper to run a function on all entities of all patches in parallel
    auto for_each_element = [halo_distance = halo_distance](const auto& entity_sets, auto apply){
      auto patches = std::distance(entity_sets.begin(), entity_sets.end());
      auto partitions = Dune::range(std::size_t{0}, static_cast<std::size_t>(patches));
      std::for_each(
#if __cpp_lib_execution >= 201603L
        std::execution::par,
#endif
        partitions.begin(), partitions.end(), [&](auto patch_id){
        for (const auto& entity : entity_sets[patch_id])
          apply(entity, patch_id, halo_distance);
      });
    };

    for (const auto& concurrent_entity_sets : partition) {
      // clean up overlap info
      entity_owner = std::vector<std::atomic<std::size_t>>(all_mapper.size());
      entity_in_overlap = std::vector<std::atomic<bool>>(all_mapper.size());
      for (std::size_t i=0; i<all_mapper.size(); i++) {
        entity_owner[i].store(std::numeric_limits<std::size_t>::max(), mo::relaxed);
        entity_in_overlap[i].store(false, mo::relaxed);
      }

      // phase 1: assign an owner to each used sub-entity
      for_each_element(concurrent_entity_sets, mark_owner);

      // phase 2: mark the overlap based in incompatible owners
      for_each_element(concurrent_entity_sets, mark_overlap);

      // phase 3: entities with an overlap are marked as shared region
      for_each_element(concurrent_entity_sets, check_overlap);
    }
}

template<class Partition>
void writePartition(Partition partition, std::string name, std::size_t patches, std::size_t halo_distance) {

  Dune::VTKWriter<typename Partition::EntitySet> vtkwriter(partition.entitySet());

  std::vector<std::size_t> entity_patche(partition.entitySet().size(0));
  std::vector<std::size_t> entity_label(partition.entitySet().size(0));
  std::vector<std::size_t> entity_region(partition.entitySet().size(0));
  Dune::MultipleCodimMultipleGeomTypeMapper<typename Partition::EntitySet> element_mapper(partition.entitySet(), Dune::mcmgElementLayout());

  std::size_t patch_id = 0;
  std::size_t label_id = 0;
  for (const auto& patches : partition) {
    for (const auto& patch : patches) {
      for (const auto& element : patch) {
        auto index = element_mapper.index(element);
        entity_patche[index] = patch_id;
        entity_label[index] = label_id;
        entity_region[index] = (partition.haloRegion(element) == Dune::PDELab::EntitySetPartition::HaloRegion::Overlap);
      }
      patch_id++;
    }
    label_id++;
  }

  vtkwriter.addCellData(entity_patche, "patch");
  vtkwriter.addCellData(entity_label, "label");
  vtkwriter.addCellData(entity_region, "halo");
  auto key = partitionName(partition) + "_" + name;

  vtkwriter.write("partition-" + key + "_p" + std::to_string(patches) + "_h" + haloName(halo_distance), Dune::VTK::ascii);
}

template<class Partition>
void testPatchedPartition(Dune::TestSuite& test, std::string name, auto entity_set, auto patches_range, auto halo_distance_range){
  for (std::size_t patches : patches_range) {
    for (std::size_t halo_distance : halo_distance_range) {
      auto start_t = std::chrono::high_resolution_clock::now();

#if HAVE_METIS && HAVE_SCOTCH_METIS && (SCOTCH_VERSION < 7)
      if constexpr (std::is_same_v<Partition, Metis<typename Partition::EntitySet>> or std::is_same_v<Partition, MetisColored<typename Partition::EntitySet>>) {
        if (patches == 0)
          continue; // automatic choice may choose 1 patch and thus not trigger the expected error below
        if (patches != 1) {
          test.checkThrow<PartitionError>([&]{
            Partition partition(entity_set, patches, halo_distance);
          }, "Metis version checks")
            << "Expected an 'PartitionError' to be thrown on construction of " << Dune::className<Partition>() << " " << name << " [" << patches << ", " << haloName(halo_distance) << "]" <<std::endl;
          continue;
        }
      }
#endif

      if (halo_distance == all_overlap_halo_region) {
        if constexpr (std::is_same_v<Partition, IteratorSplitColored<typename Partition::EntitySet>>
#if HAVE_METIS
                      or std::is_same_v<Partition, MetisColored<typename Partition::EntitySet>>
#endif
                    ) {
          test.checkThrow<PartitionError>([&]{
            Partition partition(entity_set, patches, halo_distance);
          }, "Coloring with all overlap") << "Expected an 'PartitionError' to be thrown";
          continue;
        }
      }

      Partition partition(entity_set, patches, halo_distance);
      auto time = std::chrono::high_resolution_clock::now() - start_t;
      std::cout  << std::chrono::duration<double>(time).count() << " s to build " << partitionName(partition) << " " << name << " [" << patches << ", " << haloName(halo_distance) << "]" << std::endl;
      testPartitionOverlapRegion(test, partition, halo_distance);
      writePartition(partition, name, patches, halo_distance);
    }
  }

};

template<class GridView>
void testPartitions(Dune::TestSuite& test, std::string name, GridView gv){
  auto makeArray = [](auto... vals){ return std::array{static_cast<std::size_t>(vals)...}; };
  auto patches_range = makeArray(0, 1, 20);
  auto halo_distance_range = makeArray(0, 1, all_interior_halo_region, all_overlap_halo_region);
  testPatchedPartition<IteratorSplit<GridView>>         (test, name, gv, patches_range, halo_distance_range);
  testPatchedPartition<IteratorSplitColored<GridView>>  (test, name, gv, patches_range, halo_distance_range);
#if HAVE_METIS
  testPatchedPartition<Metis<GridView>>                 (test, name, gv, patches_range, halo_distance_range);
  testPatchedPartition<MetisColored<GridView>>          (test, name, gv, patches_range, halo_distance_range);
#endif

  Dune::PDELab::EntitySetPartition::Identity<GridView> identity(gv);
  testPartitionOverlapRegion(test, identity, Dune::PDELab::EntitySetPartition::all_interior_halo_region);
}

int main(int argc, char **argv) {

  Dune::MPIHelper::instance(argc, argv);

  Dune::TestSuite test("PartitionTestSuite");

  {
    YaspUnitSquare grid;
    grid.globalRefine(5);
    testPartitions(test, "YaspUnitSquare", grid.leafGridView());
  }

  {
    auto grid = TriangulatedLDomainMaker<Dune::UGGrid<2>>::create();
    grid->globalRefine(4);
    testPartitions(test, "TriangulatedLDomain", grid->leafGridView());
  }

  {
    auto grid = UnitTetrahedronMaker<Dune::UGGrid<3>>::create();
    grid->globalRefine(4);
    testPartitions(test, "UnitTetrahedron", grid->leafGridView());
  }

  return test.exit();
}
