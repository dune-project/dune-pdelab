#ifndef DUNE_PDELAB_TEST_BASIS_ORDERING_HH
#define DUNE_PDELAB_TEST_BASIS_ORDERING_HH

#include <dune/pdelab/concepts/basis.hh>
#include <dune/pdelab/common/tree_traversal.hh>

#include <dune/common/exceptions.hh>

#include <gtest/gtest.h>

#include <unordered_map>
#include <set>
#include <sstream>
#include <ostream>
#include <mutex>


namespace Dune::PDELab::inline Experimental {

template<Concept::Basis Basis>
void test_ordering(const Basis& basis) {

  //!@note this test relies that the direction of the multi-index on the local
  // view is root-to-leaf as in dune-functions
  // (legacy dune-pdelab grid function spaces are reversed in the wrapper)

  bool failed = false;
  std::stringstream out;

  auto is_dense = [](auto){return true;};

  auto dim = basis.dimension();
  out << "Basis Dimension: " << dim << std::endl;

  using LocalView = typename Basis::LocalView;
  // using LocalIndexSet = typename Basis::LocalIndexSet;

  LocalView lbasis1 = basis.localView();
  [[maybe_unused]] LocalView lbasis2 = basis.localView();
  // LocalIndexSet lindex_set = basis.localIndexSet();

  using CI = typename Basis::MultiIndex;
  using SizeType = typename CI::value_type;
  std::unordered_map<CI, std::set<SizeType>> prefix_set;

  // loop over the whole range in the entity set
  for(auto&& entity : elements(basis.entitySet())) {

    // inform the local basis about the current entity
    lbasis1.bind(entity);

    // // only check locks if available and sync mark is not used
    // if constexpr (Concept::Lockable<LocalView>) {
    //   // locks are only block if the lock actually is needed
    //   if (lbasis1.memoryRegion() == Dune::Assembler::MemoryRegion::Shared) {
    //     lbasis2.bind(entity);
    //     assert(lbasis1.memoryRegion() == lbasis2.memoryRegion());
    //     // first local basiss tries to adquire the lock
    //     std::unique_lock lock1{lbasis1, std::try_to_lock};
    //     // this code should be sequential, so try_lock() should never fail!
    //     if (not lock1.owns_lock())
    //       DUNE_THROW(Dune::ParallelError,
    //         "Basiss should not fail to lock on purely sequential loops");
    //     // second basis also tries to lock
    //     std::unique_lock lock2{lbasis2, std::try_to_lock};
    //     // since first lock has not released his ownership something is wrong if second basis is also locked
    //     if (lock2.owns_lock())
    //       DUNE_THROW(Dune::ParallelError,
    //         "Two local basiss bound to the same entity should not be able to be locked at the same time");

    //     lbasis2.unbind();
    //   }
    // }

    // iterate over all the nodes on the basis
    Dune::PDELab::forEachLeafNode(lbasis1.tree(), [&](const auto& node, auto& path){
      // visit all degrees of freedom
      for (std::size_t dof = 0; dof < node.size(); ++dof) {
        // get the container index for the local degree of freedom
        CI prefix = node.index(dof);
        out << "  Index: " << prefix << std::endl;

        // loop all possible prefixes for the current container index
        while (prefix.size() != 0) {
          // get outer container index in the prefix
          SizeType block_index = back(prefix);
           // because indices are 0-index based, all indices should be less than the whole basis dimension
          failed = (dim < block_index);
          // reduce current prefix by one
          prefix.pop_back();
          // add block index to the prefix set
          prefix_set[prefix].insert(block_index);
          out << "    Prefix: " << prefix << std::endl;
          out << "      BlockIndex: " << block_index << std::endl;
          out << "      Current size: " << prefix_set[prefix].size()
                    << std::endl;
          out << "      Current range: [" << *prefix_set[prefix].begin()
                    << ", " << *prefix_set[prefix].rbegin() << "]" << std::endl;
          out << "      Size: " << basis.size(prefix) << std::endl;
        }
      }
    });

    lbasis1.unbind();
  }

  // we have collected all possible prefixes, time to test against basis.size(reverse(prefix))
  for (const auto &[prefix, set] : prefix_set) {
    assert(not set.empty());
    auto front = *set.begin();
    auto back = *set.rbegin();
    out << "Suffix: " << prefix << "   --->   Range: [" << front << ", "
              << (back + 1) << ")" << std::endl;

    // check if pattern is dense: all the sets for all prefixes should be full
    bool dense_pattern = true; // whether prefix sampled a dense pattern
    if (((back + 1) != set.size()) or (front != 0)) {
      dense_pattern = false;
      out << "  CI prefix '" << prefix
                << "' has no dense pattern: " << std::endl;
      for (auto i : set)
        out << "    * " << i << std::endl;
    }

    auto size = basis.size(prefix);
    if ((back + 1) != size) {
      dense_pattern = false;
      out << "  CI prefix has size '" << size
                << "' but it does not match size of sampled pattern '"
                << set.size() << "'" << std::endl;
    }
    if (not dense_pattern) {
      if (is_dense(prefix)) {
        failed |= true;
      } else {
        // if ordering is not dense, we at least check that we are not out of bounds
        failed |= (back + 1) > size;
      }
    }
  }

  SCOPED_TRACE(out.str());

  EXPECT_FALSE(failed) << "CI vs CI Sufix size information does not match";
}

} // end namespace Dune::PDELab::inline Experimental


#endif // DUNE_PDELAB_TEST_BASIS_ORDERING_HH
