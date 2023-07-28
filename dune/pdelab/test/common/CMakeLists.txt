dune_add_test(
  SOURCES test-multiindex.cc
  LINK_LIBRARIES ${GTEST_MAIN_LIBRARIES})

dune_add_test(
  SOURCES test-tree.cc)

dune_add_test(
  SOURCES test-container-entry.cc
  LINK_LIBRARIES ${GTEST_MAIN_LIBRARIES})

dune_add_test(
  SOURCES test-for-each.cc
  LINK_LIBRARIES ${GTEST_MAIN_LIBRARIES})

dune_add_test(
  SOURCES test-property.cc
  LINK_LIBRARIES ${GTEST_MAIN_LIBRARIES})

dune_add_test(
  SOURCES test-property-tree.cc
  LINK_LIBRARIES ${GTEST_MAIN_LIBRARIES})
