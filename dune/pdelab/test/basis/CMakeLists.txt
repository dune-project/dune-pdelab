add_subdirectory(constraints)
add_subdirectory(ordering)
add_subdirectory(wrapper)

add_custom_target(test-basis-ordering)

# make a test for each of the basis fixture list
foreach(fixture LeafBasisFixtures ArrayBasisFixtures VectorBasisFixtures TupleBasisFixtures)

  string(TOLOWER ${fixture} fixture_id)

  dune_add_test(
    NAME test-basis-ordering-${fixture_id}
    COMPILE_FLAGS -DBASIS_FIXTURES=${fixture}
    SOURCES test-ordering.cc)
  gtest_discover_tests(test-basis-ordering-${fixture_id})
  add_dependencies(test-basis-ordering test-basis-ordering-${fixture_id})

endforeach()
