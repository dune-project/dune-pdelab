cmake_minimum_required(VERSION 3.6)

dune_add_test(
  SOURCES test-multiindex.cc
  LINK_LIBRARIES GTest::gtest_main)
