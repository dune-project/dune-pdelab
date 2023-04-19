if (NOT GTest_FOUND)
  find_package(GTest)

  if (NOT GTest_FOUND)
    include(FetchContent)
    FetchContent_Declare(
      googletest
      GIT_REPOSITORY https://github.com/google/googletest.git
      GIT_TAG v1.13.0
    )

    FetchContent_MakeAvailable(googletest)
    include(GoogleTest)

    set_target_properties(gtest gtest_main gmock gmock_main PROPERTIES EXCLUDE_FROM_ALL TRUE)
    if(${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.18)
      set(CMAKE_GTEST_DISCOVER_TESTS_DISCOVERY_MODE PRE_TEST)
    endif()
  endif()
endif()
