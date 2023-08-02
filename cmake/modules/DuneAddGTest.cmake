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
    option(BUILD_GMOCK "Builds the googlemock subproject" ON)
    set(BUILD_GMOCK OFF)
    option(INSTALL_GTEST "Enable installation of googletest. (Projects embedding googletest may want to turn this OFF.)" ON)
    set(INSTALL_GTEST OFF)
    include(GoogleTest)

    set_target_properties(${GTEST_BOTH_LIBRARIES} PROPERTIES EXCLUDE_FROM_ALL TRUE)
    if(${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.18)
      set(CMAKE_GTEST_DISCOVER_TESTS_DISCOVERY_MODE PRE_TEST)
    endif()
  endif()
endif()
