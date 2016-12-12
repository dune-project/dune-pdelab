# Variables used by this module:
#
# EIGEN3_INCLUDE_DIR  Path to the directory with the Eigen3 headers
#
# Variables defined by this module:

# EIGEN3_FOUND  True if Eigen3 found and usable

# look for Eigen3 CMake test
# if it is not shipped with the CMake installation, look in the directory of Eigen3 headers
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${EIGEN3_INCLUDE_DIR}/cmake)
find_package(Eigen3)

set(HAVE_EIGEN ${EIGEN3_FOUND})

if(EIGEN3_FOUND)
  dune_register_package_flags(COMPILE_DEFINITIONS "ENABLE_EIGEN=1"
                              INCLUDE_DIRS "${EIGEN3_INCLUDE_DIR}")
endif(EIGEN3_FOUND)

function(add_dune_eigen_flags)
  if(EIGEN3_FOUND)
    cmake_parse_arguments(ADD_EIGEN "SOURCE_ONLY;OBJECT" "" "" ${ARGN})
    if(ADD_EIGEN_SOURCE_ONLY)
      set(_prefix SOURCE)
    else()
      set(_prefix TARGET)
    endif()

    include_directories(${EIGEN3_INCLUDE_DIR})
    set_property(${_prefix} ${ADD_EIGEN_UNPARSED_ARGUMENTS}
      APPEND PROPERTY
      COMPILE_DEFINITIONS "ENABLE_EIGEN=1")
  endif(EIGEN3_FOUND)
endfunction(add_dune_eigen_flags)
