find_package(Eigen3)

set(HAVE_EIGEN ${EIGEN3_FOUND})

if(EIGEN3_FOUND)
  include_directories(${EIGEN3_INCLUDE_DIR})
endif(EIGEN3_FOUND)
