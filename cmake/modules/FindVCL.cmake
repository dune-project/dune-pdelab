# Variables used by this module:
#
# VCL_ROOT  Path to the directory with the VCL headers
#
# Variables defined by this module:

# VCL_FOUND        True if VCL found and usable
# VCL_INCLUDE_DIR  Path to the VCL include directories

# look for the header file "vectorclass.h"
find_path(VCL_INCLUDE_DIR
  NAMES "vectorclass.h"
  PATHS "${VCL_ROOT}" "${VCL_ROOT}/include"
  )

# compile a simple program using VCL
include(CMakePushCheckState)
include(CheckCXXSourceCompiles)
cmake_push_check_state()
set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${VCL_INCLUDE_DIR})
set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} -march=native -O3")
check_cxx_source_compiles("
#include \"vectorclass.h\"
int main() {
  Vec4i a(10,11,12,13);
  Vec4i b(20,21,22,23);

  Vec4i c = a + b;

  int d = 0;
  d += horizontal_add(a);
  d += horizontal_add(b);

  return 0;
}
" VCL_COMPILE_TEST)
cmake_pop_check_state()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  "VCL"
  DEFAULT_MSG
  VCL_INCLUDE_DIR
  VCL_COMPILE_TEST
  )

# set HAVE_VCL
set(HAVE_VCL ${VCL_FOUND})

# register all VCL related flags
if(VCL_FOUND)
  dune_register_package_flags(INCLUDE_DIRS "${VCL_INCLUDE_DIR}")
endif(VCL_FOUND)
