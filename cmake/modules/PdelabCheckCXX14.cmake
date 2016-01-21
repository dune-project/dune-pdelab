# .. cmake_module::
#
#    Module that checks for C++14 support.
#
#    The behaviour of this module can be modified by the following variable:
#
#    :ref:`PDELAB_DISABLE_CXX_VERSION_CHECK`
#       Disable checking for std=c++14
#
#    This module internally sets the following variables, which are then
#    exported into the config.h of the current dune module.
#
#    :code:`PDELAB_HAVE_CXX14`
#       The module sets this variable to TRUE if your compiler supports C++14.
#       You can use this variable in your CMakeLists.txt files to guard targets
#       that should only be built in C++14 mode.
#
#
# .. cmake_variable:: PDELAB_DISABLE_CXX_VERSION_CHECK
#
#    You may set this variable to TRUE to disable checking for
#    std=c++14.
#
# .. cmake_function:: pdelab_require_cxx14
#
#    This function checks whether your compiler supports C++14 and aborts with
#    an error if it doesn't.
#

include(CMakePushCheckState)
cmake_push_check_state()

# test for C++14 flag
if(NOT PDELAB_DISABLE_CXX_VERSION_CHECK)
  # try to use compiler flag -std=c++14
  include(TestCXXAcceptsFlag)
  check_cxx_accepts_flag("-std=c++14" PDELAB_CXX_FLAG_CXX14)

  include(CheckCXXSourceCompiles)
  cmake_push_check_state()
  set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} -std=c++14 ")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++14 ")
  set(CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_MINSIZEREL} -std=c++14 ")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -std=c++14 ")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -std=c++14 ")
  check_cxx_source_compiles("
      #include <memory>

      int main() {
        std::make_unique<int>();
      }
    " PDELAB_CXX_LIB_SUPPORTS_CXX14)
  cmake_pop_check_state()
endif()

if(PDELAB_CXX_FLAG_CXX14 AND PDELAB_CXX_LIB_SUPPORTS_CXX14)
  set(PDELAB_CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} -std=c++14")
  set(PDELAB_CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 ")
  set(PDELAB_CXX_STD14_FLAGS "-std=c++14")
  set(PDELAB_HAVE_CXX14 TRUE)
else()
  set(PDELAB_HAVE_CXX14 FALSE)
endif()

cmake_pop_check_state()


function(pdelab_require_cxx14)
  if(NOT PDELAB_HAVE_CXX14)
    message(FATAL_ERROR "This module requires a compiler with support for C++14, which you compiler does not have.")
  endif()
endfunction()
