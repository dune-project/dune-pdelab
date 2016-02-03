# .. cmake_module::
#
#    Module that checks for C++14 support.
#
#    This module internally sets the following variables, which are then
#    exported into the config.h of the current dune module.
#
#    :code:`PDELAB_HAVE_CXX14`
#       The module sets this variable to TRUE if your compiler supports C++14.
#       You can use this variable in your CMakeLists.txt files to guard targets
#       that should only be built in C++14 mode.
#
# .. cmake_function:: pdelab_require_cxx14
#
#    This function checks whether your compiler supports C++14 and aborts with
#    an error if it doesn't.
#

include(CMakePushCheckState)

if(CXX_FLAG_CXX14 AND CXX_LIB_SUPPORTS_CXX14)
  set(PDELAB_HAVE_CXX14 TRUE)
else()
  set(PDELAB_HAVE_CXX14 FALSE)
endif()

function(pdelab_require_cxx14)
  if(NOT PDELAB_HAVE_CXX14)
    message(FATAL_ERROR "This module requires a compiler with support for C++14, which you compiler does not have.")
  endif()
endfunction()
