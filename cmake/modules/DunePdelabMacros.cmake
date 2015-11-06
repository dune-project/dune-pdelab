include(DunePdelabTestMacros)
include(UsePETSc)
include(UseEigen)

function(add_dune_petsc_flags)
  if(PETSC_FOUND)
    cmake_parse_arguments(ADD_PETSC "SOURCE_ONLY;OBJECT" "" "" ${ARGN})
    if(ADD_PETSC_SOURCE_ONLY)
      set(_prefix SOURCE)
      set(_source_only SOURCE_ONLY)
      include_directories(${PETSC_INCLUDES})
    else()
      set(_prefix TARGET)
      if(ADD_PETSC_OBJECT)
        set(_prefix TARGET)
      else(ADD_PETSC_OBJECT)
        foreach(_target ${ADD_PETSC_UNPARSED_ARGUMENTS})
          target_link_libraries(${_target} ${PETSC_LIBRARIES})
        endforeach(_target ${ADD_PETSC_UNPARSED_ARGUMENTS})
      endif(ADD_PETSC_OBJECT)
      include_directories(${PETSC_INCLUDES})
    endif()

    set_property(${_prefix} ${ADD_PETSC_UNPARSED_ARGUMENTS} APPEND PROPERTY COMPILE_DEFINITIONS ENABLE_PETSC ${PETSC_DEFINITIONS})
    if(NOT (ADD_PETSC_SOURCE_ONLY OR ADD_PETSC_OBJECT))
      set_property(${_prefix} ${ADD_PETSC_UNPARSED_ARGUMENTS} APPEND PROPERTY LINK_LIBRARIES ${PETSC_LIBRARIES})
    endif(NOT (ADD_PETSC_SOURCE_ONLY OR ADD_PETSC_OBJECT))

  endif(PETSC_FOUND)
endfunction(add_dune_petsc_flags)
