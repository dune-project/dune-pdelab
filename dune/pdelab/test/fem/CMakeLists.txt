set(orders 0 1 2 3)
set(dimensions 2 3)
set(geometries simplex cube)

# Generate generic variables to be use in definitions
foreach(element ${elements})
  string(TOUPPER ${element} element_${upper})
  set(${element}_factory ${element_${upper}}FEMFactory)
  set(${element}_use_factory USE_${element_${upper}}_FEM_FACTORY)
endforeach(element ${elements})

# raviart thomas is different, we therefore overwrite the values with the
# correct ones.
set(rt_factory RTBDMFEMFactory)
set(rt_use_factory USE_RT_BDM_FEM_FACTORY)
set(rt_factory_class_definition
  FEM_FACTORY_FEM_CLASS=Dune::PDELab::RaviartThomasLocalFiniteElementMap)

foreach(geo ${geometries})
  string(TOUPPER ${geo} ${geo}_upper)
  list(APPEND GEOMETRIES ${GEO})
endforeach(geo ${geometries})


set(TESTS)


# Raviart-Thomas
# RT0 and RT1 in 2D for simplex and cube

foreach(geometry ${geometries})
  foreach(order 0 1)
    set(target testrt${order}${geometry}2dfem)
    add_executable(${target} testfemwithgfs.cc)
    set_property(TARGET ${target} APPEND PROPERTY COMPILE_DEFINITIONS
      FEM_DIM=2 FEM_${${geometry}_upper}
      USE_RT_BDM_FEM_FACTORY
      FEM_FACTORY=RTBDMFEMFactory
      FEM_FACTORY_FEM_CLASS=Dune::PDELab::RaviartThomasLocalFiniteElementMap
      FEM_FACTORY_ORDER=${order}
      GRIDSDIR="${CMAKE_SOURCE_DIR}/grids" )
    add_dune_gmp_flags(${target})
    add_dune_alugrid_flags(${target})
    list(APPEND TESTS ${target})
  endforeach(order 0 1)
endforeach(geometry ${geometries})


# RT2 in 2D for cube

set(target testrt2cube2dfem)
add_executable(${target} testfemwithgfs.cc)
set_property(TARGET ${target} APPEND PROPERTY COMPILE_DEFINITIONS
  FEM_DIM=2 FEM_SIMPLEX_CUBE
  USE_RT_BDM_FEM_FACTORY
  FEM_FACTORY=RTBDMFEMFactory
  FEM_FACTORY_FEM_CLASS=Dune::PDELab::RaviartThomasLocalFiniteElementMap
  FEM_FACTORY_ORDER=2
  GRIDSDIR="${CMAKE_SOURCE_DIR}/grids")
add_dune_gmp_flags(${target})
add_dune_alugrid_flags(${target})
list(APPEND TESTS ${target})


# RT 0 and RT1 in 3D for cube

foreach(order 0 1)
  set(target testrt${order}cube3dfem)
  add_executable(${target} testfemwithgfs.cc)
  set_property(TARGET ${target} APPEND PROPERTY COMPILE_DEFINITIONS
    FEM_DIM=3 FEM_${cube_upper}
    USE_RT_BDM_FEM_FACTORY
    FEM_FACTORY=RTBDMFEMFactory
    FEM_FACTORY_FEM_CLASS=Dune::PDELab::RaviartThomasLocalFiniteElementMap
    FEM_FACTORY_ORDER=${order}
    GRIDSDIR="${CMAKE_SOURCE_DIR}/grids" )
  add_dune_gmp_flags(${target})
  add_dune_alugrid_flags(${target})
  list(APPEND TESTS ${target})
endforeach(order 0 1)


# BDM 1 in 2D for simplex and cube

foreach(geometry ${geometries})
  set(target testbdm1${geometry}2dfem)
  add_executable(${target} testfemwithgfs.cc)
  set_property(TARGET ${target} APPEND PROPERTY COMPILE_DEFINITIONS
    FEM_DIM=2 FEM_${${geometry}_upper}
    USE_RT_BDM_FEM_FACTORY
    FEM_FACTORY=RTBDMFEMFactory
    FEM_FACTORY_FEM_CLASS=Dune::PDELab::BrezziDouglasMariniLocalFiniteElementMap
    FEM_FACTORY_ORDER=1
    GRIDSDIR="${CMAKE_SOURCE_DIR}/grids" )
  add_dune_gmp_flags(${target})
  add_dune_alugrid_flags(${target})
  list(APPEND TESTS ${target})
endforeach(geometry ${geometries})


# OPB order 0-3 in 2D and 3D for simplex and cube

foreach(dimension ${dimensions})
  foreach(geometry ${geometries})
    foreach(order ${orders})
      set(target testopb${order}${geometry}${dimension}dfem)
      add_executable(${target} testfemwithgfs.cc)
      set_property(TARGET ${target} APPEND PROPERTY COMPILE_DEFINITIONS
        FEM_DIM=${dimension} FEM_${${geometry}_upper}
        USE_OPB_FEM_FACTORY
        FEM_FACTORY=OPBFEMFactory
        FEM_FACTORY_ORDER=${order}
        GRIDSDIR="${CMAKE_SOURCE_DIR}/grids" )
      add_dune_gmp_flags(${target})
      add_dune_alugrid_flags(${target})
      list(APPEND TESTS ${target})
    endforeach(order 0 1)
  endforeach(geometry ${geometries})
endforeach(dimension ${dimensions})


# Pk order 0-3 in 2D and 3D for simplex

foreach(dimension ${dimensions})
  foreach(order ${orders})
    set(target testp${order}${geometry}${dimension}dfem)
    add_executable(${target} testfemwithgfs.cc)
    set_property(TARGET ${target} APPEND PROPERTY COMPILE_DEFINITIONS
      FEM_DIM=${dimension} FEM_SIMPLEX
      USE_PK_FEM_FACTORY
      FEM_FACTORY=PKFEMFactory
      FEM_FACTORY_ORDER=${order}
      GRIDSDIR="${CMAKE_SOURCE_DIR}/grids" )
    add_dune_gmp_flags(${target})
    add_dune_alugrid_flags(${target})
    list(APPEND TESTS ${target})
  endforeach(order 0 1)
endforeach(dimension ${dimensions})


# Rannacher-Turek in 2D and 3D for cube

foreach(dimension ${dimensions})
  set(target testrannacherturekcube${dimension}dfem)
  add_executable(${target} testfemwithgfs.cc)
  set_property(TARGET ${target} APPEND PROPERTY COMPILE_DEFINITIONS
    FEM_DIM=${dimension} FEM_CUBE
    USE_RANNACHER_TUREK_FEM_FACTORY
    FEM_FACTORY=RannacherTurekFEMFactory
    GRIDSDIR="${CMAKE_SOURCE_DIR}/grids" )
  add_dune_gmp_flags(${target})
  add_dune_alugrid_flags(${target})
  list(APPEND TESTS ${target})
endforeach(dimension ${dimensions})


# We do not want want to build the tests during make all
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${TESTS})
