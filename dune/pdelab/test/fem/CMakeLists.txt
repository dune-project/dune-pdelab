set(TESTS testrt0simplex2dfem testrt1simplex2dfem testrt0cube2dfem testrt1cube2dfem testrt2cube2dfem testrt0cube3dfem testbdm1simplex2dfem testbdm1cube2dfem testopb0cube2dfem testopb1cube2dfem testopb2cube2dfem testopb3cube2dfem  testopb0simplex2dfem testopb1simplex2dfem testopb2simplex2dfem

set(orders 0 1 2 3)
set(dimensions 2 3)
set(geometries simplex cube)

set(elements rt opb pk)

# Generate generic variables to be use in definitions
foreach(element ${elements})
  string(TOUPPER ${element} element_${upper})
  set(${element}_factory ${element_${upper}}FEMFACTORY)
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

# Tests using rtbmfem.hh
foreach(element ${elements})
  foreach(dimension ${dimensions})
    foreach(geometry ${geometries})
      foreach(order 0 1)
        set(target testrt${order}${geometry}${dimension}dfem)
        add_executable(${target} testfemwithgfs.cc)
        set_property(${target} APPEND PROPERTY COMPILE_DEFINITIONS
          FEM_DIM=${dimension} FEM_${${geometry}_upper}
          ${${element}_use_factory}
          FEM_FACTORY=${${element}_factory}
          ${${element}_factory_class_definition}
          FEM_FACTORY_ORDER=${order})
        add_dune_gmp_flags(${target})
        list(APPEND TESTS ${target})
      endforeach(order 0 1)
    endforeach(geometry ${geometries})
  endforeach(dimension ${dimensions})
endforeach(element ${elements})

# We do not want want to build the tests during make all
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${TEST})
