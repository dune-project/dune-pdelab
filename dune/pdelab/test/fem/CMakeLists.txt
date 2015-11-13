set(orders 0 1 2 3)
set(dimensions 2 3)
set(geometries simplex cube)
# store upper case strings for the geometries, too.
foreach(geo ${geometries})
  string(TOUPPER ${geo} ${geo}_upper)
endforeach(geo ${geometries})

# Raviart-Thomas
# RT0 and RT1 in 2D for simplex and cube
foreach(geometry ${geometries})
  foreach(order 0 1)
    dune_add_test(NAME testrt${order}${geometry}2dfem
                  SOURCES testfemwithgfs.cc
                  COMPILE_DEFINITIONS FEM_DIM=2
                      FEM_${${geometry}_upper}
                      USE_RT_BDM_FEM_FACTORY
                      FEM_FACTORY=RTBDMFEMFactory
                      FEM_FACTORY_FEM_CLASS=Dune::PDELab::RaviartThomasLocalFiniteElementMap
                      FEM_FACTORY_ORDER=${order}
                      GRIDSDIR="${CMAKE_CURRENT_SOURCE_DIR}/grids"
    )
  endforeach()
endforeach()

# RT2 in 2D for cube
dune_add_test(NAME testrt2cube2dfem
              SOURCES testfemwithgfs.cc
              COMPILE_DEFINITIONS
                  FEM_DIM=2
                  FEM_CUBE
                  USE_RT_BDM_FEM_FACTORY
                  FEM_FACTORY=RTBDMFEMFactory
                  FEM_FACTORY_FEM_CLASS=Dune::PDELab::RaviartThomasLocalFiniteElementMap
                  FEM_FACTORY_ORDER=2
                  GRIDSDIR="${CMAKE_CURRENT_SOURCE_DIR}/grids"
)

# RT 0 and RT1 in 3D for cube
foreach(order 0 1)
  dune_add_test(NAME testrt${order}cube3dfem
                SOURCES testfemwithgfs.cc
                COMPILE_DEFINITIONS
                    FEM_DIM=3
                    FEM_CUBE
                    USE_RT_BDM_FEM_FACTORY
                    FEM_FACTORY=RTBDMFEMFactory
                    FEM_FACTORY_FEM_CLASS=Dune::PDELab::RaviartThomasLocalFiniteElementMap
                    FEM_FACTORY_ORDER=${order}
                    GRIDSDIR="${CMAKE_CURRENT_SOURCE_DIR}/grids"
  )
endforeach()

# BDM 1 in 2D for simplex and cube
foreach(geometry ${geometries})
  dune_add_test(NAME testbdm1${geometry}2dfem
                SOURCES testfemwithgfs.cc
                COMPILE_DEFINITIONS
                      FEM_DIM=2
                      FEM_${${geometry}_upper}
                      USE_RT_BDM_FEM_FACTORY
                      FEM_FACTORY=RTBDMFEMFactory
                      FEM_FACTORY_FEM_CLASS=Dune::PDELab::BrezziDouglasMariniLocalFiniteElementMap
                      FEM_FACTORY_ORDER=1
                      GRIDSDIR="${CMAKE_CURRENT_SOURCE_DIR}/grids"
   )
endforeach()

# OPB order 0-3 in 2D and 3D for simplex and cube
foreach(dimension ${dimensions})
  foreach(geometry ${geometries})
    foreach(order ${orders})
      dune_add_test(NAME testopb${order}${geometry}${dimension}dfem
                    SOURCES testfemwithgfs.cc
                    COMPILE_DEFINITIONS
                        FEM_DIM=${dimension}
                        FEM_${${geometry}_upper}
                        USE_OPB_FEM_FACTORY
                        FEM_FACTORY=OPBFEMFactory
                        FEM_FACTORY_ORDER=${order}
                        GRIDSDIR="${CMAKE_CURRENT_SOURCE_DIR}/grids")
      endforeach()
  endforeach()
endforeach()

# Pk order 0-3 in 2D and 3D for simplex
foreach(dimension ${dimensions})
  foreach(order ${orders})
    dune_add_test(NAME testp${order}${dimension}dfem
                  SOURCES testfemwithgfs.cc
                  COMPILE_DEFINITIONS
                      FEM_DIM=${dimension}
                      FEM_${${geometry}_upper}
                      USE_PK_FEM_FACTORY
                      FEM_FACTORY=PKFEMFactory
                      FEM_FACTORY_ORDER=${order}
                      GRIDSDIR="${CMAKE_CURRENT_SOURCE_DIR}/grids"
    )
  endforeach()
endforeach()

# Rannacher-Turek in 2D and 3D for cube
foreach(dimension ${dimensions})
  dune_add_test(NAME testrannacherturekcube${dimension}dfem
                SOURCES testfemwithgfs.cc
                COMPILE_DEFINITIONS
                    FEM_DIM=${dimension}
                    FEM_CUBE
                    USE_RANNACHER_TUREK_FEM_FACTORY
                    FEM_FACTORY=RannacherTurekFEMFactory
                    GRIDSDIR="${CMAKE_CURRENT_SOURCE_DIR}/grids"
  )
endforeach()
