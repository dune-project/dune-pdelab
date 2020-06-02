#=========================
# Linear Matrix-Free Tests
#=========================

dune_add_test(NAME matrix_free_linear_partial_matrix_free_fastdg
              SOURCES matrix_free_linear.cc
              COMPILE_DEFINITIONS PARTIAL_MATRIX_FREE FASTDG)

dune_add_test(NAME matrix_free_linear_partial_matrix_free_nonfast
              SOURCES matrix_free_linear.cc
              COMPILE_DEFINITIONS PARTIAL_MATRIX_FREE)

dune_add_test(NAME matrix_free_linear_fully_matrix_free_fastdg
              SOURCES matrix_free_linear.cc
              COMPILE_DEFINITIONS FULLY_MATRIX_FREE FASTDG)

dune_add_test(NAME matrix_free_linear_fully_matrix_free_nonfast
              SOURCES matrix_free_linear.cc
              COMPILE_DEFINITIONS FULLY_MATRIX_FREE)

# SOR only works with FASTDG
dune_add_test(NAME matrix_free_linear_matrix_free_sor_fastdg
              SOURCES matrix_free_linear.cc
              COMPILE_DEFINITIONS SOR FASTDG)

dune_symlink_to_source_files(FILES matrix_free_linear.ini)

#============================
# Nonlinear Matrix-Free Tests
#============================

add_executable(matrix_free_nonlinear_matrix_based matrix_free_nonlinear.cc)
target_compile_definitions(matrix_free_nonlinear_matrix_based PRIVATE MATRIX_BASED)

add_executable(matrix_free_nonlinear_partial_matrix_free matrix_free_nonlinear.cc)
target_compile_definitions(matrix_free_nonlinear_partial_matrix_free PRIVATE PARTIAL_MATRIX_FREE)

add_executable(matrix_free_nonlinear_fully_matrix_free matrix_free_nonlinear.cc)
target_compile_definitions(matrix_free_nonlinear_fully_matrix_free PRIVATE FULLY_MATRIX_FREE)

add_executable(matrix_free_nonlinear_matrix_based_sor matrix_free_nonlinear.cc)
target_compile_definitions(matrix_free_nonlinear_matrix_based_sor PRIVATE MATRIX_BASED_SOR)

add_executable(matrix_free_nonlinear_matrix_free_sor matrix_free_nonlinear.cc)
target_compile_definitions(matrix_free_nonlinear_matrix_free_sor PRIVATE MATRIX_FREE_SOR)

dune_symlink_to_source_files(FILES matrix_free_nonlinear.ini)
