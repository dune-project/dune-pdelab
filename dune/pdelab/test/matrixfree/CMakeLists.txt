#=========================
# Linear Matrix-Free Tests
#=========================

add_executable(matrix_free_linear_matrix_based matrix_free_linear.cc)
target_compile_definitions(matrix_free_linear_matrix_based PRIVATE MATRIX_BASED)

add_executable(matrix_free_linear_partial_matrix_free matrix_free_linear.cc)
target_compile_definitions(matrix_free_linear_partial_matrix_free PRIVATE PARTIAL_MATRIX_FREE)

add_executable(matrix_free_linear_fully_matrix_free matrix_free_linear.cc)
target_compile_definitions(matrix_free_linear_fully_matrix_free PRIVATE FULLY_MATRIX_FREE)

add_executable(matrix_free_linear_matrix_based_sor matrix_free_linear.cc)
target_compile_definitions(matrix_free_linear_matrix_based_sor PRIVATE MATRIX_BASED_SOR)

add_executable(matrix_free_linear_matrix_free_sor matrix_free_linear.cc)
target_compile_definitions(matrix_free_linear_matrix_free_sor PRIVATE MATRIX_FREE_SOR)

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

dune_symlink_to_source_files(FILES matrix_free_nonlinear.ini)