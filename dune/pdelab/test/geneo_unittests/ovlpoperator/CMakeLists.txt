dune_add_test(NAME test-geneo-ovlpoperator
              SOURCES test-geneo-ovlpoperator.cc
              MPI_RANKS 2
              TIMEOUT 300)

dune_symlink_to_source_files(FILES test-geneo-ovlpoperator.ini
                                   referencematrix_0
                                   referencematrix_1)
