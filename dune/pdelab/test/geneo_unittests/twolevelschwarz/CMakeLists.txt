dune_add_test(NAME test-geneo-twolevelschwarz
              SOURCES test-geneo-twolevelschwarz.cc
              MPI_RANKS 2
              TIMEOUT 300)

dune_symlink_to_source_files(FILES test-geneo-twolevelschwarz.ini)
