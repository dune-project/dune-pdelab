dune_add_test(NAME test-geneo-partitionofunity
              SOURCES test-geneo-partitionofunity.cc
              MPI_RANKS 2
              TIMEOUT 300)

dune_symlink_to_source_files(FILES test-geneo-partitionofunity.ini)
