set(modules
  DunePdelabMacros.cmake
  FindCXXUnorderedContainers.cmake
  FindPackageMultipass.cmake
  FindPETSc.cmake
  ResolveCompilerPaths.cmake
  UsePETSc.cmake)
install(FILES ${modules} DESTINATION ${DUNE_INSTALL_MODULEDIR})
