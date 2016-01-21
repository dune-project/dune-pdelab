set(modules
  CorrectWindowsPaths.cmake
  DunePdelabMacros.cmake
  DunePdelabTestMacros.cmake
  FindCXXUnorderedContainers.cmake
  FindPackageMultipass.cmake
  FindPETSc.cmake
  PdelabCheckCXX14.cmake
  ResolveCompilerPaths.cmake
  UseEigen.cmake
  UsePETSc.cmake)
install(FILES ${modules} DESTINATION ${DUNE_INSTALL_MODULEDIR})
