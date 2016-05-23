set(modules
  AddVCLFlags.cmake
  CorrectWindowsPaths.cmake
  DunePdelabMacros.cmake
  DunePdelabTestMacros.cmake
  FindPackageMultipass.cmake
  FindPETSc.cmake
  FindVCL.cmka
  ResolveCompilerPaths.cmake
  UseEigen.cmake
  UsePETSc.cmake)
install(FILES ${modules} DESTINATION ${DUNE_INSTALL_MODULEDIR})
