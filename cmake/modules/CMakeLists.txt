set(modules
  DunePdelabMacros.cmake
  FindCXXDecltype.cmake
  FindCXXInitializerList.cmake
  FindCXXUnorderedContainers.cmake
  FindGCCTypeof.cmake
  FindPackageMultipass.cmake
  FindPETSc.cmake
  ResolveCompilerPaths.cmake
  UsePETSc.cmake)
install(FILES ${modules} DESTINATION ${DUNE_INSTALL_MODULEDIR})
