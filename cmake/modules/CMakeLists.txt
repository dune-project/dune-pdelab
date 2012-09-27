set(modules DunePdelabMacros.cmake FindPackageMultipass.cmake  FindPETSc.cmake  ResolveCompilerPaths.cmake UsePETSc.cmake)
install(FILES ${modules} DESTINATION ${DUNE_INSTALL_MODULEDIR})
