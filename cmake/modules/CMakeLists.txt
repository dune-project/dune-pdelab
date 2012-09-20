set(modules DunePdelabMacros.cmake FindPackageMultipass.cmake  FindPETSc.cmake  ResolveCompilerPaths.cmake)
install(FILES ${modules} DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/cmake/modules)
