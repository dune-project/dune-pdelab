# creates and export a Perfetto::SDK target

if (NOT TARGET Perfetto::SDK)
  include(FetchContent)
  message("-- Declaring Perfetto")
  FetchContent_Declare(
    perfetto
    GIT_REPOSITORY https://android.googlesource.com/platform/external/perfetto
    GIT_TAG        v32.1
  )

  # configure perfetto targets
  if(NOT perfetto_POPULATED)
    message("-- Populating Perfetto")
    FetchContent_Populate(perfetto)

    dune_add_library(perfetto-sdk STATIC
      SOURCES "${perfetto_SOURCE_DIR}/sdk/perfetto.cc"
      EXPORT_NAME Perfetto::SDK)

    target_link_libraries(perfetto-sdk PUBLIC $<$<BOOL:${Threads_FOUND}>:Threads::Threads>)
    target_compile_definitions(perfetto-sdk INTERFACE HAVE_PERFETTO)
    install(DIRECTORY ${perfetto_SOURCE_DIR}/sdk DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/external/perfetto-sdk)
    target_include_directories(perfetto-sdk
      PUBLIC
        $<BUILD_INTERFACE:${perfetto_SOURCE_DIR}/sdk>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/dune/external/perfetto-sdk/sdk>)
  endif()
endif()
