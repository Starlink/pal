# cspell:disable
cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

set(VER_MAJ 0) # major
set(VER_MIN 9) # minor
set(VER_PAT 10) # patch
set(VER_SS ${VER_MAJ}.${VER_MIN}.${VER_PAT})
message("Version: ${VER_SS}")

project(cx_pal VERSION ${VER_SS})

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
add_definitions(-DPACKAGE_VERSION="${VER_SS}")
add_definitions(-DPACKAGE_VERSION_INTEGER=${VER_MAJ})

set(CMAKE_BUILD_TYPE Release)

find_package(PkgConfig)
pkg_check_modules(Erfa REQUIRED erfa>=1.7)

include_directories(${PROJECT_SOURCE_DIR}/)

file(GLOB sources ${PROJECT_SOURCE_DIR}/*.c)
file(GLOB headers ${PROJECT_SOURCE_DIR}/*.h)

set(LINK_LIBS erfa)

add_library(pal SHARED ${sources} ${headers})
target_link_libraries(pal ${LINK_LIBS})
set_target_properties(pal PROPERTIES VERSION ${VER_SS} SOVERSION
                                                          ${VER_MAJ})
install(FILES ${headers} DESTINATION include/star)
install(
  TARGETS pal
  RUNTIME DESTINATION bin/
  LIBRARY DESTINATION lib/
  ARCHIVE DESTINATION lib/)
