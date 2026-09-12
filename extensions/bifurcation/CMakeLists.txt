cmake_minimum_required(VERSION 3.21)
project(nekrs-bifurcation LANGUAGES CXX)

option(BIF_ENABLE_TRILINOS "Build the NOX/LOCA/Belos/Anasazi adapter" OFF)
option(BIF_ENABLE_NEKRS "Build the nekRS state/flow-map adapter" OFF)
option(BIF_BUILD_APP "Build the nekrs-bif reference executable" ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_library(bifurcation STATIC
  src/core.cpp
  src/finite_difference.cpp
  src/event_tracking.cpp
  src/hopf_localization.cpp
  src/checkpoint.cpp)
target_include_directories(bifurcation PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

if(BIF_BUILD_APP AND NOT BIF_ENABLE_NEKRS)
  add_executable(nekrs-bif app/nekrs_bif.cpp)
  target_link_libraries(nekrs-bif PRIVATE bifurcation)
endif()

if(BIF_ENABLE_NEKRS)
  if(BIF_BUILD_APP AND NOT BIF_ENABLE_TRILINOS)
    message(FATAL_ERROR "BIF_ENABLE_NEKRS with BIF_BUILD_APP requires BIF_ENABLE_TRILINOS for LOCA/Belos/Anasazi")
  endif()
  if(NOT TARGET nekrs-lib)
    message(FATAL_ERROR "BIF_ENABLE_NEKRS requires the extension to be added to a nekRS build or NEKRS_LIB_TARGET to be supplied")
  endif()
  add_library(bifurcation_nekrs STATIC
    adapters/nekrs/session.cpp
    adapters/nekrs/api.cpp)
  target_include_directories(bifurcation_nekrs PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/adapters/nekrs
    ${CMAKE_SOURCE_DIR}/src/lib)
  target_link_libraries(bifurcation_nekrs PUBLIC bifurcation nekrs-lib)
  if(BIF_BUILD_APP)
    add_executable(nekrs-bif app/nekrs_bif_nekrs.cpp)
    target_include_directories(nekrs-bif PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR}/adapters/nekrs
      ${CMAKE_CURRENT_SOURCE_DIR}/adapters/trilinos
      ${CMAKE_SOURCE_DIR}/src/lib
      ${CMAKE_SOURCE_DIR}/src/platform
      ${CMAKE_SOURCE_DIR}/src/platform/utils)
    target_link_libraries(nekrs-bif PRIVATE bifurcation_nekrs bifurcation_trilinos)
  endif()
endif()

if(BIF_ENABLE_TRILINOS)
  # The development tree supplied with nekRS exports package configs from
  # <build>/cmake_packages/NOX; an installed tree exports the same NOX config
  # below its prefix.  NOX::all_libs carries the verified Teuchos/Tpetra,
  # Belos, Anasazi and LOCA closure.
  find_package(NOX CONFIG REQUIRED
    PATHS "${Trilinos_DIR}/cmake_packages/NOX" "${Trilinos_DIR}/NOX"
    NO_DEFAULT_PATH)
  add_library(bifurcation_trilinos STATIC
    adapters/trilinos/source_api_probe.cpp
    adapters/trilinos/nox_group.cpp
    adapters/trilinos/loca_runner.cpp
    adapters/trilinos/linear_solver.cpp
    adapters/trilinos/map_eigensolver.cpp)
  target_include_directories(bifurcation_trilinos PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/adapters/trilinos)
  target_link_libraries(bifurcation_trilinos PUBLIC bifurcation NOX::all_libs)
  # Match the standard recorded by the supplied development build; a mismatched
  # MPI/compiler/standard ABI is rejected by the build recipe rather than hidden.
  set_property(TARGET bifurcation_trilinos PROPERTY CXX_STANDARD 20)
  set_property(TARGET bifurcation_trilinos PROPERTY CXX_STANDARD_REQUIRED ON)
  target_compile_definitions(bifurcation_trilinos PUBLIC BIF_WITH_TRILINOS=1)
endif()

install(TARGETS bifurcation
        ARCHIVE DESTINATION lib
        INCLUDES DESTINATION include)
install(DIRECTORY include/ DESTINATION include)
if(BIF_BUILD_APP AND TARGET nekrs-bif)
  install(TARGETS nekrs-bif RUNTIME DESTINATION bin)
endif()
