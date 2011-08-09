cmake_minimum_required(VERSION 2.6)

project(exportpred)

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set(exportpred_VERSION_MAJOR 2)
set(exportpred_VERSION_MINOR 0)
set(exportpred_VERSION_PATCH 0)

set(EXPORTPRED_VERSION ${exportpred_VERSION_MAJOR}.${exportpred_VERSION_MINOR}.${exportpred_VERSION_PATCH})
set(VERSION ${EXPORTPRED_VERSION})

option(BUILD_COVERAGE                  "Compile with gcov"                   OFF)
option(BUILD_SHARED_LIBS               "Compile libghmm as shared"           ON)
option(EXPORTPRED_DEBUG                "Compile in debug code"               OFF)

if (MSVC)
  # For MSVC, CMake sets certain flags to defaults we want to override.
  # This replacement code is taken from sample in the CMake Wiki at
  # http://www.cmake.org/Wiki/CMake_FAQ#Dynamic_Replace.
  foreach (flag_var
           CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
           CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO)
    # In hermetic build environments, tests may not have access to MS runtime
    # DLLs, so this replaces /MD (CRT libraries in DLLs) with /MT (static CRT
    # libraries).
    string(REPLACE "/MD" "-MT" ${flag_var} "${${flag_var}}")
    # We prefer more strict warning checking for building Google Test.
    # Replaces /W3 with /W4 in defaults.
    string(REPLACE "/W3" "-W4" ${flag_var} "${${flag_var}}")
  endforeach()
endif()

if (BUILD_COVERAGE)
  set(BUILD_SHARED_LIBS OFF)
  set(CMAKE_BUILD_TYPE DEBUG)
  IF (CMAKE_COMPILER_IS_GNUCC)
    SET(CMAKE_C_FLAGS "-g -O0 -Wall -fprofile-arcs -ftest-coverage")
    SET(CMAKE_CXX_FLAGS "-g -O0 -Wall -fprofile-arcs -ftest-coverage")
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fprofile-arcs -ftest-coverage -lgcov")
  ENDIF()
endif(BUILD_COVERAGE)

if(WIN32)
  set(BUILD_SHARED_LIBS OFF) # until everything is exported
  add_definitions(-D_USE_MATH_DEFINES)
  add_definitions(-DNOMINMAX)
endif(WIN32)

add_subdirectory(external/gtest-1.6.0)

configure_file (
  "${exportpred_SOURCE_DIR}/include/GHMM/cmake-config.h.in"
  "${exportpred_BINARY_DIR}/include/GHMM/config.h"
)

include_directories(${exportpred_BINARY_DIR}/include)

add_subdirectory(lib) 
add_subdirectory(include) 
add_subdirectory(src) 

include(CMakeCPack.cmake)
include(CMakeCTest.cmake)
